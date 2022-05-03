/**CALCSPHERCONVHESSIAN Calculate the Hessian (a matrix of second partial
*            derivatives) for a monostatic or bistatic range and a
*            spherical direction measurement in 3D, ignoring atmospheric
*            effects, with respect to Cartesian position. This type of
*            Jacobian is useful when performing tracking using Cartesian-
*            converted measurements where the clutter density is specified
*            in the measurement coordinate system, not the converted
*            measurement coordinate system.
*
*INPUTS: zSpher A 3X1 point in range and azimuth and angle in the format
*           [range;azimuth;angle], where the angles are given in radians.
* systemType An optional parameter specifying the axes from which the
*           angles are measured in radians. Possible values are
*           0 (The default if omitted) Azimuth is measured 
*             counterclockwise from the x-axis in the x-y plane. Elevation
*             is measured up from the x-y plane (towards the z-axis). This
*             is consistent with common spherical coordinate systems for
*             specifying longitude (azimuth) and geocentric latitude
*             (elevation).
*           1 Azimuth is measured counterclockwise from the z-axis in the
*             z-x plane. Elevation is measured up from the z-x plane
*             (towards the y-axis). This is consistent with some spherical
*             coordinate systems that use the z axis as the boresight
*             direction of the radar.
*           2 This is the same as 0 except instead of being given
*             elevation, one desires the angle away from the z-axis, which
*             is (pi/2-elevation).
*           3 This is the same as 0 except azimuth is measured clockwise
*             from the y-axis in the x-y plane instead of counterclockwise
*             from the x-axis. This coordinate system often arises when
*             given "bearings" in a local East-North-Up coordinate system,
*             where the bearing directions are measured East of North.
* useHalfRange An optional boolean value specifying whether the bistatic
*           (round-trip) range value has been divided by two. This normally
*           comes up when operating in monostatic mode (the most common
*           type of spherical coordinate system), so that the range
*           reported is a one-way range (or just half a bistatic range).
*           The default if this parameter is not provided is false if lTx
*           is provided and true if it is are omitted (monostatic). If no
*           range values are provided, an empty matrix can be passed.
*       lTx The 3X1 [x;y;z] location vector of the transmitter in global
*           Cartesian coordinates. If this parameter is omitted or an
*           empty matrix is passed, then the transmitter is assumed to be
*           at the origin.
*       lRx The 3X1 [x;y;z] location vector of the receiver in Cartesian
*           coordinates. If this parameter is omitted or an empty matrix
*           is passed, then the receiver is assumed to be at the origin.
*         M A 3X3 rotation matrix to go from the alignment of the global
*           coordinate system to that at the receiver. The z-axis of the
*           local coordinate system of the receiver is the pointing
*           direction of the receiver. If omitted or an empty matrix is
*           passed, then it is assumed that the local coordinate system is
*           aligned with the global and M=eye(3) --the identity matrix is
*           used. 
*
*OUTPUTS: H The 3X3X3 Hessian matrix, where H(:,:,1) is the Hessian for the
*           range component, H(:,:,2) is the Hessian for the azimuth
*           component, and H(:,:,3) is the Hessian for the elevation
*           component. The ordering of the derivatives in each matrix is:
*                  [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
*                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
*                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
*          note that each matrix is symmetric (i.e.
*                   d^2/(dydx)=d^2/(dxdy) ).
*
*The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*H=calcSpherConvHessian(zSpher,systemType,useHalfRange,lTx,lRx,M)
*
*July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    const double *zSpher, *lRx, *lTx, *M;
    const double zeros[]={0,0,0};
    const double identMat[]={1,0,0,0,1,0,0,0,1};
    size_t systemType;
    bool useHalfRange;
    mxArray *HMATLAB;
    double *H;

    if(nrhs<1||nrhs>6) {
        mexErrMsgTxt("Wrong number of inputs.");
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    if(mxGetM(prhs[0])!=3||mxIsEmpty(prhs[0])||mxGetN(prhs[0])>1) {
        mexErrMsgTxt("x has the wrong dimensionality.");
    }

    checkRealDoubleArray(prhs[0]);
    zSpher=mxGetDoubles(prhs[0]);
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        systemType=getSizeTFromMatlab(prhs[1]);
    } else {
        systemType=0;
    }
    
    if(systemType!=0&&systemType!=1&&systemType!=2&&systemType!=3) {
        mexErrMsgTxt("Unknown system type specified.");
    }
    
    if(nrhs<3||mxIsEmpty(prhs[2])) {
        if(nrhs<4||mxIsEmpty(prhs[3])) {//If lTx is omitted
            useHalfRange=true;
        } else {
            useHalfRange=false; 
        }
    } else {
        useHalfRange=getBoolFromMatlab(prhs[2]);
    }
    
    //Allocate space for the return values
    {
        mwSize dims[3];
  
        dims[0]=3;
        dims[1]=3;
        dims[2]=3;
            
        HMATLAB=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        H=mxGetDoubles(HMATLAB);
    }
    
    if(nrhs>3) {
        if(nrhs>3&&!mxIsEmpty(prhs[3])) {
            if(mxGetM(prhs[3])!=3||mxGetN(prhs[3])!=1) {
                mexErrMsgTxt("lTx has the wrong dimensionality.");
            }
            checkRealDoubleArray(prhs[3]);

            lTx=mxGetDoubles(prhs[3]);
        } else {
            lTx=zeros;
        }

        if(nrhs>4&&!mxIsEmpty(prhs[4])) {
            if(mxGetM(prhs[4])!=3||mxGetN(prhs[4])!=1) {
                mexErrMsgTxt("lRx has the wrong dimensionality.");
            }
            checkRealDoubleArray(prhs[4]);

            lRx=mxGetDoubles(prhs[4]);
        } else {
            lRx=zeros;
        }

        if(nrhs>5&&!mxIsEmpty(prhs[5])) {
            if(mxGetM(prhs[5])!=3||mxGetN(prhs[5])!=3) {
                mexErrMsgTxt("M has the wrong dimensionality.");
            }

            checkRealDoubleArray(prhs[5]);
            M=mxGetDoubles(prhs[5]);
        } else {
            M=identMat;
        }
    
        calcSpherConvHessianGenCPP(H,zSpher,systemType,useHalfRange,lTx,lRx,M);
    } else {
        calcSpherConvHessianCPP(H,zSpher,systemType,useHalfRange);
    }

    plhs[0]=HMATLAB;
    
}

/*LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
