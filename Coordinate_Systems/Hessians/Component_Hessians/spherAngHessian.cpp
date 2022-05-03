/**SPHERANGHESSIAN Determine the Hessian matrix (a matrix of second partial
*          derivatives) of a 3D spherical azimuth and elevation measurement
*          with respect to 3D position. Relativity and atmospheric effects
*          are not taken into account.
*
*INPUTS: xG The 3XN position vectors in the global coordinate system, each
*          with [x;y;z] components.
* systemType An optional parameter specifying the axes from which the
*          angles are measured in radians. Possible vaues are
*          0 (The default if omitted) Azimuth is measured counterclockwise
*            from the x-axis in the x-y plane. Elevation is measured up
*            from the x-y plane (towards the z-axis). This is consistent
*            with common spherical coordinate systems for specifying
*            longitude (azimuth) and geocentric latitude (elevation).
*          1 Azimuth is measured counterclockwise from the z-axis in the
*            z-x plane. Elevation is measured up from the z-x plane
*            (towards the y-axis). This is consistent with some spherical
*            coordinate systems that use the z axis as the boresight
*            direction of the radar.
*          2 This is the same as 0 except instead of being given
*            elevation, one desires the angle away from the z-axis, which
*            is (pi/2-elevation).
*          3 This is the same as 0 except azimuth is measured clockwise
*            from the y-axis in the x-y plane instead of counterclockwise
*            from the x-axis. This coordinate system often arises when
*            given "bearings" in a local East-North-Up coordinate system,
*            where the bearing directions are measured East of North.
*      lRx The 3X1 position vector of the receiver. If omitted, the
*          receiver is placed at the origin. All vectors in x are assumed
*          to be from the same receiver.
*        M A 3X3 rotation matrix from the global Coordinate system to the
*          orientation of the coordinate system at the receiver. If
*          omitted, it is assumed to be the identity matrix. All vectors in
*          x are assumed to have the same rotation matrix.
*
*OUTPUTS: HTotal A 3X3X2XN matrix such that HTotal(:,:,1,i) is the Hessian
*          matrix with respect to the azimuth component and HTotal(:,:,2,i)
*          is the Hessian matrix with respect to the elevation component.
*          The elements in the matrices for each component/ point are
*          ordered [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
*                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
*                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
*          note that each matrix is symmetric (i.e.
*           d^2/(dydx)=d^2/(dxdy) ).
*
*More details are given in the native Matlab implementation.
*
*The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*HTotal=spherAngHessian(xG,systemType,lRx,M);
*
*June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    const double *xG, *lRx, *M;
    const double zeros[]={0,0,0};
    const double identMat[]={1,0,0,0,1,0,0,0,1};
    mxArray *HMATLAB;
    double *H;
    size_t curPoint;
    
    size_t N, systemType;

    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Wrong number of inputs.");
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    if(mxGetM(prhs[0])!=3||mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("xG has the wrong dimensionality.");
    }
    N=mxGetN(prhs[0]);
    checkRealDoubleArray(prhs[0]);
    xG=mxGetDoubles(prhs[0]);
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        systemType=getSizeTFromMatlab(prhs[1]);
    } else {
        systemType=0;
    }
    
    if(systemType!=0&&systemType!=1&&systemType!=2&&systemType!=3) {
        mexErrMsgTxt("Unknown system type specified.");
    }

    //Allocate space for the return values
    {
        mwSize dims[4];
  
        dims[0]=3;
        dims[1]=3;
        dims[2]=2;
        dims[3]=N;
            
        HMATLAB=mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
        H=mxGetDoubles(HMATLAB);
    }
    
    if(nrhs>2) {
        if(!mxIsEmpty(prhs[2])) {
            if(mxGetM(prhs[2])!=3||mxGetN(prhs[2])!=1) {
                mexErrMsgTxt("lRx has the wrong dimensionality.");
            }
            checkRealDoubleArray(prhs[2]);

            lRx=mxGetDoubles(prhs[2]);
        } else {
            lRx=zeros;
        }

        if(nrhs>3&&!mxIsEmpty(prhs[3])) {
            if(mxGetM(prhs[3])!=3||mxGetN(prhs[3])!=3) {
                mexErrMsgTxt("M has the wrong dimensionality.");
            }

            checkRealDoubleArray(prhs[3]);
            M=mxGetDoubles(prhs[3]);
        } else {
            M=identMat;
        }
        
        for(curPoint=0;curPoint<N;curPoint++) {
            spherAngHessianGenCPP(H,xG,systemType,lRx,M);

            xG+=3;
            H+=18;
        }
    } else {
        for(curPoint=0;curPoint<N;curPoint++) {
            spherAngHessianCPP(H,xG,systemType);

            xG+=3;
            H+=18;
        }
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
