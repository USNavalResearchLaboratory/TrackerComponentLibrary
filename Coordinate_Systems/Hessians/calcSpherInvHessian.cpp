/**CALCSPHERINVHESSIAN Determine the Hessian matrix (a matrix of second
*           partial derivatives) of a 3D Cartesian point with respect to
*           monostatic spherical range, azimuth, and elevation components.
*
*INPUTS: z The 3XN position vectors in the global spherical coordinate
*          system, each with [range;Az;El] components.
* systemType An optional parameter specifying the axis from which the
*          angles are measured in radians. Possible values are
*          0 (The default if omitted) Azimuth is measured 
*            counterclockwise from the x-axis in the x-y plane. Elevation
*            is measured up from the x-y plane (towards the z-axis). This
*            is consistent with common spherical coordinate systems for
*            specifying longitude (azimuth) and geocentric latitude
*            (elevation).
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
*
*OUTPUTS: HTotal The 3X3X3XN set of Hessian matrices, where H(:,:,1,i) is
*           the Hessian for the x component for the ith point in z,
*           H(:,:,2) is the Hessian for the y component, and H(:,:,3) is
*           the Hessian for the z component. The ordering of the
*           derivatives in each matrix is:
*                  [d^2/(drdr),  d^2/(drdAz),  d^2/(drdEl);
*                   d^2/(dAzdr), d^2/(dAzdAz), d^2/(dAzdEl);
*                   d^2/(dEldr), d^2/(dEldAz), d^2/(dEldE;)];
*           note that each matrix is symmetric (i.e. 
*           d^2/(dAzdr)=d^2/(drdAz) ).
*
*More details are given in the native Matlab implementation.
*
*The algorithm can be compiled for use in Matlab using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*HTotal=calcSpherInvHessian(z,systemType);
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
    const double *z;
    size_t N, systemType;
    mxArray *HMATLAB;
    double *H;
    size_t curPoint;
    
    if(nrhs<1||nrhs>2) {
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
    z=mxGetDoubles(prhs[0]);
    
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
        dims[2]=3;
        dims[3]=N;
            
        HMATLAB=mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
        H=mxGetDoubles(HMATLAB);
    }
    
    for(curPoint=0;curPoint<N;curPoint++) {
        calcSpherInvHessianCPP(H,z,systemType);

        z+=3;
        H+=27;
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
