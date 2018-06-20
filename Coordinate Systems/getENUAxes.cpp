/**GETENUAXES A C++ version of a function to find the East-North-Up unit
*             vectors at the given point as well as the magnitudes of the
*             derivatives of a position vector with respect to latitude,
*             longitude and height.
*
*INPUTS:   plhPoint The point at which the axes are to be found given in
%                   terms of [latitude;longitude] (with an assumed
%                   ellipsoidal height of 0) or [latitude;longitude;height]
%                   with the geodetic latitude and longitude in radians and
%                   the height in meters. The latitude should be between
%                   -pi/2 and pi/2. The height does not change the unit
%                   direction vectors u but it does change c.
*      justVertical An optional parameter. If this is given and is
*                   true, then u and c only for the Up direction will be
*                   returned. The default is false.
*           a       The semi-major axis of the reference ellipsoid. If
*                   this argument is omitted, the value in
*                   Constants.WGS84SemiMajorAxis is used.
*           f       The flattening factor of the reference ellipsoid. If
*                   this argument is omitted, the value in
*                   Constants.WGS84Flattening is used.
*
*OUTPUTS:   u       u(:,1), u(:,2) and u(:,3) are respectively the East,
*                   North and Up unit vectors.
*           c       c(1), c(2) and c(3) are the respective magnitudes of
*                   the derivative of the Cartesian position with respect
*                   to latitude, longitude and height.
*
*If only the Up vector is desired, include the fourth parameter set to
*true. Further details on the use of the function can be found in the
*Matlab implemntation of getENUAxes.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[u,c]=getENUAxes(plhPoint,justVertical,a,f);
*
*February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double a, f, plhPoint[3];
    bool justVertical;
    mxArray *uMATLAB, *cMATLAB;
    double *u, *c;
    
    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    //Get the latitude, longitude, height point.
    checkRealDoubleArray(prhs[0]);
    {
    size_t mDim=mxGetM(prhs[0]);
    size_t nDim=mxGetN(prhs[0]);
    double *thePoint=reinterpret_cast<double*>(mxGetData(prhs[0]));
    //If a point with ellipsoidal height is provided.
    if((mDim==3&&nDim==1)||(mDim==1&&nDim==3)) {
        plhPoint[0]=thePoint[0];
        plhPoint[1]=thePoint[1];
        plhPoint[2]=thePoint[2];
    } else {
        //If a point without ellipsoidal height is provided, set the height
        //to zero.
        if((mDim==2&&nDim==1)||(mDim==1&&nDim==2)) {
            plhPoint[0]=thePoint[0];
            plhPoint[1]=thePoint[1];
            plhPoint[2]=0;
        }else {
            mexErrMsgTxt("The point has the wrong dimensionality.");
            return;
        }
    }
    }

    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        justVertical=getBoolFromMatlab(prhs[1]);
    }
    else {
        justVertical=false;
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        a=getDoubleFromMatlab(prhs[2]);
    } else {//Load the default value if none is supplied.
        a=getScalarMatlabClassConst("Constants","WGS84SemiMajorAxis");
    }
    
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        f=getDoubleFromMatlab(prhs[3]);
    } else {//Load the default value if none is supplied.
        f=getScalarMatlabClassConst("Constants","WGS84Flattening");
    }
        
    //Allocate the return variables.
    if(justVertical==false) {
        uMATLAB=mxCreateDoubleMatrix(3, 3,mxREAL);
        cMATLAB=mxCreateDoubleMatrix(3, 1,mxREAL);
    } else {
        uMATLAB=mxCreateDoubleMatrix(3, 1,mxREAL);
        cMATLAB=mxCreateDoubleMatrix(1, 1,mxREAL);
    }
    u=reinterpret_cast<double*>(mxGetData(uMATLAB));
    c=reinterpret_cast<double*>(mxGetData(cMATLAB));
    
    //Compute the values to return.
    getENUAxesCPP(u,c,plhPoint,justVertical,a,f);
    
    //Set the return values.
    plhs[0]=uMATLAB;
    if(nlhs>1) {
        plhs[1]=cMATLAB;
    } else {
        mxDestroyArray(cMATLAB);
    }
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
