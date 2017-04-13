/**SPHERANGGRADIENTCPP Determine the gradient of a 3D spherical azimuth and
*          elevation measurement with respect to 3D position. Higher order
*          gradient terms are not provided and are zero. Relativity and
*          atmospheric effects are not taken into account.
*
*INPUTS: x The 3X1 target position vector in the global coordinate system
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
*      lRx The 3X1  position vector of the receiver. If omitted, the
*          receiver is placed at the origin.
*        M A 3X3 rotation matrix from the global Coordinate system to the
*          orientation of the coordinate system at the receiver. If
*          omitted, it is assumed to be the identity matrix.
*
*OUTPUTS: J A 2X3 Jacobian matrix where the rows are [azimuth;elevation] in
*           that order and the columns take the derivative of the rows
*           component with respect to [x,y,z] in that order.
*
*The derivatives can be computed in a straightforward manner from
*the basic relation between spherical and Cartesian coordinates, which is
*given in Chapter 14.4.4.1 of [1], among other sources.
*
*Note that singularities exist at the poles; that is when the elevation is
*+/-(pi/2).
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*J=spherAngGradient(xG,systemType,lRx,M);
*
*REFERENCES:
*[1] R. L. Duncombe, "Computational techniques," in Explanatory Supplement
*    to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
*    Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013.
*
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double *xG, *lRx, *M;
    double lRxLocal[3], MLocal[9];
    int systemType;
    mxArray *retMat;
    double *retData;
    
    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxGetM(prhs[0])!=3||mxGetN(prhs[0])!=1) {
       mexErrMsgTxt("The position has the wrong dimensionality.");
       return;
    }

    checkRealDoubleArray(prhs[0]);
    xG=reinterpret_cast<double*>(mxGetData(prhs[0]));
    
    if(nrhs<2||mxIsEmpty(prhs[1])) {
        systemType=0;
    } else {
        systemType=getIntFromMatlab(prhs[1]);
    }
    
    if(nrhs<3||mxIsEmpty(prhs[2])) {
        lRx=lRxLocal;
        lRx[0]=0;
        lRx[1]=0;
        lRx[2]=0;
    } else {
        if(mxGetM(prhs[2])!=3||mxGetN(prhs[2])!=1) {
            mexErrMsgTxt("The receiver location has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[2]);
        lRx=reinterpret_cast<double*>(mxGetData(prhs[2]));
    }
    
    if(nrhs<4||mxIsEmpty(prhs[3])) {
        M=MLocal;
        M[0]=1;
        M[1]=0;
        M[2]=0;
        
        M[3]=0;
        M[4]=1;
        M[5]=0;
        
        M[6]=0;
        M[7]=0;
        M[8]=1;
    } else {
        if(mxGetM(prhs[3])!=3||mxGetN(prhs[3])!=3) {
            mexErrMsgTxt("The receiver location has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[3]);
        M=reinterpret_cast<double*>(mxGetData(prhs[3]));
    }
    
    retMat=mxCreateDoubleMatrix(2,3,mxREAL);
    retData=reinterpret_cast<double*>(mxGetData(retMat));
    
    spherAngGradientGenCPP(retData,xG,systemType,lRx, M);
    
    plhs[0]=retMat;
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
