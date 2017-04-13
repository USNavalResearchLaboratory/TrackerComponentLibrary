/**RANGEGRADIENT Determine the gradient of a 2D or 3D bistatic range
*           measurement with respect to position (gradient components for
*           velocity etc. are zero and are not provided).Atmospheric and
*           other propagation effects are not taken into account.
*
*INPUTS: x A numPosDimX1 target position vector of the form [x;y] or
*          [x;y;z].
* useHalfRange A boolean value specifying whether the bistatic range value
*          should be divided by two. This normally comes up when operating
*          in monostatic mode, so that the range reported is
*          a one-way range. The default if this parameter is not provided
*          is false.
*      lTx The 3X1 (in 3D) or 2X1 (in 2D) position vector of the
*          transmitter. If this parameter is omitted or an empty matrix is
*          passed, then the transmitter is assumed to be at the origin.
*      lRx The 3X1 (in 3D) or 2X1 (in 2D) position vector of the receiver.
*          If this parameter is omitted or an empty matrix is passed, then
*          the receiver is assumed to be at the origin.
*
*OUTPUTS: J A 1XnumPosDim gradient of the bistatic range with derivatives
*           taken with respect to components [x,y,z] in 3D or [x,y] in 2D
*           in that order.
*
*Gradients with respect to bistatic range are discussed in [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*J=rangeGradient(x,useHalfRange,lTx,lRx)
*
*REFERENCES:
*[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
*    bistatic measurements," IEEE Aerospace and Electronic Systems
*    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
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
    size_t numRows;
    double *x;
    bool useHalfRange;
    double *lTxLocal=NULL;
    double *lRxLocal=NULL;
    double *lTx, *lRx;
    double *retData;
    double tempSpace[4];
    mxArray *retMat;

    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    numRows=mxGetM(prhs[0]);
    
    if(numRows>4||mxGetN(prhs[0])!=1) {
       mexErrMsgTxt("The position has the wrong dimensionality.");
       return;
    }

    checkRealDoubleArray(prhs[0]);
    x=reinterpret_cast<double*>(mxGetData(prhs[0]));
    
    if(nrhs<2||mxIsEmpty(prhs[1])) {
        useHalfRange=false;
    } else {
        useHalfRange=getBoolFromMatlab(prhs[1]);
    }
    
    if(nrhs<3||mxIsEmpty(prhs[2])) {
        //Allocate a vector of the correct size of doubles.
        lTxLocal=reinterpret_cast<double*>(mxCalloc(numRows,sizeof(double)));
        lTx=lTxLocal;
    } else {
        if(mxGetM(prhs[2])!=numRows||mxGetN(prhs[0])!=1) {
            mexErrMsgTxt("The transmitter location has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[2]);
        lTx=reinterpret_cast<double*>(mxGetData(prhs[2]));
    }
    
    if(nrhs<4||mxIsEmpty(prhs[3])) {
        //Allocate a vector of the correct size of doubles.
        lRxLocal=reinterpret_cast<double*>(mxCalloc(numRows,sizeof(double)));
        lRx=lRxLocal;
    } else {
        if(mxGetM(prhs[3])!=numRows||mxGetN(prhs[3])!=1) {
            mexErrMsgTxt("The transmitter location has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[3]);
        lRx=reinterpret_cast<double*>(mxGetData(prhs[3]));
    }
    
    retMat=mxCreateDoubleMatrix(1,numRows,mxREAL);
    retData=reinterpret_cast<double*>(mxGetData(retMat));

    rangeGradientCPP(numRows,retData,tempSpace,x,useHalfRange,lTx,lRx);
    
    plhs[0]=retMat;
    //Free temporary memory, used.
    if(lTxLocal!=NULL) {
        mxFree(lTxLocal);
    }
    
    if(lRxLocal!=NULL) {
        mxFree(lRxLocal);
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
