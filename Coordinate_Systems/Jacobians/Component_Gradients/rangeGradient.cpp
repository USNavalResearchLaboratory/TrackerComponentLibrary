/**RANGEGRADIENT Determine the gradient of a 2D or 3D bistatic range
*           measurement with respect to position. See the comments ot the
*           Matlab version for more details.
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
    size_t numRows,N,i;
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
    N=mxGetN(prhs[0]);
    
    if(numRows>4||mxIsEmpty(prhs[0])) {
       mexErrMsgTxt("The position has the wrong dimensionality.");
       return;
    }

    checkRealDoubleArray(prhs[0]);
    x=mxGetDoubles(prhs[0]);
    
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
        if(mxGetM(prhs[2])!=numRows||mxGetN(prhs[2])!=1) {
            mexErrMsgTxt("The transmitter location has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[2]);
        lTx=mxGetDoubles(prhs[2]);
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
        lRx=mxGetDoubles(prhs[3]);
    }
    
    {
        mwSize dims[3];
        dims[0]=1;
        dims[1]=numRows;
        dims[2]=N;
        
        retMat=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    }
    
    retData=mxGetDoubles(retMat);

    for(i=0;i<N;i++) {
        rangeGradientCPP(numRows,retData,tempSpace,x,useHalfRange,lTx,lRx);
        retData+=numRows;
        x=x+numRows;
    }
    
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
