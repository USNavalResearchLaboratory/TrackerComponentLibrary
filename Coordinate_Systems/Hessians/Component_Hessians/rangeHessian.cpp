/**RANGEHESSIAN Determine the Hessian matrix (a matrix of second-order
*          partial derivatives) of a bistatic range rate measurement with
*          respect to Cartesian position. Atmospheric and other propagation
*          effects are not taken into account. More details are given in
*          the native Matlab implementation.
*
*The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
* H=rangeHessian(x,useHalfRange,lTx,lRx);
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
    size_t numDim,N;
    double *H;
    const double *x, *lTx, *lRx;
    const double zeros[]={0,0,0};
    bool useHalfRange;
    mxArray *HMATLAB;
    size_t curPoint, numDimIncr;

    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Wrong number of inputs.");
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    numDim=mxGetM(prhs[0]);
    
    if(numDim!=1&&numDim!=2&&numDim!=3) {
        mexErrMsgTxt("x has an unsupported dimensionality.");
    }

    //The number of elements in each H matrix.
    numDimIncr=numDim*numDim;

    if(mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("x should not be empty.");
    }

    N=mxGetN(prhs[0]);
    checkRealDoubleArray(prhs[0]);
    x=mxGetDoubles(prhs[0]);
    
    if(nrhs>1) {
        useHalfRange=getBoolFromMatlab(prhs[1]);
    } else {
        useHalfRange=false;
    }
    
    //Allocate space for the return values
    {
        mwSize dims[3];
  
        dims[0]=numDim;
        dims[1]=numDim;
        dims[2]=N;
            
        HMATLAB=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        H=mxGetDoubles(HMATLAB);
    }
    
    if(nrhs>2) {
        if(!mxIsEmpty(prhs[2])) {
            if(mxGetM(prhs[2])!=numDim||mxGetN(prhs[2])!=1) {
                mexErrMsgTxt("lRx has the wrong dimensionality.");
            }

            checkRealDoubleArray(prhs[2]);
            lTx=mxGetDoubles(prhs[2]);
        } else {
            lTx=zeros;
        }

        if(nrhs>3&&!mxIsEmpty(prhs[3])) {
            if(mxGetM(prhs[3])!=numDim||mxGetN(prhs[3])!=1) {
                mexErrMsgTxt("lRx has the wrong dimensionality.");
            }

            checkRealDoubleArray(prhs[3]);
            lRx=mxGetDoubles(prhs[3]);
        } else {
            lRx=zeros;
        }

        for(curPoint=0;curPoint<N;curPoint++) {
            rangeHessianGenCPP(numDim,H,x,useHalfRange,lTx,lRx); 

            x+=numDim;
            H+=numDimIncr;
        }
    } else {
        //If it is monostatic.
        for(curPoint=0;curPoint<N;curPoint++) {
            rangeHessianCPP(numDim,H,x,useHalfRange); 

            x+=numDim;
            H+=numDimIncr;
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
