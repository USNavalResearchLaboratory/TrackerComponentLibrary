/**GETNEXTTUPLE This function produces the next tuple in a counting series.
 *          This can be used as a method of implementing nested loops. Each
 *          index of the tuple counts from 0 to the value in maxVals and
 *          resets. This is essentially a function for counting using
 *          different bases for each digit.
 *
 *See the comments to the Matlab implemention for more details on the
 *algorithm.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * tuple=getNextTuple(param1,maxVals,firstIsMostSig)
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool firstIsMostSig=true;
    size_t numDim;
    size_t *tuple;
    size_t *maxVals;

    if(nrhs>3||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("The previous tuple cannot be empty.");
        return;
    }
    
    //If the first tuple is requested.
    if(nrhs==1) {
        const size_t numEls=getSizeTFromMatlab(prhs[0]);
        mxArray *retMat=mxCreateDoubleMatrix(numEls, 1, mxREAL);
        double *retVals=mxGetDoubles(retMat);
        size_t i;
        
        for(i=0;i<numEls;i++) {
            retVals[i]=0;
        }
        
        plhs[0]=retMat;
        return;
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        firstIsMostSig=getBoolFromMatlab(prhs[2]);
    }
    
    //Get the first two inputs
    tuple=copySizeTArrayFromMatlab(prhs[0],&numDim);
    {
        size_t arrayLen;
        maxVals=copySizeTArrayFromMatlab(prhs[1], &arrayLen);
        
        if(arrayLen<numDim) {
            mexErrMsgTxt("The dimensionalities of the first two inputs are not the same.");
            return;
        }
    }
    
    bool isPastEnd=getNextTupleCPP(numDim,tuple,maxVals,firstIsMostSig);
    
    if(isPastEnd) {
        //Return an empty matrix.
        plhs[0]=mxCreateDoubleMatrix(0, 0,mxREAL);
    } else {
        plhs[0]=sizeTMat2MatlabDoubles(tuple,numDim, 1);
    }

    mxFree(tuple);
    mxFree(maxVals);
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
