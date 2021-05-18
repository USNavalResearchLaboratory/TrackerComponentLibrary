/**GENALLTUPLES Generate all tuples in a counting series. Each index of the
 *              tuple counts from 0 to the corresponding value in maxVals.
 *              This is just a method of counting uses different bases for
 *              each digit.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * theTuples=genAllTuples(maxVals,firstIsMostSig)
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "combinatorialFuns.hpp"
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool firstIsMostSig=true;
    
    if(nrhs>2||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        firstIsMostSig=getBoolFromMatlab(prhs[1]);
    }
    
    size_t N;
    size_t *maxVals=copySizeTArrayFromMatlab(prhs[0], &N);
    const size_t bytesNeeded=bufferSize4AllTuples<size_t>(N, maxVals); 
    //Alllocate space for the matrix.
    size_t *theTuples=static_cast<size_t *>(mxMalloc(bytesNeeded));
    const size_t tupleCount=genAllTuplesCPP(theTuples,N,maxVals, firstIsMostSig);
    
    plhs[0]=sizeTMat2MatlabDoubles(theTuples,N,tupleCount);
    mxFree(theTuples);
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

