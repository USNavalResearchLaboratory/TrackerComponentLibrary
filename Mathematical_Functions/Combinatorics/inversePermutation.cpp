/**INVERSEPERMUTATION Given a permutation of the number 1:N, return the
*          inverse of the permutation. If something were rearranged with
*          the original permutation, the inverse permutation puts it back
*          into the original order. See the Matlab implementation for
*          more comments.
*
* The algorithm can be compiled for use in Matlab using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* xPerm=inversePermutation(xPerm,algorithm)
*
*October 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. */

#include "MexValidation.h"
/*This header is required by Matlab*/
#include "mex.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t algorithm=0;
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>2){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }

    if(mxIsEmpty(prhs[0])){
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        algorithm=getSizeTFromMatlab(prhs[1]);
    }

    checkRealArray(prhs[0]);
    size_t N;
    ptrdiff_t *xPerm=copyPtrDiffTArrayFromMatlab(prhs[0],&N);

    //Make sure that the numbers 0->N-1 are preset in xPerm exactly once.
    bool *temp=(bool*)mxCalloc(N,sizeof(bool));//All initialized to 0
    for(size_t k=0;k<N;k++) {
        const ptrdiff_t curVal=xPerm[k]-1; 
        if(curVal>=0||curVal<static_cast<ptrdiff_t>(N)) {
            temp[curVal]=true;
        } else {
             mexErrMsgTxt("xPerm contains a value that is too large.");
        }
    }
    //Sum up how many values are present.
    size_t sumVal=0;
    for(size_t k=0;k<N;k++) {
        sumVal+=temp[k];
    }
    if(sumVal!=N) {
        mexErrMsgTxt("xPerm does not contain all of the values from 1 to N.");
    }

    //If here, we know that xPerm contains the values from 0 to N-1 and
    //thus no issues with reading past the end of arrays should arise.
    mxFree(temp);
    switch(algorithm)
    {
        case 0://Algorithm I in [1]. An in-place algorithm. It relies on
               //The indexation starting from 1, not 0.
            inversePermutationI(xPerm,N,true);
            break;
        case 1://Algorithm J in [1]. An in-place algorithm.
            inversePermutationJ(xPerm,N,true);
            break;
        case 2://A non-in-place algorithm.
        {
            //We must allocate something to hold the result.
            ptrdiff_t *xInvPerm=(ptrdiff_t*)mxMalloc(N*sizeof(ptrdiff_t));
            
            inversePermutationExpl(xInvPerm, xPerm,N,true);
            mxFree(xPerm);
            xPerm=xInvPerm;
            break;
        }
        default:
            mexErrMsgTxt("Unknown algorithm specified.");
    }

    //Return the result as Matlab doubles.
    plhs[0]=ptrDiffTMat2MatlabDoubles(xPerm,mxGetM(prhs[0]), mxGetN(prhs[0]));
    mxFree(xPerm);
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
