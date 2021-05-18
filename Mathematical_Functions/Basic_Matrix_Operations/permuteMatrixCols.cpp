/**PERMUTEMATRIXCOLS This function is equivalent in Matlab to V=V(:,idx),
 *  where idx is a permutation vector. That is, it rearranges the columns
 *  of v according to the specified permutation. See the Matlab
 *  implementation for more details. This is just a matrix version of the
 *  function permuteVector. Currently, only matrices of doubles are
 *  supported.
 *
%INPUTS: V A numRowsXN matrix.
%      idx A length N permutation vector. This is a permutation of the
%          values 1:N.
%
%OUTPUTS: V The same as the input V but with the columns permuted
%           according to idx.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 * V=permuteMatrixCols(V,idx)
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "matrix.h"
/*This is for input validation*/
#include "MexValidation.h"
//For copy.
#include <algorithm>

//For permuteMatrixColsCPP
#include "basicMatOpsCPP.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxIsEmpty(prhs[0])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//v=[];
        return;
    }

    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    const size_t nDim1=mxGetNumberOfDimensions(prhs[0]);
    const size_t nDim2=mxGetNumberOfDimensions(prhs[1]);
    
    if(nDim1!=2) {
        mexErrMsgTxt("V has the wrong dimensionality.");
        return;
    }
    
    if(nDim2!=2) {
        mexErrMsgTxt("idx has the wrong dimensionality.");
        return;
    }
    
    const size_t numRows=mxGetM(prhs[0]);
    size_t numCols=mxGetN(prhs[0]);
        
    if(mxGetM(prhs[1])!=1&&mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("idx has the wrong dimensionality.");
        return;
    }
    
    if(numCols!=mxGetNumberOfElements(prhs[1])) {
        mexErrMsgTxt("idx must have the same number of elements as V has columns.");
        return;
    }
    
    const double *const V=mxGetDoubles(prhs[0]);
    size_t *idx=copySizeTArrayFromMatlab(prhs[1],&numCols);
    
    for(size_t i=0;i<numCols;i++) {
        if(idx[i]<1||idx[i]>numCols) {
            mexErrMsgTxt("The values in idx must be from 1 to the length of idx.");
            return;
        }
    }
    
    //Switch the indexation to C/C++ indexation.
    for(size_t k=0;k<numCols;k++) {
        idx[k]--;
    }
    
    //Copy v for the return value.
    mxArray *VMATLAB=mxCreateDoubleMatrix(numRows, numCols, mxREAL);
    double *VOut=mxGetDoubles(VMATLAB);
    std::copy(V,V+numRows*numCols,VOut);
    
    //Do the actual permutation.
    permuteMatrixColsCPP(numRows,numCols, VOut,idx);
    
    plhs[0]=VMATLAB;
    mxFree(idx);
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
