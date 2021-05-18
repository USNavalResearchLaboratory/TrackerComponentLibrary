/**PERMUTEVECTOR This function is equivalent in Matlab to v=v(idx), where
 *  idx is a permutation vector. That is, it rearranges the elements of v
 *  according to the specified permutation. See the Matlab implementation
 *  for more details. Currently, only matrices of doubles are supported.
 *
 *INPUTS: v A 1XN or NX1 vector.
 *      idx A length N permutation vector. This is a permutation of the
 *         values 1:N.
 *
 *OUTPUTS: v The same as the input v but with the elements permuted
 *          according to idx.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 * v=permuteVector(v,idx)
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
//For permuteVectorCPP
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
        mexErrMsgTxt("v has the wrong dimensionality.");
        return;
    }
    
    if(nDim2!=2) {
        mexErrMsgTxt("idx has the wrong dimensionality.");
        return;
    }
    
    const size_t numRows=mxGetM(prhs[0]);
    const size_t numCols=mxGetN(prhs[0]);
    
    if(numRows!=1&&numCols!=1) {
        mexErrMsgTxt("v has the wrong dimensionality.");
        return;
    }
    
    if(mxGetM(prhs[1])!=1&&mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("idx has the wrong dimensionality.");
        return;
    }
    
    size_t numEls=mxGetNumberOfElements(prhs[0]);

    if(numEls!=mxGetNumberOfElements(prhs[1])) {
        mexErrMsgTxt("The lengths of v and idx do not agree.");
        return;
    }
    
    const double *const v=mxGetDoubles(prhs[0]);
    size_t *idx=copySizeTArrayFromMatlab(prhs[1],&numEls);
    
    for(size_t i=0;i<numEls;i++) {
        if(idx[i]<1||idx[i]>numEls) {
            mexErrMsgTxt("The values in idx must be from 1 to the length of idx.");
            return;
        }
    }
    
    //Switch the indexation to C/C++ indexation.
    for(size_t k=0;k<numEls;k++) {
        idx[k]--;
    }
    
    //Copy v for the return value.
    mxArray *vMATLAB=mxCreateDoubleMatrix(numRows, numCols, mxREAL);
    double *vOut=mxGetDoubles(vMATLAB);
    std::copy(v,v+numEls,vOut);
    
    //Do the actual permutation.
    permuteVectorCPP(numEls,vOut,idx);
    
    plhs[0]=vMATLAB;
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
