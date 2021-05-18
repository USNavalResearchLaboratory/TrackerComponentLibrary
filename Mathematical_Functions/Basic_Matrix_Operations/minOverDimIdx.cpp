/**MINOVERDIMIDX Given an S-dimensional matrix C, find the minimum value in
 *      the matrix when the index in dimension selDim is fixed to selIdx.
 *      For example, if C is 5 dimensional and selDim=4this function
 *      evaluates the equivalent of min(vec(C(:,:,:,selIdx,:))). This
 *      function can be useful when the dimensionality of the matrix in
 *      question can vary. See the Matlab implementation for more details.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *minVal=minOverDimIdx(C,selDim,selIdx)
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//For min
#include <algorithm>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    if(nrhs!=3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }
    
    checkRealDoubleHypermatrix(prhs[0]);
    
    const size_t S=mxGetNumberOfDimensions(prhs[0]);
    const size_t * const nDims=mxGetDimensions(prhs[0]);
    const double * const C=mxGetDoubles(prhs[0]);
    size_t selDim=getSizeTFromMatlab(prhs[1]);
    size_t selIdx=getSizeTFromMatlab(prhs[2]);//Indexation from 0.
    
    if(selDim>S||selDim<1) {
        mexErrMsgTxt("selDim is invalid.");
        return;  
    }
    //Indexation from 0 instead of 1.
    selDim--;
    
    if(selIdx<1||selIdx>nDims[selDim]) {
        mexErrMsgTxt("selIdx is invalid.");
        return;  
    }
    
    //Indexation from 0 instead of 1.
    selIdx--;
        
    double minVal;
    
    if(S==1) {
        minVal=C[selIdx];
    } else if(selDim==0) {
        const size_t startDim=nDims[0];
        size_t endDim=nDims[1];
        for(size_t k=2;k<S;k++) {
            endDim*=nDims[k];
        }

        minVal=C[selIdx];
        size_t curIdx=selIdx;
        for(size_t k=1;k<endDim;k++) {
            curIdx+=startDim;
            minVal=std::min(minVal,C[curIdx]);
        }
    } else if(selDim==S-1) {
        size_t startDim=nDims[0];
        for(size_t k=1;k<S-1;k++) {
            startDim*=nDims[k];
        }
        
        size_t offset=startDim*selIdx;
        minVal=C[offset];

        for(size_t k=1;k<startDim;k++) {
            offset++;
            minVal=std::min(minVal,C[offset]);
        }
    } else {
        size_t startDim=nDims[0];
        for(size_t k=1;k<selDim;k++) {
            startDim*=nDims[k];
        }
        size_t endDim=nDims[selDim+1];
        for(size_t k=selDim+2;k<S;k++) {
            endDim*=nDims[k];
        }
        
        const size_t offsetFixed=startDim*selIdx;
        const size_t stepEnd=startDim*nDims[selDim];
        minVal=C[offsetFixed];
        
        size_t outerOffset=offsetFixed;
        for(size_t idx2=0;idx2<endDim;idx2++) {
            for(size_t idx1=0;idx1<startDim;idx1++) {
                minVal=std::min(minVal,C[outerOffset+idx1]);
            }
            outerOffset+=stepEnd;
        }
    }
    
    plhs[0]=mxCreateDoubleScalar(minVal);
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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
*/
