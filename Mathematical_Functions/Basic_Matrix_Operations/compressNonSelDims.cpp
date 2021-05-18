/**COMPRESSNONSELDIMS Given a hypermatrix, reduce the number of dimensions
 *      keeping the selected dimensions as distrinct dimensions. For
 *      example, if selDim=3, then a matrix of the form
 *      C(:,:,selDim,:,:,:) (where we marked the selected dimensions), will
 *      be reduced to the form C(:,selDim,:). See the Matlab implementation
 *      for more details.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * [C,newSelDims]=compressNonSelDims(C,selDims,justDims)
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "basicMatOpsCPP.hpp"
//For max.
#include <algorithm>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool justDims=false;
    size_t *nDims;
    size_t S;
    
    if(nrhs>3||nrhs<2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        justDims=getBoolFromMatlab(prhs[2]);
    }
    
    if(mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("The C input cannot be empty.");
        return;
    }

    if(justDims) {
        nDims=copySizeTArrayFromMatlab(prhs[0],&S);
    } else {
        S=mxGetNumberOfDimensions(prhs[0]);
        const size_t two=2;
        const size_t *const nDimsTemp=mxGetDimensions(prhs[0]);
        
        //2, because Matlab requires that all things are 2D, even if there
        //is a singleton dimension.
        nDims=new size_t[std::max(S,two)];
        std::copy(nDimsTemp,nDimsTemp+S,nDims);
    }
    
    mxArray *selDimsRetVal;
    size_t numNewDims;

    if(mxIsLogical(prhs[1])) {
        if(mxGetNumberOfElements(prhs[1])!=S) {
            mexErrMsgTxt("If selDims is logical, it must have the same number of elements as the number of dimensions of C.");
            return;
        }

        selDimsRetVal=mxDuplicateArray(prhs[1]);
        mxLogical *selDims=mxGetLogicals(selDimsRetVal);

        compressNonSelDimsBool(S,nDims,selDims,numNewDims);
        
        //Adjust to the updated length of selDimsRetVal.
        mxSetM(selDimsRetVal,1);
        mxSetN(selDimsRetVal,numNewDims);
    } else {
        //If here, we were passed a number of indices, not a boolean
        //array. Assume they are doubles, since that is the default
        //format in Matlab.
        const size_t numSel=mxGetNumberOfElements(prhs[1]);
        if(numSel>S) {
            mexErrMsgTxt("The C input cannot be empty.");
            return;
        }

        selDimsRetVal=mxDuplicateArray(prhs[1]);
        double *selDims=mxGetDoubles(selDimsRetVal);
        for(size_t k=0;k<numSel;k++) {
            if(selDims[k]<1||selDims[k]>S) {
                mexErrMsgTxt("selDims is invalid.");
                return;
            }
            selDims[k]--;//Convert to indexation from 0.
        }

        compressNonSelDimsCPP(S,nDims,selDims,numSel,numNewDims);

        //Convert to Matlab indexation.
        for(size_t k=0;k<numSel;k++) {
            selDims[k]++;
        }
    }

    if(numNewDims==1) {
        //Matlab requires everything to be at least two dimensions, so add
        //a trailing singleton dimension.
        numNewDims=2;
        nDims[1]=1;
    }
    
    if(justDims) {
        plhs[0]=sizeTMat2MatlabDoubles(nDims,1,numNewDims);
        mxFree(nDims);
    } else {
        mxArray *CMat=mxDuplicateArray(prhs[0]);
        mxSetDimensions(CMat, nDims, numNewDims);
        //We cannot free nDims here; it is taken, not copied.
        plhs[0]=CMat;
    }

    if(nlhs>1) {
        plhs[1]=selDimsRetVal;
    }else {
        mxDestroyArray(selDimsRetVal);
    }
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
