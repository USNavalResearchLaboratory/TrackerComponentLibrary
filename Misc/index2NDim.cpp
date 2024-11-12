/**INDEX2NDIM A C++ implementation of a function such that given the number
 *           of dimensions of a multidimensional array and the linear index
 *           of an item in the array, the set of indices that addresses the
 *           point in the array is found. This function is like the ind2sub
 *           function, but it returns a single vector when dealing with
 *           multidimensional arrays rather than returning multiple
 *           outputs. This function is useful if one wishes to write a
 *           function where the number of dimensions needed in a
 *           hypermatrix is not known a priori, such as when handling
 *           multivariate polynomials. Inputs for this implementation must
 *           be doubles.
 *
 *INPUTS: dims A 1XnumDim or numDimX1 vector holding the size of each of
 *            the nDim dimensions. The values are >=1. Alternatively, if
 *            isProdVals=true, these can be the equivalent of 
 *            [1;cumprod(dims)] rather than the dimensions directly.
 *            Passing the product values rather than the dimensions will
 *            speed up the function.
 *        idx A numIdxX1 or 1XnumIdx vector of linear indices (starting
 *            from 1, not 0).
 * areProdVals An optional boolean parameter indicating whether dims is a
 *            set of dimensions or product values, as described above. The
 *            default if omitted or an empty matrix is passed is false.
 *
 *OUTPUTS: indices A numDimXnumIdx matrix of sets of indices corresponding
 *                 to each of the numIdx linear indices. Indexation starts
 *                 from 1, not 0. If the index is above the maximum
 *                 number of elements, an empty matrix is returned.
 *
 *The function nDim2Index is the inverse of this function.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *indices=index2NDim(dims,idx,areProdVals)
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "miscFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool areProdVals=false;
    
    if(nrhs<2||nrhs>3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        areProdVals=getBoolFromMatlab(prhs[2]);
    }

    //If an empty matrix is passed, return an empty matrix.
    if(mxIsEmpty(prhs[1])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

    size_t numDimEls=mxGetNumberOfElements(prhs[0]);
    size_t *dims=copySizeTArrayFromMatlab(prhs[0],&numDimEls);
    size_t numIdx=mxGetNumberOfElements(prhs[1]);
    size_t *idx=copySizeTArrayFromMatlab(prhs[1],&numIdx);
    size_t numDim;
    
    size_t *newIndices;
    bool isPastEnd;
    if(areProdVals==true) {
        numDim=numDimEls-1;
        size_t *maxVals=dims;
        newIndices=new size_t[numDim*numIdx];
        
        isPastEnd=index2NDimCPP(numDim, numIdx, idx, maxVals, newIndices);
    } else {
        numDim=numDimEls;
        newIndices=new size_t[numDim*numIdx];
        size_t *tempSpace=new size_t[numDim+1];
        
        isPastEnd=index2NDimByDimCPP(numDim, numIdx, idx, dims, newIndices, tempSpace);
           
        delete[] tempSpace;
    }
    
    mxFree(idx);    
    mxFree(dims);

    //Return an empty matrix if an index is beyond the end of the array.
    //Also, idx cannot be 0. That is an error.
    if(isPastEnd) {
        delete[] newIndices;
        plhs[0]=mxCreateDoubleMatrix(numDim,0,mxREAL);
        return;
    }
    
    plhs[0]=sizeTMat2MatlabDoubles(newIndices,numDim,numIdx);
    
    delete[] newIndices;
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
