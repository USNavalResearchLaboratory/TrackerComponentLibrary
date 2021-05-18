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

bool index2NDimCPP(const size_t numDim, const size_t numIdx,const size_t * const idx, const size_t * const maxVals, size_t *newIndices) {
/**INDEX2NDIMCPP Given the number of dimensions of a multidimensional array
 *           and the linear index of an item in the array (starting from 1,
 *           not 0), find the set of indices (starting from 1) that
 *           addresses the point in the array. This function is like the
 *           ind2sub function, but it returns a single vector when dealing
 *           with multidimensional arrays rather than returning multiple
 *           outputs.
 *
 *INPUTS: numDim The dimensionality of the indices.
 *        numIdx The number of indices provided that should be unranked to
 *               form tuples.
 *           idx A length numIdx array holding the indices to unrank. These
 *               start from 1, not 0.
 *       maxVals The product values related to the dimensions, as described
 *               above.
 *    newIndices A numDim*numIdx array into which the unranked tuples will
 *               be placed. These are placed one after the other.
 *
 *OUTPUTS: The return value indicates whether one is past the final indexed
 *         value. If idx is too large (or 0, which is invalid) then nothing
 *         is put in newIndices and this function returns true. Otherwise,
 *         the function returns false and newIndices holds the new set of
 *         indices.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t maxVal=maxVals[numDim];
    
    for(size_t i=0;i<numIdx;i++) {
        size_t curIndic=numDim;
        size_t curIdx=idx[i];
        
        //The index is beyond the final one or is invalid.
        if(curIdx>maxVal||curIdx==0) {
            return true;
        }
        
        do {
            curIndic--;
            //integer division.
            size_t wholeVal=(curIdx-1)/maxVals[curIndic];
            
            newIndices[i*numDim+curIndic]=wholeVal+1;
            curIdx-=wholeVal*maxVals[curIndic];
        } while(curIndic>0);
    }
    return false;
}

bool index2NDimByDimCPP(const size_t numDim, const size_t numIdx,const size_t *idx, const size_t *dims, size_t *newIndices, size_t *tempSpace) {
/**INDEX2NDIMBYDIMCPP This is the same as the index2NDimCPP function,
 *     except it is parameterized by the dimensions (dims) of the index
 *     vector (dims) rather than by the product vector that is used in
 *     index2NDimCPP. The tempSpace buffer must have enough space to hold
 *     at least numDim+1 size_t values.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    //This memory buffer is numDim+1 in size.
    size_t *maxVals=tempSpace;        
    maxVals[0]=1;
    for(size_t i=1;i<numDim+1;i++) {
        maxVals[i]=maxVals[i-1]*dims[i-1];
    }
    return index2NDimCPP(numDim, numIdx, idx, maxVals, newIndices);
}

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
