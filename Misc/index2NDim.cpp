/**INDEX2NDIM A C++ implemetation of a function such that given the number
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
 *             the nDim dimensions. The values are >=1.
 *         idx A numIdxX1 or 1XnumIdx vector of linear indices (starting
 *             from 1, not 0).
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
 *indices=index2NDim(dims,idx);
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
//For std::copy
#include <algorithm>
//For floor
#include <cmath>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double *dims, *idxIn, *idx, *indices, *maxVals;
    double dimProd=1;
    size_t i,curIndic, numIdx, numDim;
    mxArray *retMat;
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }

    numIdx=mxGetNumberOfElements(prhs[1]);
    numDim=mxGetNumberOfElements(prhs[0]);
    
    //If an empty matrix is passed, return an empty matrix.
    if(numIdx==0) {
        plhs[0]=mxCreateDoubleMatrix(numDim,0,mxREAL);
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    dims=mxGetDoubles(prhs[0]);
    idxIn=mxGetDoubles(prhs[1]);
    
    //Determine the total number of elements implied by the dimensions.
    for(i=0;i<numDim;i++) {
        dimProd*=dims[i];
    }
    
    //Return an empty matrix if an index is beyond the end of the array.
    for(i=0;i<numIdx;i++) {
        if(idxIn[i]>dimProd) {
            plhs[0]=mxCreateDoubleMatrix(numDim,0,mxREAL);
            return;
        }
    }
    
    //We have to duplicate the idx matrix to make changes to it.
    idx=reinterpret_cast<double*>(mxMalloc(numIdx*sizeof(double)));
    std::copy(idxIn,idxIn+numIdx,idx);
    
    retMat=mxCreateDoubleMatrix(numDim,numIdx,mxREAL);
    indices=mxGetDoubles(retMat);
    
    //The indices basically form a counting system where the bases of each
    //position are determined by dim. Matlab makes the first dimension the
    //least significant. Thus, finding the indices is done in the same
    //manner as one would decompose a number into digits or bits.
    maxVals=reinterpret_cast<double*>(mxMalloc(numDim*sizeof(double)));
    maxVals[0]=1;
    for(i=1;i<numDim;i++) {
        maxVals[i]=maxVals[i-1]*dims[i-1];
    }

    curIndic=numDim;
    do {
        curIndic--;
        for(i=0;i<numIdx;i++) {
            double wholeVal=floor((idx[i]-1)/maxVals[curIndic]);
            indices[i*numDim+curIndic]=wholeVal+1;
            idx[i]=idx[i]-wholeVal*maxVals[curIndic];
        }
    } while(curIndic>0);

    mxFree(idx);    
    mxFree(maxVals);
    plhs[0]=retMat;     
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
