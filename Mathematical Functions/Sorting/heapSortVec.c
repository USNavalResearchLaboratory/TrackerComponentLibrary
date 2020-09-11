/**HEAPSORTVEC A C-code (for Matlab) function to sort a vector of scalar
 *         numeric values in ascending or descending order using heap sort.
 *         Heap sort has a lower worst-case performance bound O(n*log(n))
 *         than quicksort O(n^2), but a worse average case performance due
 *         to a multiplicative constant. The best case complexity is also
 *         O(n*log(n)) Note that heapSort is not a stable sorting
 *         algorithm, meaning that the order of items having the same value
 *         might change. To sort more general data types, consider the
 *         heapSort function.
 *
 *INPUTS: a The vector that is to be sorted. This can be any of the
 *          following types in Matlab: double, single, int8, int16, int32,
 *          int64, uint8, uint16, uint32, and uint64.
 * direction A value indicating the direction in which A should be sorted.
 *          Possible values are boolean:
 *          0 (The default if omitted or an empty matrix is passed) sort in
 *            ascending order.
 *          1 Sort in descending order.
 *
 *OUTPUTS: a The input vector a with its elements sorted.
 *   idxList The indices of the original vector with respect to the sorted
 *           order.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * [a,idxList]=heapSortVec(a,direction)
 *
 *The algorithm for creating and updating the heap is generally based on
 *the class implementation described in Chapter 6.4 of [1], though classes
 *are not used here. The implementation using such a heap for sorting is
 *described in Chapter 5.2.3 of [2].
 *
 *REFERENCES:
 *[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
 *    Reading, MA: Addison-Wesley, 1999.
 *[2] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
 *    ed. Reading, MA: Addison-Wesley, 1998, vol. 3.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//For memcpy.
#include <string.h>

#include "mathFuncsC.h"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    mxArray *aMATLAB;
    size_t *idxList=NULL;
    size_t M,N,numEls;
    bool direction=false;
    
    if(nrhs>2||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>2) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("a cannot be an empty matrix.");
        return;
    }
    
    if(mxIsComplex(prhs[0])) {
        mexErrMsgTxt("a must be real.");
        return;
    }
    
    if(mxGetNumberOfDimensions(prhs[0])>2) {
        mexErrMsgTxt("a must be a vector.");
        return;
    }
    
    M=mxGetM(prhs[0]);
    N=mxGetN(prhs[0]);
    
    if(M!=1&&N!=1) {
        mexErrMsgTxt("a must be a vector.");
        return;
    }
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        direction=getBoolFromMatlab(prhs[1]);   
    }
    
    //M holds the length of the array.
    if(M<N) {
        numEls=N;
    } else {
        numEls=M;
    }
    
    if(nlhs>1) {
        idxList=mxMalloc(numEls*sizeof(size_t));
    }

    switch(mxGetClassID(prhs[0])) {
        case mxDOUBLE_CLASS:
        {
            double *a, *aOrig;
            aOrig=mxGetPr(prhs[0]);
            
            aMATLAB=mxCreateDoubleMatrix(M,N,mxREAL);
            a=mxGetPr(aMATLAB);
            
            //Copy the input vector, because it will be modified.
            memcpy(a,aOrig,sizeof(double)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCDouble(numEls,a,idxList,direction);
            } else {
                heapSortVecCDouble(numEls,a,direction);
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *a, *aOrig;
            aOrig=mxGetSingles(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxSINGLE_CLASS,mxREAL); 
            a=mxGetSingles(aMATLAB);
            
            memcpy(a,aOrig,sizeof(float)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCFloat(numEls,a,idxList,direction);
            } else {
                heapSortVecCFloat(numEls,a,direction);
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_t *a, *aOrig;
            aOrig=mxGetInt8s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxINT8_CLASS,mxREAL);  
            a=mxGetInt8s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(int8_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCInt8T(numEls,a,idxList,direction);
            } else {
                heapSortVecCInt8T(numEls,a,direction);
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_t *a, *aOrig;
            aOrig=mxGetUint8s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxUINT8_CLASS,mxREAL);
            a=mxGetUint8s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(uint8_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCUInt8T(numEls,a,idxList,direction);
            } else {
                heapSortVecCUInt8T(numEls,a,direction);
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_t *a, *aOrig;
            aOrig=mxGetInt16s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxINT16_CLASS,mxREAL);
            a=mxGetInt16s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(int16_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCInt16T(numEls,a,idxList,direction);
            } else {
                heapSortVecCInt16T(numEls,a,direction);
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_t *a, *aOrig;
            aOrig=mxGetUint16s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxUINT16_CLASS,mxREAL);
            a=mxGetUint16s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(uint16_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCUInt16T(numEls,a,idxList,direction);
            } else {
               heapSortVecCUInt16T(numEls,a,direction);
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_t *a, *aOrig;
            aOrig=mxGetInt32s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxINT32_CLASS,mxREAL);
            a=mxGetInt32s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(int32_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCInt32T(numEls,a,idxList,direction);
            } else {
                heapSortVecCInt32T(numEls,a,direction);
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_t *a, *aOrig;
            aOrig=mxGetUint32s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxUINT32_CLASS,mxREAL);
            a=mxGetUint32s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(uint32_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCUInt32T(numEls,a,idxList,direction);
            } else {
                heapSortVecCUInt32T(numEls,a,direction);
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_t *a, *aOrig;
            aOrig=mxGetInt64s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxINT64_CLASS,mxREAL);
            a=mxGetInt64s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(int64_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCInt64T(numEls,a,idxList,direction);
            } else {
                heapSortVecCInt64T(numEls,a,direction);
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_t *a, *aOrig;
            aOrig=mxGetUint64s(prhs[0]);
            
            aMATLAB=mxCreateNumericMatrix(M,N,mxUINT64_CLASS,mxREAL);
            a=mxGetUint64s(aMATLAB);
            
            memcpy(a,aOrig,sizeof(uint64_t)*numEls);
            if(nlhs>1) {
                heapSortVecIdxCUInt64T(numEls,a,idxList,direction);
            } else {
                heapSortVecCUInt64T(numEls,a,direction);
            }
            break;
        }
        case mxCHAR_CLASS:
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxLOGICAL_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("a has an unsupported data type");
            return;
    }

    plhs[0]=aMATLAB;
    if(nlhs>1) {
        size_t i;
        
        //Convert from C indices to Matlab indices.
        for(i=0;i<numEls;i++) {
            idxList[i]++;
        }

        plhs[1]=sizeTMat2MatlabDoubles(idxList,1,numEls);
        
        mxFree(idxList);
    }
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
