/**MATHFUNCSC A header file for C implementations of mathematical
 *            functions. See the files implementing each function for more
 *            details on their usage.
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MATHFUNCSC
#define MATHFUNCSC

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __cplusplus
#if __STDC_VERSION__>=199901L
#include <stdbool.h>
#else
#ifndef _bool_T
#define false 0
#define true 1
#define bool int
#define _bool_T
#endif
#endif
#endif

//Defines the size_t and ptrdiff_t types
#include <stddef.h>
//This defines the standard integer data types.
#include <stdint.h>

size_t binSearchC(size_t numInVec, double *vec, double key,int choice);

//heapSort for different data types when sorting with an index vector
void heapSortVecIdxCDouble(size_t numInHeap, double *a, size_t *idxList, const bool direction);
void heapSortVecIdxCFloat(size_t numInHeap, float *a, size_t *idxList, const bool direction);
void heapSortVecIdxCInt(size_t numInHeap, int *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUInt(size_t numInHeap, unsigned int *a, size_t *idxList, const bool direction);
void heapSortVecIdxCSizeT(size_t numInHeap, size_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCPtrDiffT(size_t numInHeap, ptrdiff_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCChar(size_t numInHeap, char *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUChar(size_t numInHeap, unsigned char *a, size_t *idxList, const bool direction);
void heapSortVecIdxCInt8T(size_t numInHeap, int8_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUInt8T(size_t numInHeap, uint8_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCInt16T(size_t numInHeap, int16_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUInt16T(size_t numInHeap, uint16_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCInt32T(size_t numInHeap, int32_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUInt32T(size_t numInHeap, uint32_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCInt64T(size_t numInHeap, int64_t *a, size_t *idxList, const bool direction);
void heapSortVecIdxCUInt64T(size_t numInHeap, uint64_t *a, size_t *idxList, const bool direction);

//heapSort for different data types when sorting without an index vector
void heapSortVecCDouble(size_t numInHeap, double *a, const bool direction);
void heapSortVecCFloat(size_t numInHeap, float *a, const bool direction);
void heapSortVecCInt(size_t numInHeap, int *a, const bool direction);
void heapSortVecCUInt(size_t numInHeap, unsigned int *a, const bool direction);
void heapSortVecCSizeT(size_t numInHeap, size_t *a, const bool direction);
void heapSortVecCPtrDiffT(size_t numInHeap, ptrdiff_t *a, const bool direction);
void heapSortVecCChar(size_t numInHeap, char *a, const bool direction);
void heapSortVecCUChar(size_t numInHeap, unsigned char *a, const bool direction);
void heapSortVecCInt8T(size_t numInHeap, int8_t *a, const bool direction);
void heapSortVecCUInt8T(size_t numInHeap, uint8_t *a, const bool direction);
void heapSortVecCInt16T(size_t numInHeap, int16_t *a, const bool direction);
void heapSortVecCUInt16T(size_t numInHeap, uint16_t *a, const bool direction);
void heapSortVecCInt32T(size_t numInHeap, int32_t *a, const bool direction);
void heapSortVecCUInt32T(size_t numInHeap, uint32_t *a, const bool direction);
void heapSortVecCInt64T(size_t numInHeap, int64_t *a, const bool direction);
void heapSortVecCUInt64T(size_t numInHeap, uint64_t *a, const bool direction);

#ifdef __cplusplus
}
#endif

#endif

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
