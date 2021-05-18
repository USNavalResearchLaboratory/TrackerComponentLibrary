 /**HEAPSORTVECC A set of C-only implementations of a function to sort a
 *              vector in ascending or descending order and keep track of
 *              how the indices of the function change. The bulk of the
 *              code for the function is written as macros that are filled
 *              in with the specific data type of the vector being sorted
 *              in each function. See the Matlab function heapSortVec for
 *              more details on the algorithm.
 *
 *Each of the functions, such as heapSortVecCIdxSizeT, has the same inputs.
 *The difference is only in the type of the input a, which is specified by
 *the name of the function. The inputs to the functions have the following
 *meaning:
 *numInHeap The size_t number of items to sort. This is the length of a and
 *          idxList.
 *        a The array whose elements will be sorted (in place).
 *  idxListHS An array of size_t elements that will hold the indices of the
 *          original vector with respect to the sorted order. Indexation
 *          begins at 0.
 * direction A boolean variable where 0 indicates sorting in ascending
 *          order and 1 sorting in descending order.
 *
 *For the functions without idx in their names, the input idxListHS is omitted.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncsC.h"

#define percolateDownIdxMacro(a,idxListHS,idx,numInHeap,direction) {\
    size_t tempIdxPerc;\
    size_t idxP=idx;\
\
    tempP=a[idxP];\
    tempIdxPerc=idxListHS[idxP];\
    if(direction==0) {\
        while(2*idxP<numInHeap) {\
            size_t child=2*idxP;\
\
            if(child!=(numInHeap-1)&&a[child+1]>a[child]) {\
                child++;\
            }\
\
            if(a[child]>tempP) {\
                a[idxP]=a[child];\
                idxListHS[idxP]=idxListHS[child];\
            } else {\
                break;\
            }\
\
            idxP=child;\
        }\
    } else {\
        while(2*idxP<numInHeap) {\
            size_t child=2*idxP;\
\
            if(child!=(numInHeap-1)&&a[child+1]<a[child]) {\
                child++;\
            }\
\
            if(a[child]<tempP) {\
                a[idxP]=a[child];\
                idxListHS[idxP]=idxListHS[child];\
            } else {\
                break;\
            }\
\
            idxP=child;\
        }\
    }\
\
    a[idxP]=tempP;\
    idxListHS[idxP]=tempIdxPerc;\
}

#define heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction) {\
    size_t i, idx;\
    \
    for(i=0;i<numInHeap;i++) {\
        idxListHS[i]=i;\
    }\
\
    idx=numInHeap/2+1;\
    do {\
        idx--;\
        percolateDownIdxMacro(a,idxListHS,idx,numInHeap,direction);\
    } while(idx>0);\
\
    while(numInHeap>1) {\
        size_t tempIdx;\
        temp=a[0];\
        a[0]=a[numInHeap-1];\
        a[numInHeap-1]=temp;\
\
        tempIdx=idxListHS[0];\
        idxListHS[0]=idxListHS[numInHeap-1];\
        idxListHS[numInHeap-1]=tempIdx;\
\
        numInHeap--;\
\
        percolateDownIdxMacro(a,idxListHS,0,numInHeap,direction);\
    }\
}

#define percolateDownMacro(a,idx,numInHeap,direction) {\
    size_t idxP=idx;\
\
    tempP=a[idxP];\
    if(direction==0) {\
        while(2*idxP<numInHeap) {\
            size_t child=2*idxP;\
\
            if(child!=(numInHeap-1)&&a[child+1]>a[child]) {\
                child++;\
            }\
\
            if(a[child]>tempP) {\
                a[idxP]=a[child];\
            } else {\
                break;\
            }\
\
            idxP=child;\
        }\
    } else {\
        while(2*idxP<numInHeap) {\
            size_t child=2*idxP;\
\
            if(child!=(numInHeap-1)&&a[child+1]<a[child]) {\
                child++;\
            }\
\
            if(a[child]<tempP) {\
                a[idxP]=a[child];\
            } else {\
                break;\
            }\
\
            idxP=child;\
        }\
    }\
\
    a[idxP]=tempP;\
}

#define heapSortVecCMacro(numInHeap,a,direction) {\
    size_t idx;\
\
    idx=numInHeap/2+1;\
    do {\
        idx--;\
        percolateDownMacro(a,idx,numInHeap,direction);\
    } while(idx>0);\
\
    while(numInHeap>1) {\
        temp=a[0];\
        a[0]=a[numInHeap-1];\
        a[numInHeap-1]=temp;\
\
        numInHeap--;\
\
        percolateDownMacro(a,0,numInHeap,direction);\
    }\
}

void heapSortVecIdxCDouble(size_t numInHeap, double *a, size_t *idxListHS, const bool direction) {
    double temp;
    double tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCFloat(size_t numInHeap, float *a, size_t *idxListHS, const bool direction) {
    float temp;
    float tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCInt(size_t numInHeap, int *a, size_t *idxListHS, const bool direction) {
    int temp;
    int tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUInt(size_t numInHeap, unsigned int *a, size_t *idxListHS, const bool direction) {
    unsigned int temp;
    unsigned int tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCSizeT(size_t numInHeap, size_t *a, size_t *idxListHS, const bool direction) {
    size_t temp;
    size_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCPtrDiffT(size_t numInHeap, ptrdiff_t *a, size_t *idxListHS, const bool direction) {
    ptrdiff_t temp;
    ptrdiff_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCChar(size_t numInHeap, char *a, size_t *idxListHS, const bool direction) {
    char temp;
    char tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUChar(size_t numInHeap, unsigned char *a, size_t *idxListHS, const bool direction) {
    unsigned char temp;
    unsigned char tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCInt8T(size_t numInHeap, int8_t *a, size_t *idxListHS, const bool direction) {
    int8_t temp;
    int8_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUInt8T(size_t numInHeap, uint8_t *a, size_t *idxListHS, const bool direction) {
    uint8_t temp;
    uint8_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCInt16T(size_t numInHeap, int16_t *a, size_t *idxListHS, const bool direction) {
    int16_t temp;
    int16_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUInt16T(size_t numInHeap, uint16_t *a, size_t *idxListHS, const bool direction) {
    uint16_t temp;
    uint16_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCInt32T(size_t numInHeap, int32_t *a, size_t *idxListHS, const bool direction) {
    int32_t temp;
    int32_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUInt32T(size_t numInHeap, uint32_t *a, size_t *idxListHS, const bool direction) {
    uint32_t temp;
    uint32_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCInt64T(size_t numInHeap, int64_t *a, size_t *idxListHS, const bool direction) {
    int64_t temp;
    int64_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

void heapSortVecIdxCUInt64T(size_t numInHeap, uint64_t *a, size_t *idxListHS, const bool direction) {
    uint64_t temp;
    uint64_t tempP;
    heapSortVecIdxCMacro(numInHeap,a,idxListHS,direction);
}

//////////////////////////

void heapSortVecCDouble(size_t numInHeap, double *a, const bool direction) {
    double temp;
    double tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCFloat(size_t numInHeap, float *a, const bool direction) {
    float temp;
    float tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCInt(size_t numInHeap, int *a, const bool direction) {
    int temp;
    int tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUInt(size_t numInHeap, unsigned int *a, const bool direction) {
    unsigned int temp;
    unsigned int tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCSizeT(size_t numInHeap, size_t *a, const bool direction) {
    size_t temp;
    size_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCPtrDiffT(size_t numInHeap, ptrdiff_t *a, const bool direction) {
    ptrdiff_t temp;
    ptrdiff_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCChar(size_t numInHeap, char *a, const bool direction) {
    char temp;
    char tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUChar(size_t numInHeap, unsigned char *a, const bool direction) {
    unsigned char temp;
    unsigned char tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCInt8T(size_t numInHeap, int8_t *a, const bool direction) {
    int8_t temp;
    int8_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUInt8T(size_t numInHeap, uint8_t *a, const bool direction) {
    uint8_t temp;
    uint8_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCInt16T(size_t numInHeap, int16_t *a, const bool direction) {
    int16_t temp;
    int16_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUInt16T(size_t numInHeap, uint16_t *a, const bool direction) {
    uint16_t temp;
    uint16_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCInt32T(size_t numInHeap, int32_t *a, const bool direction) {
    int32_t temp;
    int32_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUInt32T(size_t numInHeap, uint32_t *a, const bool direction) {
    uint32_t temp;
    uint32_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCInt64T(size_t numInHeap, int64_t *a, const bool direction) {
    int64_t temp;
    int64_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
}

void heapSortVecCUInt64T(size_t numInHeap, uint64_t *a, const bool direction) {
    uint64_t temp;
    uint64_t tempP;
    heapSortVecCMacro(numInHeap,a,direction);
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
