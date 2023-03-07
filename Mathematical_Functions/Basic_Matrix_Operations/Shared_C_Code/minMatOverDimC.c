/*MINMATOVERDIMC A set of C-only implementations of a function to obtain
*                the hypermatrix that is result of minimizing a larger
*                hypermatrix over a particular dimension. This is
*                essentially min(C,[],minIdx) in matlab. The bulk of the
*                code for the function is written as a macro that is filled
*                in with the specific data type of the hypermatrix being
*                minimized. See the Matlab function minMatOverDim for
*                details.
*
*Each of the functions, such as heapSortVecCIdxSizeT, has the same
*inputs. The difference is only in the type of the inputs M and C, which
*are specified by the name of the function. The inputs to the functions
*have the following meaning:
* S The number of dimensions of matrix C and the number of elements in
*   nDims. S>1.
* nVals A length-S array where nVals[i] is the size of the ith dimension
*   of the hypermatrix C.
* M The matrix M into which the results of the minimization are placed.
*   This cannot be the same as C. The number of elements is the product
*   of everything in nVals except nVals[minIdx]. Values are stored by
*   column, not row. If minIdx>=S, then C will just be copied into M.
* C The original hypermatrix on which the minimization will be performed.
*   Values are stored by column, as in Fortran and Matlab.
* minIdx The dimension over which the minimization will be performed.
*    minIdx.=0.
*March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "basicMatOps.h"

#define minMatOverDimCMacro(S,nVals,M,C,minIdx) {\
    const size_t nMinIdx=nVals[minIdx];\
    size_t totalNumElsM, incrMinIdx, incrBigStep,CStartIdx, curMCol, curEl;\
    size_t i;\
\
    if(minIdx>=S) {\
        totalNumElsM=prodVectorSizeT(nVals,S);\
\
        for(i=0;i<totalNumElsM;i++) {\
            M[i]=C[i];\
        }\
        return;\
    }\
\
    incrMinIdx=prodVectorSizeT(nVals,minIdx);\
    incrBigStep=incrMinIdx*nVals[minIdx];\
\
    totalNumElsM=incrMinIdx*prodVectorSizeT(nVals+minIdx+1,S-minIdx-1);\
\
    CStartIdx=0;\
    curMCol=0;\
    for(curEl=0;curEl<totalNumElsM;curEl++) {\
        size_t curIdx=CStartIdx;\
\
        minVal=C[curIdx];\
        curIdx+=incrMinIdx;\
        for(i=1;i<nMinIdx;i++) {\
            curVal=C[curIdx];\
\
            if(curVal<minVal) {\
                minVal=curVal;\
            }\
            curIdx+=incrMinIdx;\
        }\
        M[curEl]=minVal;\
\
        if(((curEl+1)%incrMinIdx)==0) {\
            curMCol++;\
\
            CStartIdx=incrBigStep*curMCol;\
        } else {\
            CStartIdx++;\
        }\
    }\
}

void minMatOverDimCDouble(const size_t S,const size_t *nVals,double *  M,const double *C,const size_t minIdx) {
    double minVal;
    double curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCFloat(const size_t S,const size_t *nVals,float *  M,const float *C,const size_t minIdx) {
    float minVal;
    float curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCInt(const size_t S,const size_t *nVals,int *  M,const int *C,const size_t minIdx) {
    int minVal;
    int curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUInt(const size_t S,const size_t *nVals,unsigned int *  M,const unsigned int *C,const size_t minIdx) {
    unsigned int minVal;
    unsigned int curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCSizeT(const size_t S,const size_t *nVals,size_t *  M,const size_t *C,const size_t minIdx) {
    size_t minVal;
    size_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCPtrDiffT(const size_t S,const size_t *nVals,ptrdiff_t *  M,const ptrdiff_t *C,const size_t minIdx) {
    ptrdiff_t minVal;
    ptrdiff_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCChar(const size_t S,const size_t *nVals,char *  M,const char *C,const size_t minIdx) {
    char minVal;
    char curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUChar(const size_t S,const size_t *nVals,unsigned char *  M,const unsigned char *C,const size_t minIdx) {
    unsigned char minVal;
    unsigned char curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCInt8T(const size_t S,const size_t *nVals,int8_t *  M,const int8_t *C,const size_t minIdx) {
    int8_t minVal;
    int8_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUInt8T(const size_t S,const size_t *nVals,uint8_t *  M,const uint8_t *C,const size_t minIdx) {
    uint8_t minVal;
    uint8_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCInt16T(const size_t S,const size_t *nVals,int16_t *  M,const int16_t *C,const size_t minIdx) {
    int16_t minVal;
    int16_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUInt16T(const size_t S,const size_t *nVals,uint16_t *  M,const uint16_t *C,const size_t minIdx) {
    uint16_t minVal;
    uint16_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCInt32T(const size_t S,const size_t *nVals,int32_t *  M,const int32_t *C,const size_t minIdx) {
    int32_t minVal;
    int32_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUInt32T(const size_t S,const size_t *nVals,uint32_t *  M,const uint32_t *C,const size_t minIdx) {
    uint32_t minVal;
    uint32_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCInt64T(const size_t S,const size_t *nVals,int64_t *  M,const int64_t *C,const size_t minIdx) {
    int64_t minVal;
    int64_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
}

void minMatOverDimCUInt64T(const size_t S,const size_t *nVals,uint64_t *  M,const uint64_t *C,const size_t minIdx) {
    uint64_t minVal;
    uint64_t curVal;
    
    minMatOverDimCMacro(S,nVals,M,C,minIdx);
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
