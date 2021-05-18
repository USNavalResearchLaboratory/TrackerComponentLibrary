/**BASICMATOPS Prototypes for C Implementations of basic operations with
 *             matrices. See the comments in the C and Matlab
 *             implementations for descriptions of the functionality.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef BASICMATOPS
#define BASICMATOPS

#ifdef __cplusplus
extern "C"
{
#endif

//Defines the size_t and ptrdiff_t types
#include <stddef.h>

//This defines the standard integer data types.
#include <stdint.h>

//The following functions are implemented in permuteMatrix.c.
size_t permuteMatrixCBufferSize(const size_t S);
void permuteMatrixC(const size_t S,const size_t *nVals, size_t *  nValsNew, void* CPermMat, const void *C, const size_t CElSize, void *tempBuffer, const size_t *order);
void permute2DimsC(const size_t *nDims,size_t *  nDimsNew, void * CPermMat,const void * COrigMat,const size_t CElSize,const size_t *dimsOrder);
void permute3DimsC(const size_t *nDims,size_t *  nValsNew,void * CPermMat,const void * COrigMat, const size_t CElSize,const size_t *dimsOrder);

//The following functions are implemented in minMatOverDimC.c
void minMatOverDimCDouble(const size_t S,const size_t *nVals,double *  M,const double *C,const size_t minIdx);
void minMatOverDimCFloat(const size_t S,const size_t *nVals,float *  M,const float *C,const size_t minIdx);
void minMatOverDimCInt(const size_t S,const size_t *nVals,int *  M,const int *C,const size_t minIdx);
void minMatOverDimCUInt(const size_t S,const size_t *nVals,unsigned int *  M,const unsigned int *C,const size_t minIdx);
void minMatOverDimCSizeT(const size_t S,const size_t *nVals,size_t *  M,const size_t *C,const size_t minIdx);
void minMatOverDimCPtrDiffT(const size_t S,const size_t *nVals,ptrdiff_t *  M,const ptrdiff_t *C,const size_t minIdx);
void minMatOverDimCChar(const size_t S,const size_t *nVals,char *  M,const char *C,const size_t minIdx);
void minMatOverDimCUChar(const size_t S,const size_t *nVals,unsigned char *  M,const unsigned char *C,const size_t minIdx);
void minMatOverDimCInt8T(const size_t S,const size_t *nVals,int8_t *  M,const int8_t *C,const size_t minIdx);
void minMatOverDimCUInt8T(const size_t S,const size_t *nVals,uint8_t *  M,const uint8_t *C,const size_t minIdx);
void minMatOverDimCInt16T(const size_t S,const size_t *nVals,int16_t *  M,const int16_t *C,const size_t minIdx);
void minMatOverDimCUInt16T(const size_t S,const size_t *nVals,uint16_t *  M,const uint16_t *C,const size_t minIdx);
void minMatOverDimCInt32T(const size_t S,const size_t *nVals,int32_t *  M,const int32_t *C,const size_t minIdx);
void minMatOverDimCUInt32T(const size_t S,const size_t *nVals,uint32_t *  M,const uint32_t *C,const size_t minIdx);
void minMatOverDimCInt64T(const size_t S,const size_t *nVals,int64_t *  M,const int64_t *C,const size_t minIdx);
void minMatOverDimCUInt64T(const size_t S,const size_t *nVals,uint64_t *  M,const uint64_t *C,const size_t minIdx);

//The following functions are implemented in basicMatOps.c and performing
//basic matrix operations that are simple to perform in Matlab, but that
//can be tedious when programming in C. There are no Matlab interfaces for
//these functions.
void permuteRowsPtrDiffT(const size_t n1,const size_t n2,ptrdiff_t *  XPerm,const ptrdiff_t *  XOrig,const size_t *rowOrder);
void identMatD(const size_t numDims, double *  I);
double vecMinD(const double*vec, const size_t numEls);
double sumVectorD(const double *vector, const size_t numEls);
size_t sumVectorSizeT(const size_t *vector, const size_t numEls);
size_t prodVectorSizeT(const size_t *vector, const size_t numEls);
size_t accessHyperMatSizeT(const size_t *prodTerms, const size_t *C,const size_t *idx,const size_t S);
void addVec2MatEndMin(double *d,size_t *dIdx,const double *C,const double *u,const size_t numElsInD, const size_t numElsInU);
void fixedIdxMatSub(const size_t S, const size_t *nVals,const double *M, double *  C,const size_t dim, const size_t idx);
size_t sub2Idx2D(const size_t idx1,const size_t idx2,const size_t dim1);
size_t sub2Idx3D(const size_t idx1,const size_t idx2,const size_t idx3,const size_t dim1,const size_t dim2);

#ifdef __cplusplus
} // extern "C"
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
