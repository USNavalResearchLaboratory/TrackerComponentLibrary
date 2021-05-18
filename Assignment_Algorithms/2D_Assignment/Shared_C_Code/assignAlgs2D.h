/**ASSIGNALGS2D This is a header for C language functions to solve variants
 *              of the 2D assignment problem. The inputs of the specific
 *              functions are described in their implementation files,
 *              which are assign2DMissedDetectC.c, assign2DFullC.c, and
 *              assign2DC.c. Additionally, functions describing the amount
 *              of memory needed for the aforementioned functions are
 *              implemented and commented inline in this header as
 *              assign2DCBufferSize, assign2DMissedDetectCBufferSize,
 *              assign2DFullCBufferSize, and assign2DFullCAltBufferSize.
 * 
 *Better understanding of the algorithms can usually be obtained from
 *looking at the Matlab implementations.
 *
 *January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef ASSIGNALGS2D
#define ASSIGNALGS2D

#ifdef __cplusplus
extern "C"
{
#endif

//This is needed for the bool type to be defined in C.
//C99 has stdbool, earlier versions do not.
#ifndef __cplusplus
#if __STDC_VERSION__>=199901L
#include <stdbool.h>
#else
//This is for the bool that some functions might define
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

size_t assign2DSimpBufferSize(size_t numRow,size_t numCol);
bool assign2DSimp(const bool maximize, const double *  CIn, double *  gain, ptrdiff_t *  col4RowOrig, ptrdiff_t *  row4ColOrig, void *tempBuffer, double *uCols, double *vRows, size_t numRow, size_t numCol);

size_t assign2DCBufferSize(const size_t numRow, const size_t numCol);
bool assign2DC(const bool maximize, double * C, double *  gain, ptrdiff_t *  col4row, ptrdiff_t *  row4col, void *tempBuffer, double *  u, double *  v, const size_t numRow, const size_t numCol);
double assign2DCBasic(const double *C, ptrdiff_t *  col4row, ptrdiff_t *  row4col, void *tempBuffer, double * u, double *  v, const size_t numRow, const size_t numCol);

size_t assign2DMissedDetectCBufferSize(const size_t numRowsTrue, const size_t numCol);
bool assign2DMissedDetectC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  tuples, void *  tempBuffer, double *  u, double *  v, const size_t numRowsTrue, const size_t numCol);
double assign2DCMissedDetectBasic(const double *C, ptrdiff_t *  tuples, void *tempBuffer, double *  u, double *  v, const size_t numRowsTrue, const size_t numCol);

size_t assign2DFullCBufferSize(const size_t numRowsTrue, const size_t numColsTrue);
bool assign2DFullC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  tuples,size_t *  numTuples,void *  tempBuffer,double *  u,double *  v,const size_t numRowsTrue, const size_t numColsTrue);
double assign2DFullCBasic(const double *C,ptrdiff_t *  tuples, size_t *  numTuples, void *tempBuffer, double *  u, double *  v, const bool hasUnconstTuple, const double zeroOffset, const size_t numRowsTrue,const size_t numColsTrue);

size_t assign2DFullCAltBufferSize(const size_t numRow, const size_t numCol);
bool assign2DFullCAlt(const bool maximize, const double *C, double *  gain, ptrdiff_t *tuples, size_t *  numTuples,void *tempBuffer, double *u, double *v, size_t numRow, size_t numCol);

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
