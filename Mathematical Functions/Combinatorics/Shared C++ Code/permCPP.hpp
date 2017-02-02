/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
#ifndef PERMCPP
#define PERMCPP

//Defines the size_t type
#include <stddef.h>

double permSquareCPP(const double *A, const size_t n, double *buffer);
double permCPP(const double *A, const size_t numRow, const size_t numCol, size_t *buffer);
double permCPPSkip(const double *A, const size_t numRowsTotal,const size_t * rows2Keep, const size_t *cols2Keep, const size_t numRowsKept, const size_t numColsKept,size_t *buffer);
double SigmaSSkip(const double *A,size_t *curComb,const size_t *rows2Keep, const size_t *cols2Keep,const size_t r,const size_t numRowsTotal,const size_t numRowsKept,const size_t numColsKept);
double SigmaS(const double *A,size_t *curComb,const size_t r,const size_t numRow,const size_t numCol);

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
