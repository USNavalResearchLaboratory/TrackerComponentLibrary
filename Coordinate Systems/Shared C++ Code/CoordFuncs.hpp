/**COORDFUNCS A header file for C++ implementations of coordinate
 *            conversion functions as well as Jacobians and Hessians
 *            related to coordinate conversions. See the files implementing
 *            each function for more details on their usage. Most of the
 *            documentation for the algorithms is with the corresponding
 *            Matlab implementations.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef COORDFUNCS
#define COORDFUNCS
#include <stddef.h>

void spher2CartCPP(double *cartPoint,const double *point,size_t systemType);
void spher2CartGenCPP(double *retData,const double *point,size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M);
void Cart2SphereGenCPP(double *retData,const double *cartPoints,const size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M);
void spher2CartNoRangeCPP(double *retData,const double *point,const size_t systemType,const double *M);
void getENUAxesCPP(double *u, double *c,const double *plhPoint,const bool justVertical,const double a,const double f);
void getEllipsHarmAxesCPP(double *u, double *c,const double *pointHarmon,const double E);
void Cart2EllipsHarmonCPP(double *pointsHarmon,const double *cartPoints, const size_t numPoints, const double E);

void ruv2CartGenCPP(double *retData,const double *z,const bool useHalfRange,const double *zTx,const double *zRx,const double *M, bool hasW);
void Cart2RuvGenCPP(double *retData,const double *points,bool useHalfRange,double *zTx,double *zRx,double *M,bool includeW);

double getRangeRate2DCPP(const double *points,bool useHalfRange,const double *xTx,const double *xRx);
double getRangeRate3DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx);

void rangeGradientCPP(const size_t numRows,double *J,double *tempSpace,const double *point,const bool useHalfRange,const double *lTx,const double *lRx);
void spherAngGradientGenCPP(double *retMat,const double *xG,const size_t systemType,const double *lRx, const double *M);
void calcSpherConvJacobGenCPP(double *J,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherConvJacobCPP(double *J, const double *point, const size_t systemType);
void calcSpherInvJacobCPP(double *J,const double *z,const size_t systemType);

void rangeHessianCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange);
void rangeHessianGenCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange,const double *lTx,const double *lRx);
void spherAngHessianCPP(double *HTotal,const double *xG,const size_t systemType);
void spherAngHessianGenCPP(double *HTotal,const double *xG,const size_t systemType,const double *lRx,const double *M);
void calcSpherInvHessianCPP(double *HTotal,const double *z,const size_t systemType);
void spherAngHessianGenCPP(double *H,const double *x,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherHessianCPP(double *H,const double *x,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherConvHessianCPP(double *H, const double *point, const size_t systemType, const bool useHalfRange);
void calcSpherConvHessianGenCPP(double *H,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);

#endif

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
