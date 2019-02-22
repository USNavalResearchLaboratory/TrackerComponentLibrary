/**MATHFUNCS  A header file for C++ implementations of mathematical
 *            functions. See the files implementing each function for more
 *            details on their usage.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MATHFUNCSCPP
#define MATHFUNCSCPP

#include <stddef.h>
#include "CountingClusterSetCPP.hpp"
#include <complex>

size_t findFirstMaxCPP(const double *arr, const size_t arrayLen);

void spherHarmonicEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetCPP<double> &C,const CountingClusterSetCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicEvalCPPComplex(double *VReal, double *VImag, double *gradVReal, double *gradVImag, double *HessianVReal, double *HessianVImag,const CountingClusterSetCPP<double> &CReal,const CountingClusterSetCPP<double> &CImag,const CountingClusterSetCPP<double> &SReal,const CountingClusterSetCPP<double> &SImag, const double *point, const size_t numPoints, const std::complex <double> a, const std::complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicSetEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetVecCPP<double> &C,const CountingClusterSetVecCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicSetEvalCPPComplex(double *VReal, double *VImag, double *gradVReal, double *gradVImag, double *HessianVReal, double *HessianVImag,const CountingClusterSetVecCPP<double> &CReal,const CountingClusterSetVecCPP<double> &CImag,const CountingClusterSetVecCPP<double> &SReal,const CountingClusterSetVecCPP<double> &SImag, const double *point, const size_t numPoints, const std::complex <double> a, const std::complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
bool spherHarmonicCovCPP(double *sigma2, double *Sigma, const CountingClusterSetCPP<double> &CStdDev,const CountingClusterSetCPP<double> &SStdDev, const double *point, const size_t numPoints,const double a, const double c, const double scalFactor);

void NALegendreCosRatCPP(CountingClusterSetCPP<double> &PBarUVals, const double theta, const double scalFactor);
void NALegendreCosRatDerivCPP(CountingClusterSetCPP<double> &dPBarUValsdTheta, const CountingClusterSetCPP<double> &PBarUVals, const double theta);
void NALegendreCosRatDeriv2CPP(CountingClusterSetCPP<double> &d2PBarUValsdTheta2, const CountingClusterSetCPP<double> &dPBarUValsdTheta, const CountingClusterSetCPP<double> &PBarUVals, const double theta);

void normHelmHoltzCPP(CountingClusterSetCPP<double> &HBar,const double u, const double scalFactor);
void normHelmHoltzDerivCPP(CountingClusterSetCPP<double> &dHBardu,const CountingClusterSetCPP<double> &HBar);
void normHelmHoltzDeriv2CPP(CountingClusterSetCPP<double> &d2HBardu2,const CountingClusterSetCPP<double> &HBar);
void normHelmHoltzDeriv3CPP(CountingClusterSetCPP<double> &d3HBardu3,const CountingClusterSetCPP<double> &HBar);

double wrapRangeCPP(const double val2Wrap, const double minBound, const double maxBound);
double wrapRangeMirrorCPP(const double val2Wrap, const double minBound, const double maxBound);
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
