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
#include "ClusterSetCPP.hpp"

size_t findFirstMaxCPP(const double *arr, const size_t arrayLen);

void spherHarmonicEvalCPP(double *V, double *gradV, const ClusterSetCPP<double> &C,const ClusterSetCPP<double> &S, const double *point, const size_t numPoints,const double a, const double c, const double scalFactor);
bool spherHarmonicCovCPP(double *sigma2, double *Sigma, const ClusterSetCPP<double> &CStdDev,const ClusterSetCPP<double> &SStdDev, const double *point, const size_t numPoints,const double a, const double c, const double scalFactor);

void NALegendreCosRatCPP(ClusterSetCPP<double> &PBarUVals, const double theta, const double scalFactor);
void NALegendreCosRatDerivCPP(ClusterSetCPP<double> &dPBarUValsdTheta, const ClusterSetCPP<double> &PBarUVals, const double theta);
void NALegendreCosRatDeriv2CPP(ClusterSetCPP<double> &d2PBarUValsdTheta2, const ClusterSetCPP<double> &dPBarUValsdTheta, const ClusterSetCPP<double> &PBarUVals, const double theta);

void normHelmHoltzCPP(ClusterSetCPP<double> &HBar,const double u, const double scalFactor);
void normHelmHoltzDerivCPP(ClusterSetCPP<double> &dHBardu,const ClusterSetCPP<double> &HBar);
void normHelmHoltzDeriv2CPP(ClusterSetCPP<double> &d2HBardu2,const ClusterSetCPP<double> &HBar);

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
