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

/**FINDFIRSTMAXCPP Given an array of sorted vector of values in increasing
*              order, which might be full of duplicates, find the first
*              occurrence of the maximum value.
*
*INPUTS: arr A length arrayLen array of values, with possible repeats,
*            sorted in increasing order.
* arrayLen The number of  elements in arr.
*
*OUTPUTS: The return valuer is index of the first occurrence of the maximum
*         value element in arr. 
*
*January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
size_t findFirstMaxCPP(const double *arr, const size_t arrayLen);


/**POLYVALGEN Return the value of a scalar polynomial evaluated at the point
*       or points given in x. The polynomial can be specified in the format
*       y=p(1)*x^N+p(2)*x^(N-1)+...+p(N)*x+p(N+1)
*       or in the format
*       y=p(N+1)*x^N+p(N)*x^(N-1)+...+p(2)*x+p(1)
*       as specified by the direction option. This function is similar to
*       polyval, but allows one to choose between the formats without
*       calling flip to reverse the order of the elements in p.
*
*INPUTS: p A length N array of polynomial coefficients.
*        x A single point at which the coefficients should be evaluated.
* firstTermHighestOrder An boolean parameter specifying which of the two
*          polynomial formats is used. Possible values are:
*          true p(1) is the coefficient of the highest order x term.
*          false p(1) is the coefficient of the lowest order x term.
*
*OUTPUTS: The return value is the value of the polynomial in p evaluated at
*         the point in x. This is the same size as x.
*
*The function is evaluated using Horner's rule, which is described in [1].
*
*Note that this template function must be called with all 3 types
*specified. The return value is the firtat, p is the second and x is the
*third. Having 3 types allows the function to be called, for example, with
*complex and real inputs, returning a complex input.
*
*As an example, to explicitly specify that doubles are inputted and
*returned, the function is called as:
*polyValGenCPP<double,double,double>(numP,p,x,firstTermHighestOrder)
*
*REFERENCES:
*[1] Weisstein, Eric W. "Horner's Rule." From MathWorld--A Wolfram Web
*    Resource. http://mathworld.wolfram.com/HornersRule.html
*
*January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
template <class T1,class T2,class T3>
T1 polyValGenCPP(const size_t numP, const T2 *p, const T3 x, const bool firstTermHighestOrder) {
    
    if(numP==0) {//Special case.
        return static_cast<T1>(0.0);
    }
    
    T1 y;
    
    if(firstTermHighestOrder==true) {
        y=p[0];
        
        for(size_t k=1;k<numP;k++) {
            y=x*y+p[k];   
        }
    } else {
        y=p[numP-1];
        
        for(ptrdiff_t k=numP-2;k>=0;k--) {
            y=x*y+p[k];
        }
    }
    
    return y;
}

void spherHarmonicEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetCPP<double> &C,const CountingClusterSetCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicEvalCPPComplex(std::complex<double> *VRet, std::complex<double> *gradVRet, std::complex<double> *HessianVRet,const CountingClusterSetCPP<std::complex<double>> &C,const CountingClusterSetCPP<std::complex<double>> &S, const double *point, const size_t numPoints, const std::complex <double> a, const std::complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicSetEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetVecCPP<double> &C,const CountingClusterSetVecCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
void spherHarmonicSetEvalCPPComplex(std::complex<double> *VRet, std::complex<double> *gradVRet, std::complex<double> *HessianVRet,const CountingClusterSetVecCPP<std::complex<double>> &C,const CountingClusterSetVecCPP<std::complex<double>> &S, const double *point, const size_t numPoints, const std::complex <double> a, const std::complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm);
bool spherHarmonicCovCPP(double *sigma2, double *Sigma, const CountingClusterSetCPP<double> &CStdDev,const CountingClusterSetCPP<double> &SStdDev, const double *point, const size_t numPoints,const double a, const double c, const double scalFactor);

void NALegendreCosRatCPP(CountingClusterSetCPP<double> &PBarUVals, const double theta, const double scalFactor);
void NALegendreCosRatDerivCPP(CountingClusterSetCPP<double> &dPBarUValsdTheta, const CountingClusterSetCPP<double> &PBarUVals, const double theta);
void NALegendreCosRatDeriv2CPP(CountingClusterSetCPP<double> &d2PBarUValsdTheta2, const CountingClusterSetCPP<double> &dPBarUValsdTheta, const CountingClusterSetCPP<double> &PBarUVals, const double theta);

void normHelmHoltzCPP(CountingClusterSetCPP<double> &HBar,const double u, const double scalFactor);
void normHelmHoltzDerivCPP(CountingClusterSetCPP<double> &dHBardu,const CountingClusterSetCPP<double> &HBar);
void normHelmHoltzDeriv2CPP(CountingClusterSetCPP<double> &d2HBardu2,const CountingClusterSetCPP<double> &HBar);
void normHelmHoltzDeriv3CPP(CountingClusterSetCPP<double> &d3HBardu3,const CountingClusterSetCPP<double> &HBar);

/**WRAPRANGECPP An algorithm to bound values to a specific range, wrapping
*               around if the value is outside the range. 
*
*INPUTS: val2Wrap The value to wrap to the range.
*        minBound The lower scalar bound of the output parameters.
*        maxBound A value > minBound that is the upper bound to the allowable
*                 range of output values.
*
*OUTPUTS: The return value is the wrapped value in the range:
*          minBound<=val<maxBound.
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
double wrapRangeCPP(const double val2Wrap, const double minBound, const double maxBound);

/**WRAPRANGECPP An algorithm to bound values to a specific range, mirroring
*               back if the values goes beyong the end of the range.  For
*               example, if minBound=-pi/2, and maxBound=pi/2, then a value
*               that is some epsilon above pi/2 will be mapped to a point
*               some epsilon below pi/2 and a value some eps below -pi/2
*               will be mapped to a point some epsilon above -pi/2.
*
*INPUTS: val2Wrap The value to wrap to the range.
*        minBound The lower scalar bound of the output parameters.
*        maxBound A value > minBound that is the upper bound to the allowable
*                 range of output values.
*
*OUTPUTS: The return value is the mirror-wrapped value in the range:
*          minBound<=val<maxBound.
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
double wrapRangeMirrorCPP(const double val2Wrap, const double minBound, const double maxBound);

/**EVALCLENSHAWRECURSERIES Use Clenshaw's method to evaluate a series of
*           the form sum_{k=0}^Nc(k+1)*F(k) where F(k) is a function that
*          is subject to the recurrence relation
*           F(k+1)=alpha*F(k)+beta*F(k-1)
*           Clenshaw's method is numerically stabler than directly using
*           the above recursion and explicitly evaluating the sum.
*           Examples of functions satifying this type of relation are sine,
*           cosine, Legendre polynomials, and Bessel functions, among
*           others.
*
*INPUTS: numC The number of elements in c.
*           c A length numC array.
*  alphaVal, betaVal These are the coefficients in the recursion on F given
*           above. Here, they are taken to be constant.
*    F0, F1 If method=0, then F0=F(0) and F1=F(1). The first two values of
*           F are needed to start the recursion. On the other hand, if
*           method=1, then F0=F(N) and F1=F(N-1) as the recursion goes in
*           the opposite direction.
* useReverse This optionally selects the type of series to use. Possible
*           values are:
*           false Use the Clenshaw series going backwards from N such that
*             F0=F(0) and
*             F1=F(1).
*           true Use the Clenshaw series going forward from 0 such that
*             F0=F(N) and F1=F(N-1). This method is generally only
*             beneficial if F(k) is small when k is large and c(k) are
*             small when k is small.
*
*OUTPUTS: The return value is the value of the series.
*
*January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
template <class T>
T evalClenshawRecurSeriesConst(const size_t numC,const T *c,const T alphaVal,const T betaVal,const T F0,const T F1,const bool useReverse) {
//alphaVal and betaVal are constants

    T y1=0;
    T y2=0;
    if(useReverse==false) {//The standard Clenshaw recursion.
        for(size_t k=numC-1;k>0;k--) {
            const T y=alphaVal*y1+betaVal*y2+c[k];
            y2=y1;
            y1=y;
        }
        return betaVal*F0*y2+F1*y1+F0*c[0];
    } else {//The Clenshaw recursion going in the opposite direction.
        for(size_t k=0;k<numC-1;k++) {
            const T y=(y2-alphaVal*y1-c[k])/betaVal;
            y2=y1;
            y1=y;
        }
        const T FN=F0;
        const T FN1=F1;

        return c[numC-1]*FN-betaVal*FN1*y1-FN*y2;
    }
}

/**EVALCLENSHAWRECURSERIES Use Clenshaw's method to evaluate a series of
*           the form sum_{k=0}^Nc(k+1)*F(k) where F(k) is a function that
*           is subject to the recurrence relation
*           F(k+1)=alpha(k)*F(k)+beta(k)*F(k-1)
*           Clenshaw's method is numerically stabler than directly using
*           the above recursion and explicitly evaluating the sum.
*           Examples of functions satifying this type of relation are sine,
*           cosine, Legendre polynomials, and Bessel functions, among
*           others.
*
*INPUTS: numC The number of elements in c.
*           c A length numC array.
*  alphaVal, betaVal These are the handles to functions that take k (as in
*             the equation above) as an input as well as the dataAlpha and
*             dataBeta inputs to this function. They return an appropriate
*             value of alpha and beta.
*    F0, F1 If method=0, then F0=F(0) and F1=F(1). The first two values of
*           F are needed to start the recursion. On the other hand, if
*           method=1, then F0=F(N) and F1=F(N-1) as the recursion goes in
*           the opposite direction.
* useReverse This optionally selects the type of series to use. Possible
*           values are:
*           false Use the Clenshaw series going backwards from N such that
*             F0=F(0) and
*             F1=F(1).
*           true Use the Clenshaw series going forward from 0 such that
*             F0=F(N) and F1=F(N-1). This method is generally only
*             beneficial if F(k) is small when k is large and c(k) are
*             small when k is small.
* dataAlpha, dataBeta These are passed to the alphaVal and betaVal
*             functions. These are only an option in case the functions
*             need additional data. if the function doesn't use these, then
*             one can just pass nullptr for these values.
*
*OUTPUTS: The return value is the value of the series.
*
*January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
template <class T>
T evalClenshawRecurSeriesFunc(const size_t numC,const T *c,T (*alphaVal)(size_t,void*),T (*betaVal)(size_t,void*),const T F0,const T F1,const bool useReverse, void *dataAlpha, void*dataBeta) {
//alphaVal and betaVal are vectors with the same length as c.
    
    T y1=0;
    T y2=0;
    if(useReverse==false) {//The standard Clenshaw recursion.

        for(size_t k=numC-1;k>0;k--) {
            const T y=alphaVal(k,dataAlpha)*y1+betaVal(k+1,dataBeta)*y2+c[k];
            y2=y1;
            y1=y;
        }
        return betaVal(1,dataBeta)*F0*y2+F1*y1+F0*c[0];
    } else {//The Clenshaw recursion going in the opposite direction.
        for(size_t k=0;k<numC-1;k++) {
            const T y=(y2-alphaVal(k,dataAlpha)*y1-c[k])/betaVal(k+1,dataBeta);
            y2=y1;
            y1=y;
        }
        const T FN=F0;
        const T FN1=F1;

        return c[numC-1]*FN-betaVal(numC-1,dataBeta)*FN1*y1-FN*y2;
    }   
}

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
