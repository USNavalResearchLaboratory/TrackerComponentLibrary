/**CALCSPHERCONVHESSIANCPP These are C++-only implementations of functions
 *          for computing the Hessian of  a monostatic or bistatic
 *          spherical measurement with respect to 3D Cartesian position
 *          when given a spherical measurement. See the Matlab
 *          implementation for more comments.
 *
*July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>

void calcSpherConvHessianGenCPP(double *H,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M) {
    double zCart[3];

    spher2CartGenCPP(zCart,zSpher,systemType,useHalfRange,lTx,lRx,M);

    rangeHessianGenCPP(3,H,zCart,useHalfRange,lTx,lRx);
    spherAngHessianGenCPP(H+9,zCart,systemType,lRx,M);
}

void calcSpherConvHessianCPP(double *H, const double *point, const size_t systemType, const bool useHalfRange) {
    double zCart[3];

    spher2CartCPP(zCart,point,systemType);
    rangeHessianCPP(3,H,zCart,useHalfRange);
    spherAngHessianCPP(H+9,zCart,systemType);
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
