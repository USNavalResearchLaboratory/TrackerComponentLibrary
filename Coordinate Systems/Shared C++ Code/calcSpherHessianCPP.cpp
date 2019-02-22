/**CALCSPHERHESSIANCPP A C++-only implementation of a function for
 *          computing the Hessian of a monostatic or bistatic spherical
 *          measurement with respect to 3D Cartesian position.  See the
 *          Matlab equivalent for more comments.
 *
*July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>

void calcSpherHessianCPP(double *H,const double *x,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M) {
    rangeHessianGenCPP(3,H,x,useHalfRange,lTx,lRx);
    spherAngHessianGenCPP(H+9,x,systemType,lRx,M);
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
