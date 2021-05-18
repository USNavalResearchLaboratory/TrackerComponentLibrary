/*WRAPRANGECPP A direct C++ implementation of functions to wrap a value to
 *             a given range either circularly (wrapRangeCPP) or reflecting
 *             back off the end (wrapRangeMirrorCPP).
 *
 *This function is a C++ implementation of the function wrapRange in
 *Matlab. See the Matlab function wrapRange for more details on the
 *implementation.
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include <cmath>
#include "mathFuncs.hpp"

double wrapRangeCPP(const double val2Wrap, const double minBound, const double maxBound) {
    const double a=(val2Wrap-minBound);
    const double m=(maxBound-minBound);
    
    //Deal with the case where the value is in the primary interval.
    //Treating this as a special case gets rid of finite precision issues
    //that can make values very close to the upper bound of the interval
    //wrap to the bottom. 
    if(val2Wrap>=minBound&&val2Wrap<maxBound) {
        return val2Wrap;
    }

    return a-m*floor(a/m)+minBound;
}

double wrapRangeMirrorCPP(const double val2Wrap, const double minBound, const double maxBound) {
    double retVal,spreadVal=maxBound-minBound;
    const double pi=3.1415926535897932384626433832795028841971693993751;
     
    //First, shift the values so that the region of wrapping is centered
    //symmetrically around zero.
    retVal=val2Wrap-minBound-spreadVal/2.0;
    
    //Now, scale the values so that  the bounds reduce to (+/-)pi/2,
    //meaning that the spreadVal becomes pi.
    retVal=(pi/spreadVal)*retVal;
    
    //Now, take the tangent and two-quadrant inverse tangent to wrap
    //everything back where it should be. Note that Matlab properly
    //handles infinite values in the inverse tangent.
    retVal=asin(sin(retVal));
    
    //Scale everything back to the original size
    retVal=(spreadVal/pi)*retVal;
    
    //Shift the origin back to where it should be.
    return retVal+minBound+spreadVal/2.0;
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
