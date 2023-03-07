/*A C++-only implementations of a functions for computing the gradient of
 *bistatic range.  See the Matlab equivalents for more comments.
 *
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
#include "CoordFuncs.hpp"
#include <math.h>

void rangeGradientCPP(const size_t numRows,double *J,double *tempSpace,const double *point,const bool useHalfRange,const double *lTx,const double *lRx) {
    size_t i;
    double normVal;
    //This keeps track of whether or not eh transmitter is collocated with
    //the target.
    bool transmitterNotWithTarget=false;

    //deltaTx=x-lTx;
    for(i=0;i<numRows;i++) {
        tempSpace[i]=point[i]-lTx[i];

        transmitterNotWithTarget=transmitterNotWithTarget||tempSpace[i];
    }

    normVal=0;
    for(i=0;i<numRows;i++) {
        normVal=normVal+tempSpace[i]*tempSpace[i];
    }  
    normVal=sqrt(normVal);

    //deltaTx.'/norm(deltaTx)
    if(transmitterNotWithTarget) {
        for(i=0;i<numRows;i++) {
            J[i]=tempSpace[i]/normVal;
        }
    } else {
        for(i=0;i<numRows;i++) {
            J[i]=0.0;
        }
    }

    //deltaRx=x-lRx;
    for(i=0;i<numRows;i++) {
        tempSpace[i]=point[i]-lRx[i];
    }

    normVal=0;
    for(i=0;i<numRows;i++) {
        normVal=normVal+tempSpace[i]*tempSpace[i];
    }  
    normVal=sqrt(normVal);

    for(i=0;i<numRows;i++) {
        J[i]=J[i]+tempSpace[i]/normVal;
    }

    if(useHalfRange) {
        for(i=0;i<numRows;i++) {
            J[i]=J[i]/2;
        }
    }
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
