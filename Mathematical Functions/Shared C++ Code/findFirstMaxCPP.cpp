/*FINDFIRSTMAXCPP A C++ implementations of a function such that given an
 *                array of sorted vector of values in increasing
 *                order, which might be full of duplicates, it finds the 
 *                first occurrence of the maximum value.
 *
 *See the Matlab implementation findFirstMax.m for more comments on the
 *algorithm.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"

size_t findFirstMaxCPP(const double *arr,const size_t arrayLen) {
    double maxVal;
    size_t maxIdx, minIdx, midIdx;
    //If there is only one element, then the problem is trivial.
    if(arrayLen==1) {
        return 0;
    }

//Check for the simplest solution: The last element is the unique maximum.
    maxVal=arr[arrayLen-1];
    if(arr[arrayLen-2]!=maxVal){
        return arrayLen-1;
    }

//In this case, we must perform a binary search for the second highest
//element.
    maxIdx=arrayLen-2;
    minIdx=0;

    while(maxIdx!=minIdx&&maxIdx-minIdx>1) {
        midIdx=(maxIdx+minIdx)/2;
        if(arr[midIdx]==maxVal) {
           maxIdx=midIdx;
        } else {
           minIdx=midIdx;
        }
    }

    if(maxIdx==minIdx) {
        return maxIdx;
    } else {
        if(arr[minIdx]==maxVal) {
            return minIdx;
        } else {
            return maxIdx;
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
