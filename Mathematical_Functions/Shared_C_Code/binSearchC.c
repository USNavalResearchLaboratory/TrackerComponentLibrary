/**BINSEARCHC A C implementation of an algorithm to perform a binary
*             search.
*
*INPUTS: numInVec The number of elements in the vector vec.
*              vec A vector with elements sorted in increasing order.
*              key The value that one wishes to find in the vector vec.
*           choice An optional parameter that determines what is returned
*                  if key is not found.
*                  0 means return the closest value.
*                  1 means return the next lower value if there is one,
*                    otherwise return the lowest value in vec.
*                  2 means return the next higher value if there is one,
*                    otherwise return the highest value in vec.
*
*OUTPUTS: The return value is the index of the found element.
*
*December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncsC.h"
//Needed for fabs
#include <math.h>

size_t binSearchC(size_t numInVec, double *vec, double key,int choice) {
    size_t UB, LB,mid;
    double diff1,diff2;
    
    UB=numInVec-1;
    LB=0;
    
    if(vec[UB]<=key){
        return UB;
    } else if(vec[LB]>=key) {
        return LB;
    }

    while(UB!=LB) {
        mid=(UB+LB)/2;

        //If the search has reached the deepest level. 
        if(mid==UB||mid==LB) {
            //First, check to see whether either one equals the key.
            if(vec[LB]==key) {
                return LB;
            } else if(vec[UB]==key) {
                return UB;
            }

            //If the key is not in vec, then the return value depends on 
            //the choice parameter.
            switch(choice) {
                case 1:
                   return LB;
                case 2:
                   return UB;
                default://Return the closest value.
                   diff1=fabs(key-vec[LB]);
                   diff2=fabs(key-vec[UB]);
                   if(diff1<diff2) {
                       return LB;
                   } else {
                       return UB;
                   }
            }
        }              

        if(vec[mid]>key) {
            UB=mid;
            continue;
        } else if(vec[mid]<key) {
            LB=mid;
            continue;
        } else {
            return mid;
        }
    }
    return UB;
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
