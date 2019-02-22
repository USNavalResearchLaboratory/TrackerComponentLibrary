/**GETNEXTCOMBO C++ functions that returns the next combination in
*              lexicographic order given the current combination. If the
*              final combination in the sequence is passed, then the
*              return value is 1.     
*
*The function getNextComboCPP implements the algorithm whereby the values
*of the combination start at 0 and the function getNextComboCPPFromOne
*implement the algorithm for values of the function starting at 1.
*
*The algorithm is from
*C. J. Mifsud, "Algorithm 154: Combination in lexicographical order," 
*Communications of the ACM, vol. 6, no. 3 pp. 103, Mar. 1963.
*
*More details on the algorithm are in the Matlab implementation.
*
*December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "getNextComboCPP.hpp"

bool getNextComboCPP(size_t *I,const size_t n,const size_t r) {
    if(I[r-1]<n-1) {
        I[r-1]+=1;
        return false;
    } else {
        size_t j, s;
        
        for(j=r-1;j>0;j--) {
           if(I[j-1]<n+j-r-1) {
               I[j-1]=I[j-1]+1;
        
               for(s=j;s<r;s++) {
                  I[s]=I[j-1]+s-(j-1); 
               }
               return false;
           }
        }
                   
        return true;
    }
}

bool getNextComboCPPFromOne(size_t *I,const size_t n,const size_t r) {
    if(I[r-1]<n) {
        I[r-1]+=1;
        return false;
    } else {
        size_t j, s;
        
        for(j=r-1;j>0;j--) {
           if(I[j-1]<n+j-r) {
               I[j-1]=I[j-1]+1;
        
               for(s=j;s<r;s++) {
                  I[s]=I[j-1]+s-(j-1); 
               }
               return false;
           }
        }
                   
        return true;
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
