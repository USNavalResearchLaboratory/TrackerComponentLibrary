/**INVERSEPERMUTATIONCPP C++ functions for computing inverse permutations.
*    The functions inversePermutationI and inversePermutationJ are in-place
*    and the function inversePermutationExpl is not in place. An option is
*    provided as to whether the values to be permuted are 1:N or 0:(N-1).
*
*The input xPerm is the original length N permutation. If startAt1 is true
*then the numbering in xPerm starts at 1; otherwise it starts at 0.
*xPerm is modified in the in-place algorithms to hold the inverse
*permutation. Otherwise, xInvPerm is the inverse permutation.
*
*Algorithms I and J are from Chapter 1.3.3 of [1].
*
*Note that xPerm is typed ptrdiff_t. If one has a size_t array, one can use
*reinterpret_cast to pass it without making a copy. Both size_t and
*ptrdiff_t are the same size.
*
*REFERENCES:
*[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
*    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
*
*October 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "combinatorialFuns.hpp"

void inversePermutationI(ptrdiff_t *xPerm,const size_t N,const bool startAt1) {
/*Algorithm I in [1]. An in-place algorithm. It relies on the indexation
 *starting from 1, not 0, so if it starts at 0, we must increment
 *everything.*/

    //Increment to make the indexation start at 1.
    if(startAt1==false) {
        for(size_t k=0;k<N;k++) {
            xPerm[k]++;
        }
    }
    
    //Step I1
    ptrdiff_t m=static_cast<ptrdiff_t>(N);
    ptrdiff_t j=-1;
    
    while(m>0) {
        //Step I2
        ptrdiff_t i=xPerm[m-1];
        if(i>=0) {
            while(1) {
                //Step I3
                xPerm[m-1]=j;
                j=-m;
                m=i;
                i=xPerm[m-1];
    
                //Step I4
                if(i>0) {
                    continue;
                }
    
                i=j;
                break;
            }
        }
    
        //Step I5
        xPerm[m-1]=-i;
        m--;
    }

    //Decrement everything so indexation starts at 0. 
    if(startAt1==false) {
        for(size_t k=0;k<N;k++) {
            xPerm[k]--;
        }
    }
}

void inversePermutationJ(ptrdiff_t *xPerm,const size_t N,const bool startAt1) {
/*Algorithm J in [1]. An in-place algorithm. It relies on the indexation
 *starting from 1, not 0, so if it starts at 0, we must increment
 *everything.*/

    //Increment to make the indexation start at 1.
    if(startAt1==false) {
        for(size_t k=0;k<N;k++) {
            xPerm[k]++;
        }
    }

    //Step J1;
    for(size_t k=0;k<N;k++) {
        xPerm[k]*=-1;
    }

    ptrdiff_t m=static_cast<ptrdiff_t>(N);
    while(m>0) {
        //Step J2
        ptrdiff_t j=m;

        //Step J3
        ptrdiff_t i=xPerm[j-1];
        while(i>0) {
            j=i;
            i=xPerm[j-1];
        }

        //Step J4
        xPerm[j-1]=xPerm[-i-1];
        xPerm[-i-1]=m;

        //Step J5
        m--;
    }

    //Decrement everything so indexation starts at 0. 
    if(startAt1==false) {
        for(size_t k=0;k<N;k++) {
            xPerm[k]--;
        }
    }
}
  
void inversePermutationExpl(ptrdiff_t *xInvPerm, const ptrdiff_t *xPerm,const size_t N,const bool startAt1) {
//Explicitly create the inverse permutation using a non-in-place algorithm.

if(startAt1) {
    for(ptrdiff_t k=0;k<static_cast<ptrdiff_t>(N);k++) {
        xInvPerm[xPerm[k]-1]=k+1;
    }
} else {
    for(ptrdiff_t k=0;k<static_cast<ptrdiff_t>(N);k++) {
        xInvPerm[xPerm[k]]=k;
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
