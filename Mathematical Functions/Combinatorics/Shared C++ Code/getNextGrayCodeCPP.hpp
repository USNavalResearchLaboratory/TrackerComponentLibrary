/*GETNEXTGRAYCODECPP A C++ function that returns the next gray code value.
*                 A gray code is a  binary code such that only one entry
*                 changes from 0 to 1 or back each step. This function can
*                 be used to get all subsets of an n-set.
*
*INPUTS:  n    The length of the code (const size_t).
*         code A pointer to a length-n array holding the previous code that
*              should be transformed to the next code value. This
*              is modified on return. It is a template, so different
*              types can be passed.
*        nCard The number of ones in code, implicitely passed. This value
*              is modified on return (size_t&)
*            j A value, implicitely passed, that will be updated to hold
*              the index of the entry in code that was modified.
*
*RETURN VALUE: The boolean return value is 1 if the final combination was
*              returned and 0 otherwise.
*
*The algorithm is NEXSUB taken from Chapter 1 of
*A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
*and Calculators, 2nd ed. New York: Academic press, 1978.
*Gray codes are also discussed in Chapter 7.2.1.1 of
*D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 2:
*Generating all Tuples and Permutations, Upper Saddle River, NJ:
*Addison-Wesley, 2009.
*
*The function checks against n so that if an invalid value of nCard is
*passed, it will not read/ write past the end of the array, even though it
*will return a code that is not the next in the sequence.
*
*Because this is a template function, the entire function has to be
*defined in the header file.
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef NEXTGRAYCODECPP
#define NEXTGRAYCODECPP

//Defines the size_t type
#include <stddef.h>
template <typename T>

bool getNextGrayCodeCPP(const size_t n, T* code, size_t &nCard, size_t &j) {    
    j=0;
    if(nCard%2!=0) {
        do {
            j++;
        } while(code[j-1]==0&&j<(n-1));
    }
    //Equivalent to code[j]=!code[j]; We turns 0 to 1 and 1 to 0.
    code[j]=1-code[j];
    nCard=nCard+2*static_cast<size_t>(code[j])-1;
    return nCard==static_cast<size_t>(code[n-1]);
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


