/**INDEX2NDIMCPP C++ functions related to the indexation of
 *               multidimensional arrays.
 *
 *October 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Get rid of a useless warning in Visual Studio about inlining.
#ifdef _MSC_VER
#pragma warning( disable: 4711 )
#endif

#include <cstddef>  // for size_t
#include "miscFuns.hpp"

bool index2NDimCPP(const size_t numDim, const size_t numIdx,const size_t * const idx, const size_t * const maxVals, size_t *newIndices) {
/**INDEX2NDIMCPP Given the number of dimensions of a multidimensional array
 *           and the linear index of an item in the array (starting from 1,
 *           not 0), find the set of indices (starting from 1) that
 *           addresses the point in the array. This function is like the
 *           ind2sub function, but it returns a single vector when dealing
 *           with multidimensional arrays rather than returning multiple
 *           outputs.
 *
 *INPUTS: numDim The dimensionality of the indices.
 *        numIdx The number of indices provided that should be unranked to
 *               form tuples.
 *           idx A length numIdx array holding the indices to unrank. These
 *               start from 1, not 0.
 *       maxVals The product values related to the dimensions, as described
 *               above.
 *    newIndices A numDim*numIdx array into which the unranked tuples will
 *               be placed. These are placed one after the other.
 *
 *OUTPUTS: The return value indicates whether one is past the final indexed
 *         value. If idx is too large (or 0, which is invalid) then nothing
 *         is put in newIndices and this function returns true. Otherwise,
 *         the function returns false and newIndices holds the new set of
 *         indices.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t maxVal=maxVals[numDim];
    
    for(size_t i=0;i<numIdx;i++) {
        size_t curIndic=numDim;
        size_t curIdx=idx[i];
        
        //The index is beyond the final one or is invalid.
        if(curIdx>maxVal||curIdx==0) {
            return true;
        }
        
        do {
            curIndic--;
            //integer division.
            size_t wholeVal=(curIdx-1)/maxVals[curIndic];
            
            newIndices[i*numDim+curIndic]=wholeVal+1;
            curIdx-=wholeVal*maxVals[curIndic];
        } while(curIndic>0);
    }
    return false;
}

bool index2NDimByDimCPP(const size_t numDim, const size_t numIdx,const size_t *idx, const size_t *dims, size_t *newIndices, size_t *tempSpace) {
/**INDEX2NDIMBYDIMCPP This is the same as the index2NDimCPP function,
 *     except it is parameterized by the dimensions (dims) of the index
 *     vector (dims) rather than by the product vector that is used in
 *     index2NDimCPP. The tempSpace buffer must have enough space to hold
 *     at least numDim+1 size_t values.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    //This memory buffer is numDim+1 in size.
    size_t *maxVals=tempSpace;        
    maxVals[0]=1;
    for(size_t i=1;i<numDim+1;i++) {
        maxVals[i]=maxVals[i-1]*dims[i-1];
    }
    return index2NDimCPP(numDim, numIdx, idx, maxVals, newIndices);
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
