/**MISCFUNS A variety of miscellansous functions.
 *
 *October 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MISC_FUNS
#define MISC_FUND

#include<cstddef>

//See the comments in index2NDimCPP.cpp for a description of the inputs.
bool index2NDimCPP(const size_t numDim, const size_t numIdx,const size_t * const idx, const size_t * const maxVals, size_t *newIndices);
bool index2NDimByDimCPP(const size_t numDim, const size_t numIdx,const size_t *idx, const size_t *dims, size_t *newIndices, size_t *tempSpace);

#endif 

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
