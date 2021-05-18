/**COMBINATORIALFUNS A header for C++ functions related to the generation
 *                  of combinatorial objects.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef COMBINATORIALFUNS
#define COMBINATORIALFUNS

//Defines the size_t type
#include <stddef.h>
//For floor.
#include <cmath>

//Prototypes of functions for computing the matrix permanent.
double permSquareCPP(const double *A, const size_t n, double *buffer);
double permCPP(const double *A, const size_t numRow, const size_t numCol, size_t *buffer);
double permCPPSkip(const double *A, const size_t numRowsTotal,const size_t * rows2Keep, const size_t *cols2Keep, const size_t numRowsKept, const size_t numColsKept,size_t *buffer);

/**GETNEXTCOMBOCPP Get the next combination in lexicographic order given
 *             the current combination. If the final combination in the
 *             sequence is passed, then the return value is true. This
 *             function generates the combination values starting from 0,
 *             not 1.  
 *
 *INPUTS: I A length r array holding the previous combination. The first
 *          combination is I=[0;1;2;...;r-1]. The modified combination is
 *          stored in here.
 *        n The number of items from which r items are chosen for
 *          combinations. The elements of I can range from 0 to n-1.
 *        r The number of items chosen from n and the length of I.
 *
 *OUTPUTS: The modified combination is stored in I. If the final
 *         combination is lexicographic order is passed (meaning that
 *         there are no more combinations, then the return value is true,
 *         otherwise it is false.
 *
 *The algorithm is from[1]. More details on the algorithm are in the
 *Matlab implementation.
 *
 *Because this is a template function, the entire function has to be
 *defined in the header file.
 *
 *REFERENCES:
 *[1] C. J. Mifsud, "Algorithm 154: Combination in lexicographical order," 
 *    Communications of the ACM, vol. 6, no. 3 pp. 103, Mar. 1963.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
template <class T>
bool getNextComboCPP(T *I,const size_t n,const size_t r) {
    if(I[r-1]<n-1) {
        I[r-1]++;
        return false;
    } else {
        size_t j, s;
        
        for(j=r-1;j>0;j--) {
           if(I[j-1]<n+j-r-1) {
               I[j-1]++;
        
               for(s=j;s<r;s++) {
                  I[s]=I[j-1]+s-(j-1); 
               }
               return false;
           }
        }
                   
        return true;
    }
}

/**GETNEXTCOMBOCPPFROMONE Get the next combination in lexicographic order
 *             given the current combination. If the final combination in
 *             the sequence is passed, then the return value is true. This
 *             function generates the combination values starting from 1,
 *             not 0.  
 *
 *INPUTS: I A length r array holding the previous combination. The first
 *          combination is I=[1;2;3;...;r]. The modified combination is
 *          stored in here.
 *        n The number of items from which r items are chosen for
 *          combinations. The elements of I can range from 1 to n.
 *        r The number of items chosen from n and the length of I.
 *
 *OUTPUTS: The modified combination is stored in I. If the final
 *         combination is lexicographic order is passed (meaning that
 *         there are no more combinations, then the return value is true,
 *         otherwise it is false.
 *
 *The algorithm is from [1]. More details on the algorithm are in the
 *Matlab implementation.
 *
 *Because this is a template function, the entire function has to be
 *defined in the header file.
 *
 *REFERENCES:
 *[1] C. J. Mifsud, "Algorithm 154: Combination in lexicographical order," 
 *    Communications of the ACM, vol. 6, no. 3 pp. 103, Mar. 1963.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
bool getNextComboCPPFromOne(T *I,const size_t n,const size_t r) {
    if(I[r-1]<n) {
        I[r-1]++;
        return false;
    } else {
        size_t j, s;
        
        for(j=r-1;j>0;j--) {
           if(I[j-1]<n+j-r) {
               I[j-1]++;
        
               for(s=j;s<r;s++) {
                  I[s]=I[j-1]+s-(j-1); 
               }
               return false;
           }
        }
                   
        return true;
    }
}

/*GETNEXTGRAYCODECPP A C++ function that returns the next gray code value.
*                 A gray code is a  binary code such that only one entry
*                 changes from 0 to 1 or back each step. This function can
*                 be used to get all subsets of an n-set.
*
*INPUTS: n The length of the code (const size_t).
*     code A pointer to a length-n array holding the previous code that
*          should be transformed to the next code value. This is modified
*          on return. It is a template, so different types can be passed.
*    nCard The number of ones in code, implicitly passed. This value is
*          modified on return (size_t&)
*        j A value, implicitly passed, that will be updated to hold the
*          index of the entry in code that was modified.
*
*RETURN VALUE: The boolean return value is 1 if the final combination was
*              returned and 0 otherwise.
*
*The algorithm is based on NEXSUB in Chapter 1 of [1]. Gray codes are
*also discussed in Chapter 7.2.1.1 of [2].
*
*The function checks against n so that if an invalid value of nCard is
*passed, it will not read/ write past the end of the array, even though it
*will return a code that is not the next in the sequence.
*
*Because this is a template function, the entire function has to be
*defined in the header file.
*
*REFERENCES:
*[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
*    and Calculators, 2nd ed. New York: Academic press, 1978.
*[2] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 2:
*    Generating all Tuples and Permutations, Upper Saddle River, NJ:
*    Addison-Wesley, 2009.
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
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

/**GETNEXTTUPLECPP This function produces the next tuple in a counting
 *     series. This can be used as a method of implementing nested loops.
 *     Each index of the tuple counts from 0 to the value in maxVals and
 *     resets. This is essentially a function for counting using
 *     different bases for each digit. This first tuple in the series is
 *     all zeros.
 *
 *INPUTS: numDim The number of dimensions of the tuple.
 *         tuple A length numDim array of elements representing the tuple.
 *               The entries should be integer numeric values. This is
 *               modified in place.
 *       maxVals A length numDim array whose elements correspond to the
 *               maximum value that each digit of the tuple can take (>=0).
 *               The elements of maxVals correspond to the base of each
 *               digit -1.
 * firstIsMostSig This is a boolean variable indicating whether tuple[0] is
 *               the most significant digit (or whether tuple[numDim-1] is
 *               the most significant). This affects the ordering of the
 *               tuples. The default if this parameter is omitted or an
 *               empty matrix is passed is true.
 *
 *OUTPUTS: The original tuple is modified with the new tuple. If the final
 *         tuple in the sequence is passed, meaning that there is no next
 *         tuple, then the return value is true, otherwise, it is false.
 *
 *Because this is a template function, the entire function has to be
 *defined in the header file.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
bool getNextTupleCPP(const ptrdiff_t numDim,T* tuple,const size_t *maxVals,const bool firstIsMostSig) {
    ptrdiff_t curLevel;
    bool isAscending;
    //tuple is the previous tuple.

    if(firstIsMostSig==false) {
        curLevel=0;
        isAscending=true;
        while(curLevel<numDim) {
            if(curLevel<0) {
                //If we have gotten a complete new tuple.
                return false;
            } else if(isAscending) {
                //Try incrementing the order at this level.
                tuple[curLevel]++;
                if(tuple[curLevel]<=maxVals[curLevel]) {
                    isAscending=false;
                    curLevel--;
                    continue;
                } else {
                    //If the value is invalid, then just keep ascending.
                    curLevel++;
                    continue;
                }
            } else {//We are descending in the sums here.
                tuple[curLevel]=0;
                curLevel=curLevel-1;
            }
        }
    } else {
        curLevel=numDim-1;
        isAscending=true;
        while(curLevel>=0) {
            if(curLevel>=numDim) {
                //If we have gotten a complete new tuple.
                return false;
            } else if(isAscending) {
                //Try incrementing the order at this level.
                tuple[curLevel]++;
                if(tuple[curLevel]<=maxVals[curLevel]) {
                    isAscending=false;
                    curLevel++;
                    continue;
                } else {
                    //If the value is invalid, then just keep ascending.
                    curLevel--;
                    continue;
                }
            } else {//We are descending in the sums here.
                tuple[curLevel]=0;
                curLevel++;
            }
        }
    }

    //If we get here, then we have gone past the final tuple.
    return true;  
}

template <class T>
size_t numTuples(const size_t N, const T *maxVals) {
/**NUMTUPLES Get the number of tuples possible when given a length-N array
 *           of maximum values to which each digit can go up. This is the
 *           number of tuples produces by the genAllTuplesCPP function. See
 *           the comments to that function for more information.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t tupleCount=maxVals[0]+1;
    size_t i;
    
    for(i=1;i<N;i++) {
        tupleCount*=(maxVals[i]+1);
    }
    
    return tupleCount;
}

/**GENALLTUPLESCPP Generate all tuples in a counting series. Each index of
 *           the tuple counts from 0 to the corresponding value in maxVals.
 *           This is just a method of counting uses different bases for
 *           each digit.
 *
 *INPUTS: theTuples The buffer in which the tuples will be places. The
 *               amount of memory needed for this can be found using the
 *               bufferSize4AllTuples function. length-N tuples are placed
 *               one after anouther in the buffer. It can be thought of as
 *               an array of tupleCount tuples of length-N stored
 *               consecutively. The ith tuple (i>=0) begins at
 *               theTuples+N*i.
 *             N The number of dimensions of the tuples.
 *       maxVals A length N array whose elements correspond to the
 *               maximum value that each digit of the tuple can take (>=0).
 *               The elements of maxVals correspond to the base of each
 *               digit -1.
 * firstIsMostSig This is a boolean variable indicating whether tuple[0] is
 *               the most significant digit (or whether tuple[numDim-1] is
 *               the most significant). This affects the ordering of the
 *               tuples.
 *
 *OUTPUTS: The tuples are put into the tuples input. The return value is
 *         the number of tuples added.
 *
 *This function just calls getNextTupleCPP in a loop.
 *
 *Because this is a template function, the entire function has to be
 *defined in the header file.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
size_t genAllTuplesCPP(T *theTuples, const size_t N, const size_t *maxVals, const bool firstIsMostSig) {
    T* curTuple=theTuples;
    
    //Fill in the first tuple.
    for(size_t i=0;i<N;i++) {
        curTuple[i]=0;
    }
    
    const size_t tupleCount=numTuples(N,maxVals);
    
    for(size_t k=1;k<tupleCount;k++) {
        //Copy the previous tuple into the memory for the current tuple.
        for(size_t i=0;i<N;i++) {
            curTuple[i+N]=curTuple[i];
        }
        curTuple+=N;
        //Overwrite the copied tuple with the new tuple.
        getNextTupleCPP(N,curTuple,maxVals,firstIsMostSig);
    }
    
    return tupleCount;
}

/**BUFFERSIZE4TUPLES Get the amount of memory in bytes needed to hold all
 *          possible tuples of length N where each digits has the upper
 *          bound specified in the corresponding entry of maxVals.
 *
 *INPUTS: N The number of dimensions of the tuples.
 *  maxVals A length N array whose elements correspond to the maximum value
 *          that each digit of the tuple can take (>=0). The elements of
 *          maxVals correspond to the base of each digit -1.
 *
 *OUTPUTS: The returned value si the number of bytes that should be
 *         allocated to hold all of the tuples.
 *
 *This is a template function, so the type must be specified.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
size_t bufferSize4AllTuples(const size_t N, const size_t *maxVals) {
    return sizeof(T)*N*numTuples(N,maxVals);
}

/**UNRANKTUPLECPP Obtain the tuple corresponding to its order in the
 *            sequence of tuples given a set of product terms derived from
 *            the maximum values that each tuple can take. If dims length
 *            numDim is the dimensionality of each dimension, then
 *            maxProdVals=[1;cumprod(maxVals+1)]. If one doesn't wish to
 *            precompute maxProdVals, then the function unrankTupleDimsCPP
 *            can be used with just the maximum values that each digit can
 *            take.
 *
 *INPUTS: numDim The dimensionality of the tuples.
 *        numTuples The number of indices provided that should be unranked
 *               to form tuples.
 *           idx A length numTuples array holding the indices to unrank.
 *   maxProdVals The product values related to the dimensions, as described
 *               above.
 *    newTuples A numDim*numTuples array into which the unranked tuples will
 *               be placed. These are placed one after the other.
 *
 *OUTPUTS: The return value indicates whether one is past the final tuple.
 *         If idx is beyond the final valid tuple, then nothing is put in
 *         newTuples and this function returns true. Otherwise, the
 *         function returns false and newTuples holds the new tuples.
 *
 *Non-integer types can be passed for the template values. They will be
 *converted to size_t values before being converted back upon return.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T1, class T2>
bool unrankTupleCPP(const size_t numDim, const size_t numTuples,const T1 * const idx, const T2 * const maxProdVals, T1 *newTuples) {
    const size_t maxVal=static_cast<size_t>(maxProdVals[numDim]);
    
    for(size_t i=0;i<numTuples;i++) {
        size_t curIndic=numDim;
        size_t curIdx=static_cast<size_t>(idx[i]);
        
        //The index is beyond the final one.
        if(curIdx>=maxVal) {
            return true;
        }
        
        do {
            curIndic--;
            size_t curMaxProdVal=static_cast<size_t>(maxProdVals[curIndic]);
            //integer division.
            const size_t wholeVal=curIdx/curMaxProdVal;
            
            newTuples[i*numDim+curIndic]=static_cast<T1>(wholeVal);
            curIdx-=wholeVal*curMaxProdVal;
        } while(curIndic>0);
    }
    return false;
}

/**UNRANKTUPLEDIMSCPP This is the same as the unrankTupleCPP function,
 *     except it is parameterized by the maximum values of the tuples
 *     indices rather than by the product vector that is used in
 *     unrankTupleCPP. The tempSpace buffer must have enough space to hold
 *     at least numDim+1 values of type T.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template <class T>
bool unrankTupleDimsCPP(const size_t numDim, const size_t numTuples,const T *idx, const T *maxVals, T *newTuples, size_t *tempSpace) {
    //This memory buffer is numDim+1 in size.
    size_t *maxProdVals=tempSpace;        
    maxProdVals[0]=1;
    for(size_t i=1;i<numDim+1;i++) {
        maxProdVals[i]=maxProdVals[i-1]*(static_cast<size_t>(maxVals[i-1])+1);
    }
    return unrankTupleCPP(numDim, numTuples, idx, maxProdVals, newTuples);
}

/**UNRANKFLOATTUPLECPP This is the same as the unrankTupleCPP function,
 *     except the template type T is assumed to be a floating point type.
 *     Here, no integer conversions are performed.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
bool unrankFloatTupleCPP(const size_t numDim, const size_t numTuples,const T *const idx, const T * const maxProdVals, T *newTuples) {
    const double maxVal=maxProdVals[numDim];
    
    for(size_t i=0;i<numTuples;i++) {
        size_t curIndic=numDim;
        T curIdx=idx[i];
        
        //The index is beyond the final one.
        if(curIdx>=maxVal) {
            return true;
        }
        
        do {
            curIndic--;
            T curMaxProdVal=maxProdVals[curIndic];
            //integer division.
            const double wholeVal=floor(curIdx/curMaxProdVal);
            
            newTuples[i*numDim+curIndic]=wholeVal;
            curIdx-=wholeVal*curMaxProdVal;
        } while(curIndic>0);
    }
    return false;
}

/**UNRANKFLOATTUPLEDIMSCPP This is the same as the unrankTupleDimsCPP
 *     function, except the template type T is assumed to be a floating
 *     point type. Here, no integer conversions are performed. The
 *     tempSpace buffer must have space for numDim+1 elements of type T.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
bool unrankFloatTupleDimsCPP(const size_t numDim, const size_t numTuples,const T *idx, const T * const maxVals, T *newTuples, T *tempSpace) {
    //This memory buffer is numDim+1 in size.
    T *maxProdVals=tempSpace;        
    maxProdVals[0]=1;
    for(size_t i=1;i<numDim+1;i++) {
        maxProdVals[i]=maxProdVals[i-1]*(maxVals[i-1]+1);
    }
    return unrankFloatTupleCPP(numDim, numTuples, idx, maxProdVals, newTuples);
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
