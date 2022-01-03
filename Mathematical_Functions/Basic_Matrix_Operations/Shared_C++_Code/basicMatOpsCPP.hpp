/**BASICMATOPSCPP C++ functions of basic operations with matrices. See the
 *             comments in the Matlab implementations for descriptions of
 *             the functionality.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef BASICMATOPSCPP
#define BASICMATOPSCPP

//Defines the size_t and ptrdiff_t types
#include <cstddef>
//For sort.
#include <algorithm>

/**COMPRESSNONSELDIMSCPP Given a length S set of values in nDims and
 *              selDims holds indices selecting numSel of those values
 *              (0<=numSel<=S), replace nDims with a set of values where
 *              all those between ones that are selected are compressed
 *              into a single value (via multiplication). The value
 *              numNewDims holds the length of the modified nDims. Also,
 *              selDims is updated to reflect the new selected dimensions.
 *              selDims(i) on the input if selDims on the input was sorted.
 *              This algorithm relates to reshaping a matrix collapsing the
 *              non-selected dimensions.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
void compressNonSelDimsCPP(const size_t S, size_t *nDims, T *selDims,const size_t numSel,size_t &numNewDims) {          
    numNewDims=0;
    
    //Sort in ascending order.
    std::sort(selDims,selDims+numSel);

    size_t curDim=0;
    for(size_t curSelIdx=0;curSelIdx<numSel;curSelIdx++) {
        size_t curSelDim=static_cast<size_t>(selDims[curSelIdx]);
        if(curDim<curSelDim) {
            //If there is a gap between the current dimension and the
            //selected dimension, everything from curDim to curSelDim-1
            //gets compressed together.
            size_t prodVal=nDims[curDim];
            for(size_t k=curDim+1;k<curSelDim;k++) {
                prodVal*=nDims[k];
            }

            nDims[numNewDims]=prodVal;
            numNewDims++;
        }
            
        //When here, curDim=curSelDim. Thus, we add the selected dimension
        //and we also record the index in the newSelDims vector.
        nDims[numNewDims]=nDims[curSelDim];
        //Update the index of the currently selected dimension.
        selDims[curSelIdx]=static_cast<T>(numNewDims);
        numNewDims++;
        curDim=curSelDim+1;
    }

    //If there are dimensions after the last selected dimension, compress
    //them.
    if(curDim<S) {
        size_t prodVal=nDims[curDim];
        for(size_t k=curDim+1;k<S;k++) {
            prodVal*=nDims[k];
        }
        nDims[numNewDims]=prodVal;
        numNewDims++;
    }
}

/**COMPRESSNONSELDIMSBOOL Given a length S set of values in nDims and
 *              selDims holds a length S set of boolean values selecting
 *              values, replace nDims with a set of values where all those
 *              between ones that are selected are compressed into a single
 *              value (via multiplication). The first numNewDims elements
 *              of selDims are updated to mark which of the collapsed
 *              dimensions correspond to the original selected dimensions.
 *              numNewDims holds the length of the modified nDims. This
 *              algorithm relates to reshaping a matrix collapsing the non-
 *              selected dimensions.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
template <class T>
void compressNonSelDimsBool(const size_t S, size_t *nDims, T *selDims,size_t &numNewDims) { 
//This function overwrites nDims and selDims with the new values.
//newSelDims is length S.

    numNewDims=0;
    size_t curDim=0;
    while(curDim<S) {
        if(selDims[curDim]) {//If it is a fixed dimension.
            nDims[numNewDims]=nDims[curDim];
            //Overwrite the old selDims.
            selDims[numNewDims]=true;
            numNewDims++;
            curDim++;
            continue;
        } else {
            //If is not a fixed dimension, we have to scan forward until we
            //either find the next fixed dimension or the end.
            size_t prodVal=nDims[curDim];
            curDim++;
            while(curDim<S&&selDims[curDim]==false) {
                prodVal*=nDims[curDim];
                curDim++;
            }
            
            nDims[numNewDims]=prodVal;
            //This is not a selected dimension.
            selDims[numNewDims]=false;
            numNewDims++;
            continue;
        }
    }
}


template<class T>
void permuteMatrixColsCPP(const size_t numRows,const size_t numCols, T*V,const size_t * const idx) {
    for(size_t k=0;k<numCols;k++) {
        size_t kn=idx[k];
        while(kn<k) {
            kn=idx[kn];
        }
        if(kn!=k) {
            for(size_t i=0;i<numRows;i++) {
                T *vkn=V+i+numRows*kn;
                T *vk=V+i+numRows*k;
                T wr=*vkn;
                *vkn=*vk;
                *vk=wr;
            }
        }
    }    
}

template<class T>
void permuteVectorCPP(const size_t numEls, T*v,const size_t * const idx) {
    for(size_t k=0;k<numEls;k++) {
        size_t kn=idx[k];
        while(kn<k) {
            kn=idx[kn];
        }
        if(kn!=k) {
            T wr=v[kn];
            v[kn]=v[k];
            v[k]=wr;
        }
    }    
}

template<class T>
void permuteMatrixRows(const size_t numRows, const size_t numCols, T *v, const size_t * const idx) {
    for(size_t k=0;k<numCols;k++) {
        permuteVectorCPP(numRows,v,idx);
        v+=numRows;
    }
}


//The functions below implement things that are trivial in Matlab but
//complicated otherwise. They do not have an equivalent files implementing
//them in Matlab (one could do it with one line).

template<class T1, class T2>
bool anyArralElsEqual(T1 *v1, T2 *v2, const size_t numEls) {
/*ANYARRAYELSEQUAL Given two arrays of length numEls, this function returns
 *                 true if any corresponding element in both arrays is
 *                 equal.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    for(size_t i=0;i<numEls;i++) {
        if(v1[i]==v2[i]) {
            return true;
        }
    }
    return false;
}

template<class T1,class T2>
bool allArrayElsEqual(T1 *v1, T2 *v2, const size_t numEls) {
/*ALLARRAYELSEQUAL Given two arrays of length numEls, this function returns
 *                 true if all corresponding element in both arrays are
 *                 equal.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    for(size_t i=0;i<numEls;i++) {
        if(v1[i]!=v2[i]) {
            return false;
        }
    }
    return true;
}

template<class T>
bool anyNonzeroArralElsEqual(T *v1, T *v2, const size_t numEls) {
/*ANYNONZEROARRAYELSEQUAL Given two arrays of length numEls, this function
 *                 returns true if any nonzero corresponding element in
 *                 both arrays is equal.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    for(size_t i=0;i<numEls;i++) {
        if(v1[i]==0||v2[i]==0) {
            continue;
        }
        
        if(v1[i]==v2[i]) {
            return true;
        }
    }
    return false;
}

template <class T>
void copyMat3DOmitDims(T * const dest, const size_t * const sourceMatDims, const T * const source,const size_t *skipIdx) {
/**COPYMAT3DOMITDIMS This function copies a 3D matrix into a buffer,
 *      skipping a single index in each dimension while doing so. If one
 *      does not wish to skip a particular dimension, one can set the skip
 *      index in that dimension to be >=the number of elements in that
 *      dimension.
 * 
 * Using Matlab notation, but with indexation from 0
 *  instead of 1, this function performs the assignment
 *  dest(:)=vec(source([0:(i-1),(i+1):(sourceMatDims[0]-1)],[0:(j-1),(j+1):(sourceMatDims[1]-1)],[0:(k-1),(k+1):(sourceMatDims[2]-1)]))
 *  where i=skipIdx[0],j=skipIdx[1], and k=skipIdx[2];
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t source01Offset=sourceMatDims[0]*sourceMatDims[1];
    
    size_t destIdx0=0;
    for(size_t idx2=0;idx2<sourceMatDims[2];idx2++) {
        const size_t offset2=(idx2>=skipIdx[2]);
        const size_t sourceIdx2=(idx2+offset2)*source01Offset;

        if(idx2+offset2==sourceMatDims[2]) {
            break;
        }
        
        for(size_t idx1=0;idx1<sourceMatDims[1];idx1++) {
            const size_t offset1=(idx1>=skipIdx[1]);
            const size_t sourceIdx1=sourceIdx2+(idx1+offset1)*sourceMatDims[0];

            if(idx1+offset1==sourceMatDims[1]) {
                break;
            }
            
            for(size_t idx0=0;idx0<sourceMatDims[0];idx0++) {
                const size_t offset0=(idx0>=skipIdx[0]);
                const size_t sourceIdx0=sourceIdx1+idx0+offset0;
                
                if(idx0+offset0==sourceMatDims[0]) {
                    break;
                }
                
                dest[destIdx0]=source[sourceIdx0];
                destIdx0++;
            }
        }
    }
}

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
