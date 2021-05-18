/**PERMUTEMATRIXC C-only implementationa of the permuteMatrix function and
 *             a helper function. See the comments to the original
 *             functions below  as well as to permuteMatrix in Matlab for
 *             details.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "basicMatOps.h"

//For memcpy
#include <string.h>

size_t permuteMatrixCBufferSize(const size_t S) {
/**PERMUTEMATRIXCBUFFERSIZE This function returns the minimum size (in
 *          bytes) that the buffer tempBuffer for the function
 *          permuteMatrixC must be. S is the number of dimensions of the
 *         (hyper) matrix C. For S<4, this is 0 - not buffer is needed.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */

    if(S>3) {
        return 4*S*sizeof(size_t);
    } else {
        return 0;   
    }
}

void permuteMatrixC(const size_t S,const size_t *nVals, size_t *  nValsNew, void* CPermMat, const void *C,const size_t CElSize, void *tempBuffer, const size_t *order) {
/**PERMUTEMATRIXC A C implementation fo a function to perform a generalized
 *     matrix tranpose. Given an n1Xn2X...nS matrix, rearrange the
 *     dimensions of C to be in the order specified by order. All elements
 *     of order must be unique, positive, integer values from 0 to S-1.
 *     This functions in the same manner as Matlab's permute function. The
 *     permutation performed by this function is not done in-place. This
 *     type of permutation is a type of tensor matrix transpose.
 *
 *INPUTS: S The number of dimensions of the matrix. S>0.
 *    nVals A length S array holding the dimensions of the matrix to be
 *          permuted. Each dimension must be >=1.
 * nValsNew An array that will hold nDims in the permuted ordering.
 * CPermMat A pointer to the space where the permuted matrix should be
 *          stored. This cannot be the same as the input C.
 *        C A pointer to the original matrix that is to be permuted.
 *          Elements are stored by columns, then rows, etc, as in Matlab.
 *  CElSize The number of bytes in each element of C and CPermMat.
 * tempBuffer A buffer of at least permuteMatrixCBufferSize(S) bytes to
 *          hold temporary results.
 *    order A length S array holding the values from 0 to S-1 with no
 *          repeats.
 *
 *OUTPUTS: The results are placed in CPermMat with nValsNew holding the
 *         permutation of nVals.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    size_t i, iNew, totalNumEls;
    size_t newLinIdx, linIdx;
    char *  CPerm=(char*)CPermMat;

    if(S==0) {
        return;
    } else if(S==1) {
        //CPerm[0]=C[0] but for elements of varying sizes.
        memcpy(CPerm,C,CElSize);
        return;
    } else if(S==2) {
        permute2DimsC(nVals,nValsNew,CPermMat,C,CElSize,order);
        return;
    } else if(S==3) {
        permute3DimsC(nVals,nValsNew,CPermMat,C,CElSize,order);
        return;
    }
    
    totalNumEls=nVals[0];
    for(i=1;i<S;i++) {
       totalNumEls*=nVals[i];
    }
    
    if(totalNumEls==0) {
        return;
    }
    
    //Partition the buffer
    size_t *  cumProdNew=(size_t*)tempBuffer;
    tempBuffer=(void*)((size_t*)tempBuffer+S);
    
    size_t *  invPerm=(size_t*)tempBuffer;
    tempBuffer=(void*)((size_t*)tempBuffer+S);
    
    size_t *  idx=(size_t*)tempBuffer;
    tempBuffer=(void*)((size_t*)tempBuffer+S);
    
    size_t *  newIdx=(size_t*)tempBuffer;

    //Initialize the indices.
    memset(idx,0,sizeof(size_t)*S);
    memset(newIdx,0,sizeof(size_t)*S);
 
    iNew=order[0];
    nValsNew[0]=nVals[iNew];
    invPerm[iNew]=0;
    cumProdNew[0]=1;

    for(i=1;i<S;i++) {
        iNew=order[i];
        
        nValsNew[i]=nVals[iNew];
        invPerm[iNew]=i;
        
        cumProdNew[i]=cumProdNew[i-1]*nValsNew[i-1];
    }
    
    newLinIdx=0;
    linIdx=0;
    while(1) {
        size_t curLevel;
        
        //Assign the current tuple. CPerm[newLinIdx]=C[linIdx];
        memcpy(CPerm+newLinIdx*CElSize,(const char*)(C)+linIdx*CElSize,CElSize);

        //Move on to the next tuple.
        linIdx++;
        
        if(linIdx>=totalNumEls) {
            return;
        }
        //The code below is similar to that in getNextTuple, because we are
        //updating the tuples for newIdx.
        curLevel=0;
        while(1) {
            //Try incrementing the order at this level.
            idx[curLevel]++;
            if(idx[curLevel]<nVals[curLevel]) {
                iNew=invPerm[curLevel];
                newIdx[iNew]++;
                newLinIdx+=cumProdNew[iNew];

                //In Matlab, this would be idx(1:(curLevel-1))=0.
                memset(idx,0,curLevel*sizeof(size_t));
                break;
            } else {
                //If the value is invalid, then just keep ascending and
                //adjust newLinIdx.
                iNew=invPerm[curLevel];

                newLinIdx-=newIdx[iNew]*cumProdNew[iNew];
                newIdx[iNew]=0;

                curLevel++;
                continue;
            }
        }
    }
}


void permute2DimsC(const size_t *nDims,size_t *  nDimsNew, void * CPermMat,const void * COrigMat,const size_t CElSize,const size_t *dimsOrder) {
/**PERMUTE2DIMS Given a 2D matrix, COrig, copy it into CPerm while
 *          rearranging the order of the dimensions according to dimsOrder.
 *          The matrices are stored by column (the same ordering as in
 *          Matlab and Fortran).
 *
 *INPUTS: nDims A length-2 vector with the size of each dimension of the
 *              2D matrix COrig.
 *     nDimsNew A length 2 vector into which the permuted nDims is placed.
 *        CPerm The space to hold the permuted matrix. This must be at
 *              least nDims[0]*nDims[1]*sizeof(double) bytes in
 *              size. This cannot be the same as COrig.
 *        COrig The original 2D matrix whose dimensions should be permuted.
 *    dimsOrder A length-2 vector whose ordering determines how the indices
 *              of COrig will be permited. This is a permutation of the
 *              numbers 0 and 1.
 * 
 *OUTPUTS: The result is put in CPerm and nDimsNew.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    if(dimsOrder[0]==0&&dimsOrder[1]) {
        nDimsNew[0]=nDims[0];
        nDimsNew[1]=nDims[1];
        memcpy(CPermMat,COrigMat,CElSize*nDims[0]*nDims[1]);
    } else {
        const size_t n1=nDims[0];
        const size_t n2=nDims[1];
        size_t i1,i2;
        char *CPerm=(char*)CPermMat;
        const char *C=(const char*)COrigMat;
        
        nDimsNew[0]=nDims[1];
        nDimsNew[1]=nDims[0];
        
        for(i2=0;i2<n2;i2++) {
            const size_t n2Offset=n1*i2;
            for(i1=0;i1<n1;i1++) {
                memcpy(CPerm+(i2+i1*n2)*CElSize,C+(i1+n2Offset)*CElSize,CElSize);
            }
        }
    }
}

void permute3DimsC(const size_t *nDims,size_t *  nDimsNew,void * CPermMat,const void* COrigMat, const size_t CElSize,const size_t *dimsOrder) {
/**PERMUTE3DIMSC Given a 3D matrix, COrig, copy it into CPerm while
 *          rearranging the order of the dimensions according to dimsOrder.
 *          The matrices are stored by column (the same ordering as in
 *          Matlab and Fortran).
 *
 *INPUTS: nDims A length-3 vector with the size of each dimension of the
 *              3D matrix COrig.
 *     nDimsNew A length 3 vector into which the permuted nDims is placed.
 *        CPerm The space to hold the permuted matrix. This must be at
 *              least nDims[0]*nDims[1]*nDims[2]*CElSize bytes in
 *              size. This cannot be the same as COrigMat.
 *        COrig The original 3D matrix whose dimensions should be permuted.
 *      CElSize The number of bytes in each element in COrig
 *    dimsOrder A length-3 vector whose ordering determines how the indices
 *              of COrig will be permited. This is a permutation of the
 *              numbers 0, 1 and 2.
 * 
 *OUTPUTS: The result is put in CPerm and nDimsNew.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t n1=nDims[0];
    const size_t n2=nDims[1];
    const size_t n1n2=n1*n2;
    const size_t n3=nDims[2];
    size_t i1,i2,i3,newIdx[3], invPerm[3];
    char *  CPerm=(char*)CPermMat;
    const char *COrig=(const char*)COrigMat;
    
    nDimsNew[0]=nDims[dimsOrder[0]];
    nDimsNew[1]=nDims[dimsOrder[1]];
    
    const size_t n1n2New=nDimsNew[0]*nDimsNew[1];
    
    nDimsNew[2]=nDims[dimsOrder[2]];
    
    invPerm[dimsOrder[0]]=0;
    invPerm[dimsOrder[1]]=1;
    invPerm[dimsOrder[2]]=2;
    
    for(i3=0;i3<n3;i3++) {
        const size_t n1n2i3=n1n2*i3;

        newIdx[invPerm[2]]=i3;
        for(i2=0;i2<n2;i2++) {
            const size_t n1i2=n1*i2;

            newIdx[invPerm[1]]=i2;
            for(i1=0;i1<n1;i1++) {
                newIdx[invPerm[0]]=i1;

                //For the given datatype, this is assigning:
                //CPerm[newIdx[0]+nDimsNew[0]*newIdx[1]+n1n2New*newIdx[2]]=COrig[i1+n1i2+n1n2i3];
                memcpy(CPerm+(newIdx[0]+nDimsNew[0]*newIdx[1]+n1n2New*newIdx[2])*CElSize,COrig+(i1+n1i2+n1n2i3)*CElSize,CElSize);
            }
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
