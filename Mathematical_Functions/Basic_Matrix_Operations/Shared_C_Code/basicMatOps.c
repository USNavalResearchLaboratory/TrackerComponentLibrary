/**BASICMATOPS Implementations of functions for performing basic matrix
 *      operations that are simple to perform in matlab, but that can be
 *      tedious when programming in C. See the comments to the individual
 *      implementations for descriptions of the functions. This function
 *      only implements the functions in basicMatOps.h that are provided
 *      without a Matlab interface.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "basicMatOps.h"

//For memset.
#include <string.h>

/*We need this for INFINITY to be defined, but if a compiler does not
 *support C99, then it must be explicitly defined.*/
#include <math.h>

//For uint64_t and other types
#include <stdint.h>

/*If a compiler does not support INFINITY in C99, then it must be
 * explicitly defined.*/
#ifndef INFINITY
static const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)(&infVal))
#endif

//If using *NIX or Mac OS X, the fmin be defined without anything extra in
//math.h
#if __STDC_VERSION__>=199901L
#include<math.h>
#define fMin(a,b) fmin(a,b)
#else
    static double fMin(double a,double b) {
    //Defined in the C99 standard, but not supported by Microsoft.
        if(a>b) {
            return b;
        }
        else {
            return a;
        }
    }
#endif

void permuteRowsPtrDiffT(const size_t n1,const size_t n2,ptrdiff_t *  XPerm,const ptrdiff_t *  XOrig,const size_t *rowOrder) {
/**PERMUTEROWSPTRDIFFT Given an n1Xn2 matrix of ptrdiff_t values stored by
 *              column in XOrig, permute the rows of the matrix according
 *              to the ordering in rowOrder (a permutation of 0 to (n-1))
 *              and store the result in XPerm.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t i1, i2;
    
    for(i2=0;i2<n2;i2++) {
        const size_t n1i2=n1*i2;
        for(i1=0;i1<n1;i1++) {
            XPerm[rowOrder[i1]+n1i2]=XOrig[i1+n1i2];
        }
    }
}

void identMatD(const size_t numDims, double *  I) {
/**IDENTMATD Set the doubles in I to be a numDimsXnumDims identity matrix,
 *           stored column-wise. I is a numDims*numDims array of doubles.
 *  
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t i;
    
    //Set I to zero.
    memset(I,0,numDims*numDims*sizeof(double));
    
    //Set the diagonal elements to 1.
    for(i=0;i<numDims;i++) {
        I[i+numDims*i]=1;
    }
}

double vecMinD(const double*vec, const size_t numEls) {
/**VECMIND Given an array containing numEls elements, return the minimum
 *         value. If numEls=0, then 0 is returned.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    double minVal;
    size_t i;
    
    minVal=0;
    for(i=0;i<numEls;i++) {
        minVal=fMin(minVal,vec[i]);
    }
    return minVal;
}


double sumVectorD(const double *vector, const size_t numEls) {
/**SUMVECTORD  Given a vector of doubles of length numEls, sum all of the
 *             values in the vector. If numEls=0, then 0 is returned.
 *  
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    double sumVal;
    size_t i;
    
    sumVal=0;
    for(i=0;i<numEls;i++) {
        sumVal+=vector[i];
    }
    return sumVal;
}

size_t sumVectorSizeT(const size_t *vector, const size_t numEls) {
/**SUMVECTORSIZET Given a vector of size_t elements of length numEls, sum
 *                all of the values in the vector. If numEls=0, then 0 is
 *                returned.
 *  
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t sumVal;
    size_t i;
    
    sumVal=0;
    for(i=0;i<numEls;i++) {
        sumVal+=vector[i];
    }
    return sumVal;
}

size_t prodVectorSizeT(const size_t *vector, const size_t numEls) {
/**PRODVECTORD Given a vector of length numEls, multiply all of the values
 *             in the vector. If numEls=0, then 1 is returned.
 *  
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t prodVal;
    size_t i;
    
    prodVal=1;
    for(i=0;i<numEls;i++) {
        prodVal*=vector[i];
    }
    return prodVal;
}

size_t accessHyperMatSizeT(const size_t *prodTerms, const size_t *C,const size_t *idx,const size_t S) {
/**ACCESSHYPERMATPTRDIFFT Given an S-dimensional hypermatrix, get the
 *       element in the position specified by the length-S index vector idx
 *       This is like the indexation C(idx(1),idx(2),..,idx(S)) in Matlab,
 *       except indexation of everything starts at 0. This variant of the
 *       function is for C being a matrix of size_t values.
 *
 *INPUTS: prodTerms If we say that nVals[i] is the length of the ith
 *                  dimension of C (starting at 0), then
 *                  prodTerms[0]=nVals[0], prodTerms[1]=nVals[0]*nVals[1],
 *                  prodTerms[2]=nVals[0]*nVals[1]*nVals[2], etc. When
 *                  repreatedly accessing a hypermatrix, it is more
 *                  efficient to precompute all of the product terms, which
 *                  is why this matrix is passed instead of nVals.
 *                C A pointer to the hypermatrix.
 *              idx A pointer to the indices (starting at 0).
 *                S The number of dimensions of C, the length of the vector
 *                  idx.
 * 
 *OUTPUTS: The return value is the element in C given by the specified
 *         index vector idx.
 *
 *OUTPUTS: The outputs are placed in d and dIdx.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */  
    size_t i;
    size_t cumOffset;
    
    cumOffset=idx[0];
    for (i=1;i<S;i++) {
        cumOffset+=prodTerms[i-1]*idx[i];
    }
    
    return C[cumOffset];
}

void addVec2MatEndMin(double *d,size_t *dIdx,const double *C,const double *u,const size_t numElsInD, const size_t numElsInU) {
/**ADDVEC2MATENDMIN Given a numElsInDXnumElsInU matrix C, and a length
 *            numElsInU vector u, this function performs an operation
 *            equivalent to the Matlab line
 *            [d,dIdx]=min(bsxfun(@ plus,C,reshape(u,[2,numElsInU])),[],2);
 *            Except the indexation here is from 0, not 1. This function
 *            can be used with hypermatrices where one wishes to apply the
 *            same type of operation to the final dimension in that one
 *            simply collapses all elements in the preceding dimensions
 *            into the rows.           
 *
 *INPUTS: d The buffer to store the result of the minimization that will
 *          hold numElsInD elements.
 *     dIdx A buffer containing the same number of elements as d that will
 *          store the minimum index of C(curEl,:)+u(:) for each element in
 *          d. This index ranges from 0 to numElsInU-1.
 *        C The input matrix. This is stored column-wise and has dimensions
 *          numElsInDXnumElsInU
 *        u The length numElsInU vector to add to the columns.
 * numElsInD The number of rows in C.
 * numElsInU The number of columns in C.
 *
 *OUTPUTS: The outputs are placed in d and dIdx.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    size_t i;

    for(i=0;i<numElsInD;i++) {
        double minVal=(double)INFINITY;
        size_t minIdx=0;
        size_t curU;
        
        for(curU=0;curU<numElsInU;curU++) {
            const size_t curIdx=i+numElsInD*curU;
            double curVal=C[curIdx]+u[curU];

            if(curVal<minVal) {
                minVal=curVal;
                minIdx=curU;
            }
        }
        d[i]=minVal;
        dIdx[i]=minIdx;
    }
}

void fixedIdxMatSub(const size_t S, const size_t *nVals,const double *M, double *  C,const size_t dim, const size_t idx) {
/**FIXEDIDXMATSUB  This function implements the matrix operation equivalent
 *    to the Matlab line C(:,:,...,:,idx,:)=C(:,:,...,:,idx,:)-M where one
 *    index is fixed. idx is the value of the fixed index and it is taken
 *    in dimension dim (dim indexed from 0). For an S-dimensional matrix,
 *    if dim>=S (S counting from 1), then it is assumed that idx=0 and one
 *    just evaluates C=C-M
 *
 *INPUTS: S The number of dimensions of the matrix (Starting at 1).
 *    nVals A length S array of the sizes of each dimension.
 *        M The matrix M as given above, having a length being the product
 *          of all elements in nVals except nVals(dim).
 *        C The multidimensional matrix C given above. This matrix is
 *          modified to hold the return value.
 *      dim A value indicating which dimension has the fixed index. This
 *          starts from 0.
 *      idx The value of the fixed index. This must be <nVals[dim].
 *
 *OUTPUTS: The output is placed in C.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    if(dim>=S) {
        size_t i;
        size_t totalNumEls=prodVectorSizeT(nVals,S);
        
        for(i=0;i<totalNumEls;i++) {
            C[i]-=M[i];
        }
        return;
    }
    
    //n1=prod(nVals(1:(dim-1))); in Matlab
    const size_t n1=prodVectorSizeT(nVals,dim);
    //totalNumElsM=n1*prod(nVals((dim+1):S)); in Matlab
    const size_t totalNumElsM=n1*prodVectorSizeT(nVals+dim+1,S-dim-1);
    const size_t incrBigStep=n1*nVals[dim];
    size_t curEl,CStartIdx,curMCol;
    
    CStartIdx=idx*n1;
    curMCol=0;
    for(curEl=0;curEl<totalNumElsM;curEl++) {
        C[CStartIdx]-=M[curEl];
        
        if((curEl+1)%n1==0) {
            curMCol++;

            CStartIdx=incrBigStep*curMCol+idx*n1;
        } else {
            CStartIdx++;
        }
    }
}

size_t sub2Idx2D(const size_t idx1,const size_t idx2,const size_t dim1) {
/**SUB2IDX2D Convert a set of two subscripts to an index for a 2D matrix.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    return idx1+idx2*dim1;
}

size_t sub2Idx3D(const size_t idx1,const size_t idx2,const size_t idx3,const size_t dim1,const size_t dim2) {
/**SUB2IDX3D Convert a set of three subscripts to an index for a 3D matrix.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    return idx1+idx2*dim1+idx3*dim1*dim2;
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
