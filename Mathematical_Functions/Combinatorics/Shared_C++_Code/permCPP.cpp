/**PERMCPP C++ functions for computing the matrix permanent allowing for
*         rows and columns to be skipped is desired (operate on a
*         submatrix). The permanent is equivalent to calculating the
*         determininant in the standard summation manner taught in school,
*         except all of the minus signs are replaced with plus signs.
*
*The functions are described in more detail in the Matlab file perm.m and
*the interface to tese functions, perm.cpp.
*
*The functions in the file are:
*double permSquareCPP(const double *A, const size_t n,double *buffer)
*-To compute the permanent of the matrix A (stored in Matlab's format by
*row) which is nXn and where buffer is an array of at least size 2*n-1,
*which will be used for scratch space suring the computation.
*
*double permCPP(const double *A, const size_t numRow, const size_t numCol,
*               size_t *buffer)
*-To compute the permanent of a numRowXnumCol matrix. The buffer must be at
*least numCol in length.
*
*double permCPPSkip(const double *A, const size_t numRowsTotal,
*                   const size_t * rows2Keep, const size_t *cols2Keep,
*                   const size_t numRowsKept, const size_t numColsKept,
*                   size_t *buffer)
*-To compute the permanent of a numRowsKeptXnumColsKept submatrix of A. The
*array rows2Keep stores the indices of the rows of A that are kept (indexed
*from 0) and the array cols2Keep stores the indices of the columns to keep
*when computing the matrix permanent. The value numRowsTotal is the number
*of rows in the entire A matrix. The buffer must be at least numColsKept in
*length.
*
*The functions
*double SigmaS(const double *A,size_t *curComb,const size_t r,const size_t numRow,const size_t numCol);
*and 
*double SigmaSSkip(const double *A,size_t *curComb,const size_t *rows2Keep, const size_t *cols2Keep,const size_t r,const size_t numRowsTotal,const size_t numRowsKept,const size_t numColsKept);
*are subroutines used in the computations with Ryser's algorithm.
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "combinatorialFuns.hpp"

//For the fill function
#include <algorithm>

//Prototypes of helper functions used in this file (defined at the bottom of the file.
double SigmaSSkip(const double *A,size_t *curComb,const size_t *rows2Keep, const size_t *cols2Keep,const size_t r,const size_t numRowsTotal,const size_t numRowsKept,const size_t numColsKept);
double SigmaS(const double *A,size_t *curComb,const size_t r,const size_t numRow,const size_t numCol);

//Compute the matrix permanent using the efficient algorithms of Nijenhuis
//and Wilf for square matrices.
double permSquareCPP(const double *A, const size_t n, double * buffer) {
    //The buffer must be at least of size 2*n-1
    double *x,*code;
    double p,sgn;
    size_t i,j,nCard;
    bool isLast;
    
    if(n==1) {
        return *A;
    }
    
    //Allocate temporary space and split it between x and code.
    x = buffer;
    code=buffer+n;
    
    p=1;
    for(i=0;i<n;i++) {
        double sumVal=0;
        
        for(j=0;j<n;j++){
            sumVal+=A[i+n*j];
        }
        x[i]=A[i+n*(n-1)]-sumVal/2;
        p*=x[i];
    }
    
    //Set all the entries in the code array to zero. That is the first gray
    //code in the sequence.
    std::fill_n(code, n-1, 0);
   
    nCard=0;
    sgn=1;
    do {
        double z,prodVal;
        
        sgn=-sgn;
        isLast=getNextGrayCodeCPP(n-1,code,nCard,j);
        
        z=2*code[j]-1;
        prodVal=1;
        for(i=0;i<n;i++) {
            x[i]=x[i]+z*A[i+n*j];
            prodVal*=x[i];
        }
        
        p+=sgn*prodVal;
    } while(isLast==false);
    
    return 2.0*(2.0*static_cast<double>(n%2)-1.0)*p;
}

//Compute the matrix permanent using Ryser's method.
double permCPP(const double *A, const size_t numRow, const size_t numCol, size_t *buffer) {
    //The buffer must be at least size numCol.
    size_t x, *curComb;
    double retVal,binomTerm; 
    //If is assumed that numCol>=numRow

    //Empty matrices have a permanent of one by definition. 
    if(numRow==0||numCol==0) {
        return 1;
    }
    
    //Allocate space for the buffer that holds the combination of columns
    //chosen.
    curComb=buffer;
    
    binomTerm=1;
    retVal=0;
    for(x=0;x<numRow;x++) {
        retVal+=SigmaS(A,curComb,numCol-numRow+x,numRow,numCol)*binomTerm;        
        binomTerm=binomTerm*(1-numRow+numCol+x)/(1+x)*(-1);
    }
    
    return retVal;
}


//Compute the matrix permanent using Ryser's method, allowing for rows and columns to be skipped
double permCPPSkip(const double *A, const size_t numRowsTotal,const size_t * rows2Keep, const size_t *cols2Keep, const size_t numRowsKept, const size_t numColsKept, size_t *buffer) {
//The buffer must be at least numColsKept in length.
    size_t x, *curComb;
    double retVal,binomTerm; 
    //If is assumed that numCol>=numRow

    //Empty matrices have a permanent of one by definition. 
    if(numRowsKept==0||numColsKept==0) {
        return 1;
    }
    
    //Allocate space for the buffer that holds the combination of columns
    //chosen.
    curComb=buffer;
    
    binomTerm=1;
    retVal=0;
    for(x=0;x<numRowsKept;x++) {
        retVal+=SigmaSSkip(A,curComb,rows2Keep,cols2Keep,numColsKept-numRowsKept+x,numRowsTotal,numRowsKept,numColsKept)*binomTerm;        
        binomTerm=binomTerm*(1-numRowsKept+numColsKept+x)/(1+x)*(-1);
    }
    
    return retVal;
}

double SigmaSSkip(const double *A,size_t *curComb,const size_t *rows2Keep, const size_t *cols2Keep,const size_t r,const size_t numRowsTotal,const size_t numRowsKept,const size_t numColsKept) {
    //This adds up all of the possible values of S(A) where r columns of A
    //have been replaced by zeros. We shal choose the  cols2Keep-r columns
    //of A that are kept that are NOT zero.
    size_t combLen=numColsKept-r;
    size_t curRow,curCol;
    double retVal=0;
    
    //Initialize the current combination bufer.
    for(curRow=0;curRow<combLen;curRow++) {
        curComb[curRow]=curRow;
    }

    do {
        //We must compute the product of the row sums involving the columns
        //selected by the current combination.
        double rowSumProd=1;
        
        for(curRow=0;curRow<numRowsKept;curRow++) {
            double rowSum=0;
            size_t curRowKept=rows2Keep[curRow];
            
            for(curCol=0;curCol<combLen;curCol++) {
                rowSum+=A[curRowKept+cols2Keep[curComb[curCol]]*numRowsTotal];
            }
            rowSumProd*=rowSum;
        }
        
        retVal+=rowSumProd;
    //While we have not yet reached the final combination.
    }while(getNextComboCPP(curComb,numColsKept,combLen)==false);

    return retVal;
}


double SigmaS(const double *A,size_t *curComb,const size_t r,const size_t numRow,const size_t numCol) {
    //This adds up all of the possible values of S(A) where r columns of A
    //have been replaced by zeros. We shal choose the  numCol-r columns of
    //A that are NOT zero.
    size_t combLen=numCol-r;
    size_t curRow,curCol;
    double retVal=0;
    
    //Initialize the current combination bufer.
    for(curRow=0;curRow<combLen;curRow++) {
        curComb[curRow]=curRow;
    }

    do {
        //We must compute the product of the row sums involving the columns
        //selected by the current combination.
        double rowSumProd=1;
        
        for(curRow=0;curRow<numRow;curRow++) {
            double rowSum=0;
            for(curCol=0;curCol<combLen;curCol++) {
                rowSum+=A[curRow+curComb[curCol]*numRow];
            }
            rowSumProd*=rowSum;
        }
        
        retVal+=rowSumProd;
    //While we have not yet reached the final combination.
    }while(getNextComboCPP(curComb,numCol,combLen)==false);

    return retVal;
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
