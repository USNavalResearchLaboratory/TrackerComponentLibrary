/**PERM A C++ implementation of a function to calculate the matrix
*       permanent allowing for rows and columns to be skipped is desired
*       (operate on a submatrix). The permanent is equivalent to
*       calculating the determininant in the standard summation manner
*       taught in school, except all of the minus signs are replaced with
*       plus signs.
*
*INPUTS: A An mXn matrix of real doubles. If m<=n, then the standard matrix
*          permanent is found. If m>n, then the permanent of A' is found to
*          be consistent with the permanents of A and A' being equal in
*          square matrices. Empty matrices have a permanent of one by
*          definition.
* boolRowsSkip  An optional list of rows in A that should be skipped when
*          computing the matrix permanent. This is an mX1 or 1Xm boolean
*          vector where a 1 means that row should be skipped. If omitted or
*          an empty matrix is passed, no rows are skipped.
*   boolColsSkip  An optional list of columsn in A that should be skipped
*          when computing the matrix permanent. This is an nX1 or 1Xn
*          boolean vector where a 1 means that column should be skipped. If
*          omitted or an empty matrix is passed, no columns are skipped.
*
*OUTPUTS: val The matrix permanent of A.
*
* Commenbts to the algorithm are given in the comments to the Matlab
* implementation.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* val=perm(A)
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//For the swap and fill functions
#include <algorithm>
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
/*This header is required by Matlab*/
#include "mex.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t numRow, numCol;
    mxArray *AMat=NULL;
    double *A,permVal;

    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    numRow = mxGetM(prhs[0]);
    numCol = mxGetN(prhs[0]);
    
    //Empty matrices have a matrix permanent of one by definition.
    if(numRow==0||numCol==0) {
        const double retVal=1;
        plhs[0]=doubleMat2Matlab(&retVal,1,1);
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    
    //If skip lists have bee provided, then use them.
    if(nrhs>1) {
        size_t i, numRowsSkipped,numColsSkipped,numRowsKept,numColsKept;
        size_t arrayLen;
        bool *boolColsSkip, *boolRowsSkip;
        size_t *keptBuffer, *rows2Keep, *cols2Keep;
        bool didAlloc=false;
                
        if(nrhs<3||mxIsEmpty(prhs[2])) {
            boolColsSkip=reinterpret_cast<bool*>(mxMalloc(numCol*sizeof(bool)));
            //Set all of the entries to false since no column is skipped.
            std::fill_n(boolColsSkip, numCol, false);
            arrayLen=numCol;
        } else {
            boolColsSkip=copyBoolArrayFromMatlab(prhs[2], &arrayLen);
            if(arrayLen!=numCol) {
                mxFree(boolColsSkip);
                mexErrMsgTxt("The column skip list has the wrong length.");
            }
        }
        
        //Count the number of skipped columns
        numColsSkipped=0;
        for(i=0;i<numCol;i++) {
            numColsSkipped+=boolColsSkip[i];
        }
        
        if(mxIsEmpty(prhs[1])) {
            boolRowsSkip=reinterpret_cast<bool*>(mxMalloc(numRow*sizeof(bool)));
            //Set all of the entries to false since no column is skipped.
            std::fill_n(boolRowsSkip, numRow, false);
            arrayLen=numRow;
        } else {
            boolRowsSkip=copyBoolArrayFromMatlab(prhs[1], &arrayLen);
            if(arrayLen!=numRow) {
                mxFree(boolRowsSkip);
                mexErrMsgTxt("The row skip list has the wrong length.");
            }
        }
        
        //Count the number of skipped rows
        numRowsSkipped=0;
        for(i=0;i<arrayLen;i++) {
            numRowsSkipped+=boolRowsSkip[i];
        }
        
        numRowsKept=numRow-numRowsSkipped;
        numColsKept=numCol-numColsSkipped;
        
        //Empty matrices have a matrix permanent of one by definition.
        if(numRowsKept==0||numColsKept==0) {
            const double retVal=1;
            plhs[0]=doubleMat2Matlab(&retVal,1,1);
            mxFree(boolColsSkip);
            mxFree(boolRowsSkip);
            return;
        }
        
        if(numRowsKept<=numColsKept) {
            A=mxGetDoubles(prhs[0]);
        } else {
            std::swap(numRowsKept, numColsKept);
            std::swap(numRow, numCol);
            std::swap(boolRowsSkip, boolColsSkip);
            
            //This is freed using mxDestroyArray
            AMat=mxCreateDoubleMatrix(numRowsKept,numColsKept,mxREAL);
            didAlloc=true;
            mexCallMATLAB(1, &AMat, 1,  const_cast<mxArray **>(&prhs[0]), "transpose");
            A=mxGetDoubles(AMat);
        }
        
        //Set the mapping of indices of the rows in the submatrix to
        //indices in the full matrix and similarly, set the mapping of
        //columns in the submatrix to columns in the full matrix.
        //The two buffers will be allocated in one go and then freed
        //together in the end.
        keptBuffer=reinterpret_cast<size_t*>(mxMalloc(sizeof(size_t)*(numRowsKept+numColsKept)));
        rows2Keep=keptBuffer;
        cols2Keep=keptBuffer+numRowsKept;
        {//Do the mapping of the rows.
            size_t cumSkip, curIdx, curRow;
            cumSkip=0;
            curIdx=0;
            for(curRow=0;curRow<numRow;curRow++) {
                if(boolRowsSkip[curRow]==true) {
                    cumSkip++;
                } else {
                    rows2Keep[curIdx]=curIdx+cumSkip;
                    curIdx++;
                }
            }
        }
        
        {//Do the mapping of the columns.
            size_t cumSkip, curIdx, curCol;
            cumSkip=0;
            curIdx=0;
            for(curCol=0;curCol<numCol;curCol++) {
                if(boolColsSkip[curCol]==true) {
                    cumSkip++;
                } else {
                    cols2Keep[curIdx]=curIdx+cumSkip;
                    curIdx++;
                }
            } 
        }
        
        {
            size_t* buffer;
            buffer=reinterpret_cast<size_t*>(mxMalloc(numColsKept*sizeof(size_t)));
            permVal=permCPPSkip(A, numRow, rows2Keep, cols2Keep, numRowsKept, numColsKept,buffer);
            mxFree(buffer);
        }
        
        //Get rid of temporary space.
        if(didAlloc) {
            mxDestroyArray(AMat);
        }
        
        mxFree(keptBuffer);
        mxFree(boolColsSkip);
        mxFree(boolRowsSkip);
        
        //Set the return value
        plhs[0]=doubleMat2Matlab(&permVal,1,1);
        return;
    }
    
    if(numRow==numCol) {
        double *buffer;
        A=mxGetDoubles(prhs[0]);
        
        //Allocate temporary scratch space.
        buffer = reinterpret_cast<double*>(mxMalloc((2*numRow+1)*sizeof(double)));
        permVal=permSquareCPP(A,numRow,buffer);
        mxFree(buffer);
    } else {
        bool didAlloc=false;
        size_t *buffer;
        if(numRow<=numCol) {
            A=mxGetDoubles(prhs[0]);
        } else {
            std::swap(numRow, numCol); 

            //This is freed using mxDestroyArray
            AMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
            didAlloc=true;
            mexCallMATLAB(1, &AMat, 1,  const_cast<mxArray **>(&prhs[0]), "transpose");
            A=mxGetDoubles(AMat);
        }
        
        buffer=reinterpret_cast<size_t*>(mxMalloc(numCol*sizeof(size_t)));
        permVal=permCPP(A,numRow,numCol,buffer);
        mxFree(buffer);
        //Get rid of temporary space.
        if(didAlloc) {
            mxDestroyArray(AMat);
        }
    }
    
    //Set the return value
    plhs[0]=doubleMat2Matlab(&permVal,1,1);
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
