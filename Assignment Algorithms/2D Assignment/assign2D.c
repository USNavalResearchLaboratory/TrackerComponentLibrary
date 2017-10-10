/**ASSIGN2DALT A C-code (for Matlab) implementation of the shortest path
 *             assignment algorithm to solve the two-dimensional assignment
 *             problem with a rectangular cost matrix C. This
 *             implementation scans the cost matrix by row rather than by
 *             column.
 *
 *INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
 *          where the largest finite element minus the smallest element is
 *          a finite quantity (does not overflow) when performing
 *          minimization and where the smallest finite element minus the
 *          largest element is finite when performing maximization. 
 *          Forbidden assignments can be given costs of +Inf for
 *          minimization and -Inf for maximization.
 * maximize If true, the minimization problem is transformed into a
 *          maximization problem. The default if this parameter is omitted
 *          is false.
 *
 *OUTPUTS:  col4row     A numRowX1 Matlab vector where the entry in each
 *                      element is an assignment of the element in that row
 *                      to a column. 0 entries signify unassigned rows.
 *          row4col     A numColX1 vector where the entry in each element
 *                      is an assignment of the element in that column to a
 *                      row. 0 entries signify unassigned columns.
 *          gain        The sum of the values of the assigned elements in
 *                      C.
 *          u           The dual variable for the columns.
 *          v           The dual variable for the rows.
 *
 *DEPENDENCIES: string.h
 *              math.h
 *              stdint.h if math.h does not define INFINITY
 *              stdlib.h
 *              matrix.h
 *              mex.h
 *              MexValidation.h
 *
 * Note that the dual variables produced by a shortest path assignment
 * algorithm that scans by row are not interchangeable with those of a
 * shortest path assignment algorithm that scans by column. Matlab stores
 * matrices row-wise.
 *
 * The algorithm is described in more detail in the comments to the Matlab
 * implementation.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * [col4row,row4col,gain,u,v]=assign2DAlt(C,maximize)
 *
 * October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/* string.h is needed for the memset function.*/
#include <string.h>
/*We need this for INFINITY to be defined, but if a compiler does not
 *support C99, then it must be explicitly defined.*/
#include <math.h>
#ifndef INFINITY
#include<stdint.h>
const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)&infVal)
#endif
/*stdlib.h and matrix.h are needed for mxCalloc.*/
#include <stdlib.h>
#include "matrix.h"
/*This header is required by Matlab.*/
#include "mex.h"
//This makes sure that the bool type is defined and provides functions for
//handling Matlab data.
#include "MexValidation.h"

double shortestPathCFast(double *C,ptrdiff_t *col4rowPtr,ptrdiff_t *row4colPtr,double *u, double *v,size_t numRow, size_t numCol);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double gain, *C, CDelta;
    size_t numRow, numCol, i;
    ptrdiff_t *col4row, *row4col;
    mxArray *col4rowMATLAB, *row4colMATLAB, *uMATLAB, *vMATLAB, *CMat;
    bool didFlip=false;
    bool maximize=false;
    
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs==2){
        maximize=getBoolFromMatlab(prhs[1]);
    }
    
    if(nrhs>2) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>5) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    /*Verify the validity of the assignment matrix.*/
    checkRealDoubleArray(prhs[0]);
    
    /* Get the dimensions of the input data and the pointer to the matrix.
     * It is assumed that the matrix is not so large in M or N as to cause
     * an overflow when using a SIGNED integer data type.*/
    numRow = mxGetM(prhs[0]);
    numCol = mxGetN(prhs[0]);
    
    /* Transpose the matrix, if necessary, so that the number of rows is
     * >= the number of columns.*/
    if(numRow>=numCol) {
        CMat=mxDuplicateArray(prhs[0]);
    } else {
        size_t temp;
        
        temp=numCol;
        numCol=numRow;
        numRow=temp;
        
        CMat=mxCreateDoubleMatrix(numCol,numRow,mxREAL);
        mexCallMATLAB(1, &CMat, 1, (mxArray **)&prhs[0], "transpose");
        didFlip=true;
    }
    
    C = (double*)mxGetData(CMat);
    
    /* The cost matrix must have all non-negative elements for the
     * assignment algorithm to work. This forces all of the elements to be
     * positive. The delta is added back in when computing the gain in the
     * end.*/
    if(maximize==false) {
        CDelta=INFINITY;
        for(i=0;i<numRow*numCol;i++) {
            if(C[i]<CDelta)
                CDelta=C[i];
        }

        for(i=0;i<numRow*numCol;i++) {
            C[i]=C[i]-CDelta;
        }
    } else {
        CDelta=-INFINITY;
        for(i=0;i<numRow*numCol;i++) {
            if(C[i]>CDelta)
                CDelta=C[i];
        }

        for(i=0;i<numRow*numCol;i++) {
            C[i]=-C[i]+CDelta;
        }
    }
    
    CDelta=CDelta*numCol;
    
   //Allocate space for the return variables to Matlab.
    
    col4rowMATLAB= allocSignedSizeMatInMatlab(numRow, 1);
    row4colMATLAB= allocSignedSizeMatInMatlab(numCol, 1);
    
    /* These will hold the dual variable values. They are all initialized
     * to zero.*/
    uMATLAB = mxCreateNumericMatrix(numCol,1,mxDOUBLE_CLASS,mxREAL);
    vMATLAB = mxCreateNumericMatrix(numRow,1,mxDOUBLE_CLASS,mxREAL);
    
    col4row=(ptrdiff_t*)mxGetData(col4rowMATLAB);
    row4col=(ptrdiff_t*)mxGetData(row4colMATLAB);
    
    gain=shortestPathCFast(C,col4row,row4col, (double*)mxGetData(uMATLAB), (double*)mxGetData(vMATLAB), numRow, numCol);
    
    mxDestroyArray(CMat);
    /* If a transposed array was used */
    if(didFlip==true) {
        size_t temp;
        ptrdiff_t *tempIPtr;
        mxArray *tempMat;

        temp=numCol;
        numCol=numRow;
        numRow=temp;
        
        tempIPtr=row4col;
        row4col=col4row;
        col4row=tempIPtr;
        
        tempMat=row4colMATLAB;
        row4colMATLAB=col4rowMATLAB;
        col4rowMATLAB=tempMat;
        
        tempMat=uMATLAB;
        uMATLAB=vMATLAB;
        vMATLAB=tempMat;
    }
    
    //If there is no feasible solution.
    if(gain==-1){
        mxDestroyArray(col4rowMATLAB);
        mxDestroyArray(row4colMATLAB);
        
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        if(nlhs>1){
            plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);
            
            if(nlhs>2){
                plhs[2]=mxCreateDoubleScalar(gain);
                
                if(nlhs>3) {
                    plhs[3]=uMATLAB;
                    
                    if(nlhs>4){
                        plhs[4]=vMATLAB;
                    }
                }
                
            }
        }
        return;
    } else {
        
        /* Adjust for shifting that was done to make everything positive.*/
        if(maximize==false) {
            gain=gain+CDelta;

        } else {
            gain=-gain+CDelta;
        }
        
        /* Convert the C indices into indices for MATLAB.*/
        for(i=0;i<numRow;i++){
            col4row[i]=col4row[i]+1;
        }
                
        plhs[0]=col4rowMATLAB;
        if(nlhs>1){
            plhs[1]=row4colMATLAB;
            
            for(i=0;i<numCol;i++){
                row4col[i]=row4col[i]+1;
            }
            
            if(nlhs>2){
                plhs[2]=mxCreateDoubleScalar(gain);
                
                if(nlhs>3) {
                    plhs[3]=uMATLAB;
                    
                    if(nlhs>4){
                        plhs[4]=vMATLAB;
                    }
                }
            }
        }  
    }
    
    /*Free the u and v matrices if they are not returned.*/
    if(nlhs<4){
        mxDestroyArray(uMATLAB);
    }
    if(nlhs<5){
        mxDestroyArray(vMATLAB);   
    }
    
    return;
}

double shortestPathCFast(double *C, ptrdiff_t *col4row, ptrdiff_t *row4col, double *u, double *v, size_t numRow, size_t numCol) {
    double *shortestPathCost;
    size_t curRow,curCol,curUnassignedCol, *pred, *Row2Scan;
    bool *ScannedRows;//This holds 1's and 0's for which columns were scanned.
    size_t* ScannedColIdx;//This holds the INDICES of the scanned rows.
    
    /* These will hold the indices of the assigned things. row4col will be
     * initiaized with -1 values to indicate unassigned columns. The
     * col4row array does not need to be initialized, because rows will be
     * assigned from the minimum index up, so one always knows which rows
     * are unassigned.*/
    for(curRow=0;curRow<numRow;curRow++){
        col4row[curRow]=-1;
    }
 
    /* These are temporary lists that will be needed when finding the
     * shortest augmenting path.
     */
    ScannedColIdx=(size_t*)mxMalloc(numCol*sizeof(size_t));
    ScannedRows=(bool*)mxMalloc(numRow*sizeof(bool));
    Row2Scan=(size_t*)mxMalloc(numRow*sizeof(size_t));
    shortestPathCost=(double*)mxMalloc(numRow*sizeof(double));
    pred=(size_t*)mxMalloc(numRow*sizeof(size_t));
    
    for(curUnassignedCol=0;curUnassignedCol<numCol;curUnassignedCol++){
        size_t numRow2Scan,numColsScanned;
        ptrdiff_t sink;
        double delta;
/* First, find the shortest augmenting path starting at
 * curUnassignedRow.*/

        /* Mark everything as not yet scanned. A 1 will be placed in each
         * column entry as it is scanned.*/
        numColsScanned=0;
        ScannedRows=(bool*)memset(ScannedRows,0,sizeof(bool)*numRow);

        for(curRow=0;curRow<numRow;curRow++){
            Row2Scan[curRow]=curRow;
            /* Initially, the cost of the shortest path to each column is not
             * known and will be made infinite.*/
            shortestPathCost[curRow]=INFINITY;
        }

        /*All columns need to be scanned.*/
        numRow2Scan=numRow;
        /*pred will be used to keep track of the shortest path.*/
                
        /*Sink will hold the final index of the shortest augmenting path.
         *If the problem is not feasible, then sink will remain -1.*/
        sink=-1;
        delta=0;
        curCol=curUnassignedCol;

        do {
            double minVal;
            //The initialization is just to silence a warning if compiling
            //using -Wconditional-uninitialized.
            size_t curRowScan,closestRow,closestRowScan=0;
            /*Mark the current column as having been visited.*/
            ScannedColIdx[numColsScanned]=curCol;
            numColsScanned++;
            
            /*Scan all of the rows that have not already been scanned.*/
            minVal=INFINITY;
            for(curRowScan=0;curRowScan<numRow2Scan;curRowScan++) {
                double reducedCost;
                
                curRow=Row2Scan[curRowScan];
                
                reducedCost=delta+C[curRow+curCol*numRow]-u[curCol]-v[curRow];
                if(reducedCost<shortestPathCost[curRow]){
                    pred[curRow]=curCol;
                    shortestPathCost[curRow]=reducedCost;
                }

                //Find the minimum unassigned column that was scanned.
                if(shortestPathCost[curRow]<minVal){
                    minVal=shortestPathCost[curRow];
                    closestRowScan=curRowScan;
                }
            }

            if(minVal==INFINITY) {
               /* If the minimum cost column is not finite, then the
                * problem is not feasible.*/
                mxFree(ScannedColIdx);
                mxFree(ScannedRows);
                mxFree(pred);
                mxFree(Row2Scan);
                mxFree(shortestPathCost);
                return -1;
            }

            /* Change the index from the relative column index to the
             * absolute column index.*/
            closestRow=Row2Scan[closestRowScan];

            /* Add the closest column to the list of scanned columns and
             * delete it from the list of columns to scan by shifting all
             * of the items after it over by one.
             */
            ScannedRows[closestRow]=true;
            
            numRow2Scan--;//One fewer column to scan.
            for(curRow=closestRowScan;curRow<numRow2Scan;curRow++){
                Row2Scan[curRow]=Row2Scan[curRow+1];
            }
            delta=shortestPathCost[closestRow];
            //If we have reached an unassigned row.
            if(col4row[closestRow]==-1) {
                sink=(ptrdiff_t)closestRow;
            } else{
                curCol=(size_t)col4row[closestRow];
            }
        } while(sink==-1);
        
/* Next, update the dual variables.*/
        //Update the first row in the augmenting path.
        u[curUnassignedCol]=u[curUnassignedCol]+delta;

        //Update the rest of the rows in the augmenting path.
        //curRow starts from 1, not zero, so that it skips curUnassignedRow.
        for(curCol=1;curCol<numColsScanned;curCol++) {
            size_t curScannedIdx=ScannedColIdx[curCol];
            u[curScannedIdx]=u[curScannedIdx]+delta-shortestPathCost[row4col[curScannedIdx]];
        }

        //Update the columns in the augmenting path.
        for(curRow=0;curRow<numRow;curRow++){
            if(ScannedRows[curRow]==true){
                v[curRow]=v[curRow]-delta+shortestPathCost[curRow];
            }
        }

//Remove the current node from those that must be assigned.
        curRow=(size_t)sink;
        do{
            ptrdiff_t h;
            curCol=pred[curRow];
            col4row[curRow]=(ptrdiff_t)curCol;
            h=row4col[curCol];
            row4col[curCol]=(ptrdiff_t)curRow;
            curRow=(size_t)h;
        } while(curCol!=curUnassignedCol);
    }

    //Determine the gain to return
    {
        double gain=0;
        for(curCol=0;curCol<numCol;curCol++){
           gain=gain+C[curCol*numRow+(size_t)row4col[curCol]];
        }
        
        /* Now that everything has been assigned, free all of the temporary
         * lists.*/
        mxFree(ScannedColIdx);
        mxFree(ScannedRows);
        mxFree(pred);
        mxFree(Row2Scan);
        mxFree(shortestPathCost);
        return gain;
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
