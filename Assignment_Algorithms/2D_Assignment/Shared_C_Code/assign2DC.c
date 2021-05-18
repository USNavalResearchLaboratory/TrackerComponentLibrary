/**ASSIGN2DC This file contains C language functions implementing a
 *          generalized Jonker-Volgenant shortest path assignment algorithm
 *          to solve the two-dimensional assignment problem with a
 *          rectangular cost matrix C. The function assign2DC (commented
 *          below) can perform cost maximization or minimization with cost
 *          matrices having positive and/or negative elements whereas the
 *          function assign2DCBasic assumes that all elements of C are
 *          positive and that only minimization is performed. See the
 *          comments to the Matlab implementation of assign2D for more
 *          details on the algorithm.
 *
 *January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "assignAlgs2D.h"

/* string.h is needed for the memset function.*/
#include <string.h>

/*We need this for INFINITY to be defined, but if a compiler does not
 *support C99, then it must be explicitly defined.*/
#include <math.h>

//For uint8_t and other types.
#include <stdint.h>

#ifndef INFINITY
const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)&infVal)
#endif

size_t assign2DSimpBufferSize(size_t numRow,size_t numCol) {
/**ASSIGN2SIMPBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign2DSimp algorithm.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    if(numRow<numCol) {
        size_t temp=numRow;
        numRow=numCol;
        numCol=temp;
    }
    return numRow*numCol*sizeof(double)+assign2DCBufferSize(numRow,numCol);
}

bool assign2DSimp(const bool maximize, const double *  CIn, double *  gain, ptrdiff_t *  col4RowOrig, ptrdiff_t *  row4ColOrig, void *tempBuffer, double *uCols, double *vRows, size_t numRow, size_t numCol) {
/**ASSIGN2DSIMP Solve the optimization problem 
 *      min (or max) sum_{i=0}^{numRow-1}sum_{j=0}^{numCol-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numCol-1}x_{i,j}<=1 for i=0:(numRow-1)
 *      sum_{i=0}^{numRow-1}x_{i,j}=1 for j=0:(numCol-1)
 *      x_{i,j}=0 or 1.
 *      where numCol doesn't have to be >=< numRows. This function can
 *      handle C with positive and negative entries. The values in C will
 *      not be modified. This function is simpler to use than assign2DC,
 *      because it does not modify C and it does not require that
 *      numCol<=numRow.
 *
 *INPUTS: maximize A boolean value indicating whether maximization is
 *                 performed. False indicates that a minimization problem
 *                 is to be solved.
 *               C A pointer to the numRowXnumCol matrix of doubles, which
 *                 is stored by column. Thus, the entry in row i and column
 *                 j is (i+numRow*j). This matrix will not be modified.
 *            gain A pointer to the double variable that will hold the gain
 *                 (the cost of the assignment) returned by this function.
 *                 This is the sum of the values of the assigned elements
 *                 in C.
 *         col4row A pointer to a length-numRow vector of type ptrdiff_t,
 *                 in which the result of the assignment is placed. The
 *                 entry in each element is an assignment of the element in
 *                 that row to a column. 
 *         row4col A pointer to a length-numCol vector of type ptrdiff_t,
 *                 in which the result of the assignment is placed. The
 *                 entry in each element is an assignment of the element in
 *                 that column to a row. -1 entries signify unassigned
 *                 columns.
 *      tempBuffer A pointer to a buffer of memory that is at least
 *                 assign2DSimpBufferSize(numRow,numCol) in size.
 *            u, v Pointers to arrays of doubles that hold the dual
 *                 variables. u must be at least numCol in size and v at
 *                 least numRow in size.
 *  numRow, numCol The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The result of the assignment is placed in col4row, row4col and
 *         the dual variables u and v if the problem is feasible. The
 *         return value of the function is 0 if the problem is feasible and
 *         1 if the problem is infeasible. If infeasible, the other return
 *         values do not mean anything.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/ 
    
    double *u, *v;
    uint8_t *bufferPtr=(uint8_t*)tempBuffer;
    ptrdiff_t *col4row, *row4col;
    
    double *C=(double*)bufferPtr;
    bufferPtr+=numRow*numCol*sizeof(double);

    if(numRow>=numCol) {
        //Duplicate the input array.
        memcpy(C,CIn,numRow*numCol*sizeof(double));
        
        u=uCols;
        v=vRows;
        
        col4row=col4RowOrig;
        row4col=row4ColOrig;
    } else {
        size_t temp, curRow, curCol;

        //Copy the elements in Cin into C in a transposed order.
        for(curCol=0;curCol<numCol;curCol++) {
            size_t colOffset=numRow*curCol;
            
            for(curRow=0;curRow<numRow;curRow++) {
                size_t rowOffset=numCol*curRow;
                C[curCol+rowOffset]=CIn[curRow+colOffset];
            }
        }

        temp=numCol;
        numCol=numRow;
        numRow=temp;
        
        u=vRows;
        v=uCols;
        
        col4row=row4ColOrig;
        row4col=col4RowOrig;
    }
        
    //The bufferPtr buffer must have space for
    //assign2DCBufferSize(numRow,numCol) additional elements when passed to
    //this function.
    return assign2DC(maximize, C, gain, col4row, row4col, (void*)bufferPtr, u, v, numRow, numCol);
}

size_t assign2DCBufferSize(const size_t numRow, const size_t numCol) {
/**ASSIGN2DCBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign2DC and assign2DCBasic algorithms.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    return (numCol+2*numRow)*sizeof(size_t)+numRow*sizeof(double)+numRow*sizeof(bool);
}

bool assign2DC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  col4row, ptrdiff_t *  row4col, void *tempBuffer, double *  u, double *  v, const size_t numRow, const size_t numCol) {
/**ASSIGN2DC Solve the optimization problem 
 *      min (or max) sum_{i=0}^{numRow-1}sum_{j=0}^{numCol-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numCol-1}x_{i,j}<=1 for i=0:(numRow-1)
 *      sum_{i=0}^{numRow-1}x_{i,j}=1 for j=0:(numCol-1)
 *      x_{i,j}=0 or 1.
 *      where numCol<=numRow. This function can handle C with positive and
 *      negative entries. The values in C might be modified.
 *
 *INPUTS: maximize A boolean value indicating whether maximization is
 *                 performed. False indicates that a minimization problem
 *                 is to be solved.
 *               C A pointer to the numRowXnumCol matrix of doubles, which
 *                 is stored by column. Thus, the entry in row i and column
 *                 j is (i+numRow*j). This matrix might be modified.
 *            gain A pointer to the double variable that will hold the gain
 *                 (the cost of the assignment) returned by this function.
 *                 This is the sum of the values of the assigned elements
 *                 in C.
 *         col4row A pointer to a length-numRow vector of type ptrdiff_t,
 *                 in which the result of the assignment is placed. The
 *                 entry in each element is an assignment of the element in
 *                 that row to a column. 
 *         row4col A pointer to a length-numCol vector of type ptrdiff_t,
 *                 in which the result of the assignment is placed. The
 *                 entry in each element is an assignment of the element in
 *                 that column to a row. -1 entries signify unassigned
 *                 columns.
 *      tempBuffer A pointer to a buffer of memory that is at least
 *                 assign2DCBufferSize(numRow,numCol) in size.
 *            u, v Pointers to arrays of doubles that hold the dual
 *                 variables. u must be at least numCol in size and v at
 *                 least numRow in size.
 *  numRow, numCol The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The result of the assignment is placed in col4row, row4col and
 *         the dual variables u and v if the problem is feasible. The
 *         return value of the function is 0 if the problem is feasible and
 *         1 if the problem is infeasible. If infeasible, the other return
 *         values do not mean anything.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/  

    double CDelta;
    size_t i, totalNumElsInC=numRow*numCol;
    
    /* The cost matrix must have all non-negative elements for the
     * assignment algorithm to work. This forces all of the elements to be
     * positive. The delta is added back in when computing the gain in the
     * end.*/
    if(maximize==false) {
        CDelta=(double)INFINITY;
        for(i=0;i<totalNumElsInC;i++) {
            if(C[i]<CDelta)
                CDelta=C[i];
        }

        //If C is all positive, do not shift.
        if(CDelta>0) {
            CDelta=0;
        }
        
        for(i=0;i<totalNumElsInC;i++) {
            C[i]=C[i]-CDelta;
        }
    } else {
        CDelta=-(double)INFINITY;
        for(i=0;i<totalNumElsInC;i++) {
            if(C[i]>CDelta)
                CDelta=C[i];
        }
        
        //If C is all negative, do not shift.
        if(CDelta<0) {
            CDelta=0;
        }

        for(i=0;i<totalNumElsInC;i++) {
            C[i]=-C[i]+CDelta;
        }
    }
    
    CDelta=CDelta*numCol;

    (*gain)=assign2DCBasic(C, col4row, row4col, tempBuffer, u, v, numRow, numCol);

    if((*gain)==-1) {
        //If the problem is infeasible.
        return true;
    } else {
        //The problem is feasible. Adjust for the shifting of the elements
        //in C.
        
        /* Adjust for shifting that was done to make everything positive.*/
        if(maximize==false) {
            (*gain)+=CDelta;
        } else {
            for(i=0;i<numCol;i++){
                u[i]=-u[i];
            }
            for(i=0;i<numRow;i++){
                v[i]=-v[i];
            }

            (*gain)=-(*gain)+CDelta;
        }

        return false;
    }
}

double assign2DCBasic(const double *C, ptrdiff_t * col4row, ptrdiff_t *  row4col, void *tempBuffer, double *  u, double *  v, const size_t numRow, const size_t numCol) {
/**ASSIGN2DC Solve the optimization problem 
 *      min sum_{i=0}^{numRow-1}sum_{j=0}^{numCol-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numCol-1}x_{i,j}=1 for i=0:(numRow-1)
 *      sum_{i=0}^{numRow-1}x_{i,j}<=1 for j=0:(numCol-1)
 *      x_{i,j}=0 or 1.
 *      where numCol>=numRow. The elements in C must all be >=0. This
 *      function will not modify C.
 *
 **INPUTS: C A pointer to the numRowXnumCol matrix of non-negative doubles,
 *          which is stored by column. Thus, the entry in row i and column
 *          j is (i+numRow*j). The entries in this matrix must all be >=0.
 *          This matrix will not be modified.
 *  col4row A pointer to a length-numRow vector of type ptrdiff_t, in which
 *          the result of the assignment is placed. The entry in each
 *          element is an assignment of the element in that row to a column. 
 *  row4col A pointer to a length-numCol vector of type ptrdiff_t, in which
 *          the result of the assignment is placed. The entry in each
 *          element is an assignment of the element in that column to a
 *          row. -1 entries signify unassigned columns.
 * tempBuffer A pointer to a buffer of memory that is at least
 *          assign2DCBufferSize(numRow,numCol) in size.
 *     u, v Pointers to arrays of doubles that hold the dual variables. u
 *          must be at least numCol in size and v at least numRow in size.
 *  numRow, numCol The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The assignemnt is placed in col4row and row4col. The dual
 *         variables are in u and v. The return variable is the gain, that
 *         is the cost of the assignment. If a gain of -1 is returned, that
 *         indicates that the assignment problem is not feasible.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    size_t curRow,curCol,curUnassignedCol;
    
    /* These will hold the indices of the assigned things. col4row will be
     * initiaized with -1 values to indicate unassigned rows. The
     * row4col array does not need to be initialized, because columns will
     * be assigned from the minimum index up, so one always knows which
     * columns are unassigned.*/
    for(curRow=0;curRow<numRow;curRow++){
        col4row[curRow]=-1;
    }
 
    //Make sure that the u and v buffers are all zeros.
    memset(u,0,sizeof(double)*numCol);
    memset(v,0,sizeof(double)*numRow);

    /* Assign the memory in tempBuffer to the temporary arrays that are
     * needed. We assume that assigning the boolean variables as the last
     * in the sequence will avoid data alignment issues. We use the
     *  keyword, since once assign the values are disjointly
     * accessed and we want to encourage compiler optimizations.
     */
    //This holds the INDICES of the scanned columns.
    size_t *  ScannedColIdx=(size_t*)tempBuffer;//Length numCol
    tempBuffer=(void*)((uint8_t*)tempBuffer+numCol*sizeof(size_t));
    
    size_t *  pred=(size_t*)tempBuffer;//Length numRow
    tempBuffer=(void*)((uint8_t*)tempBuffer+numRow*sizeof(size_t));
    
    size_t * Row2Scan=(size_t*)tempBuffer;//Length numRow
    tempBuffer=(void*)((uint8_t*)tempBuffer+numRow*sizeof(size_t));
    
    double *  shortestPathCost=(double*)(tempBuffer);//Length numRow
    tempBuffer=(void*)((uint8_t*)tempBuffer+numRow*sizeof(double));
    
    //This holds 1's and 0's for which rows were scanned.
    bool *  ScannedRows=(bool*)(tempBuffer);//Length numRow

    for(curUnassignedCol=0;curUnassignedCol<numCol;curUnassignedCol++){
        size_t numRow2Scan,numColsScanned;
        ptrdiff_t sink;
        double delta;
/* First, find the shortest augmenting path starting at
 * curUnassignedCol.*/

        /* Mark everything as not yet scanned. A 1 will be placed in each
         * row entry as it is scanned.*/
        numColsScanned=0;
        memset(ScannedRows,0,sizeof(bool)*numRow);

        for(curRow=0;curRow<numRow;curRow++){
            Row2Scan[curRow]=curRow;

            /* Initially, the cost of the shortest path to each row is not
             * known and will be made infinite.*/
            shortestPathCost[curRow]=(double)INFINITY;
        }

        /*All rows need to be scanned.*/
        numRow2Scan=numRow;
        /*pred will be used to keep track of the shortest path.*/
                
        /*sink will hold the final index of the shortest augmenting path.
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
            minVal=(double)INFINITY;
            for(curRowScan=0;curRowScan<numRow2Scan;curRowScan++) {
                double reducedCost;
                
                curRow=Row2Scan[curRowScan];
                
                reducedCost=delta+C[curRow+curCol*numRow]-u[curCol]-v[curRow];
                if(reducedCost<shortestPathCost[curRow]){
                    pred[curRow]=curCol;
                    shortestPathCost[curRow]=reducedCost;
                }

                //Find the minimum unassigned row that was scanned.
                if(shortestPathCost[curRow]<minVal){
                    minVal=shortestPathCost[curRow];
                    closestRowScan=curRowScan;
                }
            }

            if(minVal==(double)INFINITY) {
               /* If the minimum cost row is not finite, then the
                * problem is not feasible.*/
                return -1;
            }

            /* Change the index from the relative row index to the
             * absolute row index.*/
            closestRow=Row2Scan[closestRowScan];

            /* Add the closest row to the list of scanned rows and
             * delete it from the list of rows to scan by shifting all
             * of the items after it over by one.
             */
            ScannedRows[closestRow]=true;
            
            numRow2Scan--;//One fewer rows to scan.
            for(curRow=closestRowScan;curRow<numRow2Scan;curRow++){
                Row2Scan[curRow]=Row2Scan[curRow+1];
            }
            delta=shortestPathCost[closestRow];
            //If we have reached an unassigned column.
            if(col4row[closestRow]==-1) {
                sink=(ptrdiff_t)closestRow;
            } else{
                curCol=(size_t)col4row[closestRow];
            }
        } while(sink==-1);
        
/* Next, update the dual variables.*/
        //Update the first column in the augmenting path.
        u[curUnassignedCol]+=delta;

        //Update the rest of the columns in the augmenting path.
        //curCol starts from 1, not zero, so that it skips
        //curUnassignedCol.
        for(curCol=1;curCol<numColsScanned;curCol++) {
            size_t curScannedIdx=ScannedColIdx[curCol];
            u[curScannedIdx]+=delta-shortestPathCost[row4col[curScannedIdx]];
        }

        //Update the rows in the augmenting path.
        for(curRow=0;curRow<numRow;curRow++){
            if(ScannedRows[curRow]==true){
                v[curRow]+=-delta+shortestPathCost[curRow];
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
            gain+=C[curCol*numRow+(size_t)row4col[curCol]];
        }
        
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
