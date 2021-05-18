/**ASSIGN2DFULLC This file contains C language functions to solve
 *              the type of 2D assignment that often arises as a
 *              subproblem subproblem when solving an S-dimensional
 *              assignment problem for target tracking, where one row and
 *              one column hold missed detection costs.
 *              The function assign2DFullC (commented below) can perform
 *              cost maximization or minimization with cost matrices having
 *              positive and/or negative elements whereas the function
 *              assign2DFullCBasic assumes that all elements of C are
 *              positive and that only minimization is performed. See the
 *              comments to the Matlab implementation of assign2DFull for
 *              more details on the algorithm.
 *
 *January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "assignAlgs2D.h"

/*string.h is needed for the memset and memcpy functions.*/
#include <string.h>
//For qsort
#include <stdlib.h>

/*We need this for INFINITY to be defined, but if a compiler does not
 *support C99, then it must be explicitly defined.*/
#include <math.h>

#ifndef INFINITY
#include<stdint.h>
const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)&infVal)
#endif

//If using *NIX or Mac OS X, isfinite functions should be defined without
//anything extra in math.h. Some Windows compilers might have problems with
//this.
#include<math.h>

static int compRows(const void *p1, const void *p2)  {
//COMPROWS This is the comparator function for qsort when it is passed a
//         2X1 array of ptrdiff_t elements that should be sorted by the
//         first row in ascending order. This is used in the function
//         assign2DFullCAlt.
    
    const ptrdiff_t val1=*(const ptrdiff_t *)p1;
    const ptrdiff_t val2=*(const ptrdiff_t *)p2;
    if(val1<val2) {
        return -1;
    } else if(val1==val2) {
        return 0;
    } else {
        return 1;
    }
}

size_t assign2DFullCBufferSize(const size_t numRowsTrue, const size_t numColsTrue) {
/**ASSIGN2DFULLCBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign2DFullC and assign2DFullCBasic algorithms.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t augmentedSize=(numRowsTrue-1)+(numColsTrue-1);
    
    return 2*augmentedSize*sizeof(ptrdiff_t)+3*augmentedSize*sizeof(size_t)+augmentedSize*sizeof(double)+augmentedSize*sizeof(bool);      
}

bool assign2DFullC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  tuples,size_t *  numTuples,void *  tempBuffer,double *  u,double *  v,const size_t numRowsTrue, const size_t numColsTrue) {
/**ASSIGN2DFULLC Solve the optimization problem 
 *      min (or max) sum_{i=0}^{numRowsTrue-1}sum_{j=0}^{numColsTrue-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numColsTrue-1}x_{i,j}=1 for i=1:(numRowsTrue-1)
 *      sum_{i=0}^{numRowsTrue-1}x_{i,j}=1 for j=1:(numColsTrue-1)
 *      x_{i,j}=0 or 1.
 *      This function can handle C with positive and negative entries.
 *      Rather than returning x, the indices of the nonzero x values
 *      are returned as tuples. A maximum of
 *      maxTuples=(numRowsTrue-1)+(numColsTrue-1)+1 tuples are possible.
 *
 *INPUTS: maximize A boolean value indicating whether maximization is
 *                 performed. False indicates that a minimization problem
 *                 is to be solved.
 *               C A pointer to the numRowsTrueXnumColsTrue matrix of
 *                 doubles, which is stored by column. Thus, the entry in
 *                 row i and column j is (i+numRowsTrue*j). This matrix
 *                 might be modified.
 *            gain A pointer to a double value that will hold the gain
 *                 (cost) of the optimal assignment.
 *          tuples A pointer to a 2XmaxTuples matrix of type ptrdiff_t,
 *                 which will hold the assigned tuples in the optimal
 *                 solution. As with C, the tuples are stored by column.
 *       numTuples A pointer to a size_t value that will indicate how many
 *                 tuples are stored in tuples.
 *      tempBuffer A pointer to a buffer of memory that is at least
 *                 assign2DFullCBufferSize(numRowsTrue, numColsTrue) bytes
 *                 in size.
 *            u, v Pointers to arrays of doubles that are each at least
 *                 augmentedSize long and can hold the dual variables.
 * numRowsTrue, numColsTrue The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The outputs are placed in gain, tuples, u, and v. The return
 *         value of this function is 0 if the problem is feasible (and
 *         gain, tuples, u, and v are valid, and it is 1 if the problem is
 *         infeasible and the values in gain, tuples, u, and v have no
 *         meaning. In either case, the values in tempBuffer and C might be
 *         modified.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */

    const size_t totalNumElsInC=numRowsTrue*numColsTrue;
    double zeroOffset;
    bool hasUnconstTuple;
    size_t i;

    /* The cost matrix must have all non-negative elements for the
     * assignment algorithm to work. Also, the algorithm only performs
     * maximization This forces all of the elements to be positive. The
     * delta is added back in to adjust the gain in the end.
     */
    if(maximize) {
        hasUnconstTuple=C[0]>0;
        
        zeroOffset=-(double)INFINITY;
        for(i=0;i<totalNumElsInC;i++) {
            if(C[i]>zeroOffset) {
                zeroOffset=C[i];
            }
        }
        
        //If C is all negative, do not shift.
        if(zeroOffset<0) {
            zeroOffset=0;
        }

        for(i=0;i<totalNumElsInC;i++) {
            C[i]=-C[i]+zeroOffset;
        }
    } else {
        hasUnconstTuple=C[0]<0;
        zeroOffset=(double)INFINITY;
        
        for(i=0;i<totalNumElsInC;i++) {
            if(C[i]<zeroOffset) {
                zeroOffset=C[i];
            }
        }
        
        //If C is all positive, do not shift.
        if(zeroOffset>0) {
            zeroOffset=0;
        }

        zeroOffset=-zeroOffset;
        for(i=0;i<totalNumElsInC;i++) {
            C[i]=C[i]+zeroOffset;
        }
    }
    
    //Run the assignment algorithm (minimization)
    (*gain)=assign2DFullCBasic(C,tuples,numTuples,tempBuffer,u,v,hasUnconstTuple,zeroOffset,numRowsTrue,numColsTrue);
    
    if((*gain)<0) {
        //If the problem is infeasible.
        return true;
    }

    //If the problem is feasible, adjust the gain to deal with having
    //adjusted the cost matrix C.
    if(maximize) {
        (*gain)=-(*gain)+zeroOffset*(*numTuples);
        for(i=0;i<numColsTrue-1;i++) {
            u[i]=-u[i];
        }
        
        for (i=0;i<numRowsTrue-1;i++) {
            v[i]=-v[i];
        }
    } else {
        (*gain)=(*gain)-zeroOffset*(*numTuples);
    }

    return false;
}

double assign2DFullCBasic(const double *C,ptrdiff_t *  tuples, size_t *  numTuples, void *tempBuffer, double *  u, double *  v,const bool hasUnconstTuple, const double zeroOffset, const size_t numRowsTrue,const size_t numColsTrue) {
/**ASSIGN2DFULLC Solve the optimization problem 
 *      min sum_{i=0}^{numRowsTrue-1}sum_{j=0}^{numColsTrue-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numColsTrue-1}x_{i,j}=1 for i=1:(numRowsTrue-1)
 *      sum_{i=0}^{numRowsTrue-1}x_{i,j}=1 for j=1:(numColsTrue-1)
 *      x_{i,j}=0 or 1.
 *      where all entries in C are >=0. Rather than returning x, the
 *      indices of the nonzero x values are returned as tuples. A maximum
 *      of
 *      maxTuples=(numRowsTrue-1)+(numColsTrue-1)+1 tuples are possible.
 *
 *INPUTS:  C A pointer to the numRowsTrueXnumColsTrue matrix of all non-
 *           negative doubles, which is stored by column. Thus, the entry
 *           in row i and column j is (i+numRowsTrue*j). This matrix will
 *           not be modified.
 *    tuples A pointer to a 2XmaxTuples matrix of type ptrdiff_t, which
 *           will hold the assigned tuples in the optimal solution. As with
 *           C, the tuples are stored by column.
 * numTuples A pointer to a size_t value that will indicate how many tuples
 *           are stored in tuples.
 * tempBuffer A pointer to a buffer of memory that is at least
 *           assign2DFullCBufferSize(numRowsTrue, numColsTrue) bytes in
 *           size.
 *      u, v Pointers to arrays of doubles that are each at least
 *           augmentedSize long and can hold the dual variables.
 * hasUnconstTuple A boolean variable indicating whether the unconstrained
 *           tuple C[0] (tuple (0,0)) should be assigned. This is present
 *           rather than just checking whether C[0]>0 to account for the
 *           possibility of transformed C matrices when once is actually
 *           performing cost maximization.
 * zeroOffset If the matrix C had originally been modified so that its
 *           elements are all positive,  (i.e. it used to have some
 *           negative elements), then this is the value that an element
 *           with a value of zero has been transformed to. This must be
 *           positive.
 * numRowsTrue, numColsTrue The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The outputs are placed in tuples, numTuples, u, and v. The gain
 *         (the cost of the optimal assignment) is returned. A return value
 *         that is >=0 indicated that the problem is feasible and the other
 *         outputs are correct. A negative value indicates that the problem
 *         is infeasible and the other outputs have no meaning.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/

    //This is the number of rows and columns in the virtually augmented
    //matrix.
    const size_t augmentedSize=(numRowsTrue-1)+(numColsTrue-1);
    ptrdiff_t *col4row, *row4col;
    size_t *pred, *Row2Scan;
    size_t *ScannedColIdx;//This holds the INDICES of the scanned columns.
    //This holds 1's and 0's for which columns were scanned.
    bool *ScannedRows;
    double *shortestPathCost;
    size_t curRow,curCol,curUnassignedCol;
            
    //Make sure that the u and v buffers are all zeros.
    memset(u,0,sizeof(double)*augmentedSize);
    memset(v,0,sizeof(double)*augmentedSize);
    
    /* These are temporary lists that will be needed when finding the
     * shortest augmenting path. All are length augmentedSize. col4row and
     * row4col are type ptrdiff_t and ScannedColIdx, Row2Scan, and pred are
     * all type size_t (to hold non-negative indices). shortestPathCost is
     * type double and ScannedRows holds boolean variables. We assume that
     * assigning the boolean variables as the last in the sequence will
     * avoid data alignment issues.
     */
    col4row=(ptrdiff_t*)(tempBuffer);
    row4col=col4row+augmentedSize;
    ScannedColIdx=(size_t*)(row4col+augmentedSize);
    Row2Scan=ScannedColIdx+augmentedSize;
    pred=Row2Scan+augmentedSize;
    shortestPathCost=(double*)(pred+augmentedSize);
    ScannedRows=(bool*)(shortestPathCost+augmentedSize);
    //ScannedRows is length augmentedSize.

    /* These will hold the indices of the assigned things. col4row will be
     * initiaized with -1 values to indicate unassigned rows. The
     * row4col array does not need to be initialized, because columns will
     * be assigned from the minimum index up, so one always knows which
     * columns are unassigned.*/
    for(curRow=0;curRow<augmentedSize;curRow++){
        col4row[curRow]=-1;
    }

    for(curUnassignedCol=0;curUnassignedCol<augmentedSize;curUnassignedCol++){
        size_t numRow2Scan,numColsScanned;
        ptrdiff_t sink;
        double delta;
/* First, find the shortest augmenting path starting at
 * curUnassignedRow.*/

        /* Mark everything as not yet scanned. A 1 will be placed in each
         * row entry as it is scanned.*/
        numColsScanned=0;
        memset(ScannedRows,0,sizeof(bool)*augmentedSize);

        for(curRow=0;curRow<augmentedSize;curRow++){
            Row2Scan[curRow]=curRow;
            /* Initially, the cost of the shortest path to each column is
             * not known and will be made infinite.*/
            shortestPathCost[curRow]=(double)INFINITY;
        }

        /*All rows need to be scanned.*/
        numRow2Scan=augmentedSize;
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
            minVal=(double)INFINITY;
            for(curRowScan=0;curRowScan<numRow2Scan;curRowScan++) {
                double reducedCost;
                
                curRow=Row2Scan[curRowScan];
                //If scanning past the end of the real rows.
                if(curRow>=numRowsTrue-1) {
                    if(curCol>=numColsTrue-1) {
                        //If we have scanned into the part of the virtually
                        //augmented matrix that is all zeros.
                        reducedCost=delta+zeroOffset-u[curCol]-v[curRow];
                    }else if(curRow==numRowsTrue-1+curCol) {
                        //If scanning a missed detection cost.
                        reducedCost=delta+C[(curCol+1)*numRowsTrue]-u[curCol]-v[curRow];
                    } else {
                        reducedCost=(double)INFINITY;
                    }
                } else {
                    //Here, we have not scanned past the end of the real
                    //rows. However, curCol might be past the end of the
                    //real columns.
                    if(curCol>=numColsTrue-1) {
                        if(curCol==numColsTrue-1+curRow) {
                            //If scanning a missed detection cost.
                            reducedCost=delta+C[curRow+1]-u[curCol]-v[curRow];
                        } else {
                            reducedCost=(double)INFINITY;
                        }
                    } else {
                        reducedCost=delta+C[(curRow+1)+(curCol+1)*numRowsTrue]-u[curCol]-v[curRow];
                    }
                }

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
            
            numRow2Scan--;//One fewer row to scan.
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
        for(curRow=0;curRow<augmentedSize;curRow++){
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

    //Next, fill in the tuples.  
    (*numTuples)=0;
    if(hasUnconstTuple) {
        tuples[0]=0;
        tuples[1]=0;
        (*numTuples)++;
    }
    
    for(curRow=0;curRow<numRowsTrue-1;curRow++) {
        const size_t idx=2*(*numTuples);
        curCol=(size_t)(col4row[curRow]+1);

        if(curCol>=numColsTrue) {
            tuples[idx]=(ptrdiff_t)(curRow+1);
            tuples[idx+1]=0;
        } else {
            tuples[idx]=(ptrdiff_t)(curRow+1);
            tuples[idx+1]=(ptrdiff_t)curCol;
        }

        (*numTuples)++;
    }

    for(curCol=0;curCol<numColsTrue-1;curCol++) {
         curRow=(size_t)(row4col[curCol]+1);

        if(curRow>=numRowsTrue) {
            size_t idx=2*(*numTuples);

            tuples[idx]=0;
            tuples[idx+1]=(ptrdiff_t)(curCol+1);

            (*numTuples)++;
        }
    }

    //Compute the gain that should be returned.
    {
        double gain=0;
        size_t curTuple;
        const size_t tupleNum=*numTuples;
        
        gain=0;
        for(curTuple=0;curTuple<tupleNum;curTuple++){
            const size_t idx=2*curTuple;
            curRow=(size_t)tuples[idx];
            curCol=(size_t)tuples[idx+1];
            
            gain+=C[curCol*numRowsTrue+curRow];
        }
        
        /* Now that everything has been assigned, free all of the temporary
         * space.*/
        return gain;
    }
}

size_t assign2DFullCAltBufferSize(const size_t numRow, const size_t numCol) {
/**ASSIGN2DFULLCALTBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign2DFullCAlt algorithm.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t numRowAug=(numRow-1)+(numCol-1);
    size_t maxRowCol;
    if(numRow<numCol) {
        maxRowCol=numCol;
    } else {
        maxRowCol=numRow;
    }
    
    return (numRowAug+maxRowCol)*sizeof(ptrdiff_t)+(maxRowCol+numRowAug*2)*sizeof(size_t)+(numRow*numCol+numRowAug)*sizeof(double)+numRowAug*sizeof(bool);
}

bool assign2DFullCAlt(const bool maximize, const double *C, double *  gain, ptrdiff_t *tuples, size_t *  numTuples,void *tempBuffer, double * u, double * v, size_t numRow, size_t numCol) {
/**ASSIGN2DFULLCALT Solve the optimization problem 
 *      min sum_{i=0}^{numRowsTrue-1}sum_{j=0}^{numColsTrue-1}C_{i,j}*x_{i,j}
 *      subject to
 *      sum_{j=0}^{numColsTrue-1}x_{i,j}=1 for i=1:(numRowsTrue-1)
 *      sum_{i=0}^{numRowsTrue-1}x_{i,j}=1 for j=1:(numColsTrue-1)
 *      x_{i,j}=0 or 1.
 *      where all entries in C are >=0. Rather than returning x, the
 *      indices of the nonzero x values are returned as tuples. A maximum
 *      of maxTuples=(numRowsTrue-1)+(numColsTrue-1)+1 tuples are possible.
 *      This implementation of the function cannot handle C matrices where
 *      both the first row and the first column of C contain non-finite
 *      elements, not counting C[0].
 *
 *INPUTS:  C A pointer to the numRowXnumCol matrix of all non-negative
 *           doubles, which is stored by column. Thus, the entry in row i
 *           and column j is (i+numRowsTrue*j). This matrix will not be
 *           modified. If both the first row and the first column of C
 *           contain non-finite elements, not counting C[0], then the
 *           problem will be labeled as infeasible, because this algorithm
 *           cannot handle such a situation.
 *    tuples A pointer to a 2XmaxTuples matrix of type ptrdiff_t, which
 *           will hold the assigned tuples in the optimal solution. As with
 *           C, the tuples are stored by column.
 * numTuples A pointer to a size_t value that will indicate how many tuples
 *           are stored in tuples.
 * tempBuffer A pointer to a buffer of memory that is at least
 *            assign2DFullCAltBufferSize(numRow,numCol) bytes in size.
 *      u, v Pointers to arrays of doubles that are both at at least
 *           (numRow-1)+(numCol-1) elements long.
 * numRow, numCol The number of rows and columns in the matrix C.
 *
 *OUTPUTS: The outputs are placed in gain, tuples, u, and v. The return
 *         value of this function is 0 if the problem is feasible (and
 *         gain, tuples, u, and v are valid, and it is 1 if the problem is
 *         infeasible and the values in gain, tuples, u, and v have no
 *         meaning. In either case, the values in tempBuffer and C might be
 *         modified. A problem will also be labeled as infeasible if C
 *         contains non finite values in both the first column and the
 *         first row, not counting C[0].
 *
 *January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/

    size_t numConstRow=numRow-1;
    size_t numConstCol=numCol-1;
    size_t curRow,curCol;
    bool hasUnconstTuple, isInfeasible, didFlip;
    size_t numAssignedRows,numUnconstRow,curTotalTuple;
    
    //To hold the modified C matrix. Memory for this is taken from
    //tempBuffer.
    double *CMod;
    //A pointer to the start of tempBuffer after the elements of CMod.
    void *restBuffer;
    
    if(maximize) {
        hasUnconstTuple=C[0]>0;
    } else {
        hasUnconstTuple=C[0]<0;
    }
    
/*The special case where a scalar is passed (a single completely
 *unconstrained  value), it is only assigned if it does not worsen the cost
 *compared to zero (nothing assigned).*/
    if(numConstCol==0&&numConstRow==0) {
        if(hasUnconstTuple) {
            tuples[0]=0;
            tuples[1]=0;
            *numTuples=1;
            *gain=C[0];
        } else {
            *numTuples=0;
            *gain=0;
        }
        //Return that the problem is feasible.
        return false;
    }

    CMod=(double*)tempBuffer;
    restBuffer=(void*)(CMod+numRow*numCol);

    /*This algorithm cannot handle non-finite terms in the first column
     *(not counting C[0]). If there are any, then transpose the matrix.
     *If both the first row and the first column have non-finite elements,
     *then the algorithm will ultimately incorrectly return that the
     *problem is not feasible.*/
    {
        size_t i;
        bool hasNonFinite=false;
        //This loop skips the unconstrained element.
        for(i=1;i<numRow;i++) {
            if(!isfinite(C[i])) {
                hasNonFinite=true;
                break;
            }
        }
        
        if(numRow>1&&hasNonFinite) {
            //Perform the transpose. The result is put into CMod.
            size_t temp;
            double *tempDPtr;

            //Copy the elements in prhs[0] into C in a transposed order.
            for(curCol=0;curCol<numCol;curCol++) {
                size_t colOffset=numRow*curCol;

                for(curRow=0;curRow<numRow;curRow++) {
                    size_t rowOffset=numCol*curRow;
                    CMod[curCol+rowOffset]=C[curRow+colOffset];                    
                }
            }
            
            temp=numCol;
            numCol=numRow;
            numRow=temp;
            
            temp=numConstRow;
            numConstRow=numConstCol;
            numConstCol=temp;
            
            tempDPtr=u;
            u=v;
            v=tempDPtr;
            
            didFlip=true;
        } else {
            //Just copy the input matrix into CMod.
            memcpy(CMod,C,numRow*numCol*sizeof(double));
            
            didFlip=false;
        }
    }

/*Subtract the costs of the unconstrained row from all of the constrained
 *rows. The cost of assigning a column to the unconstrained row thus
 *becomes 0 (even though we do not update it), and we can then solve the
 *problem using the assign2DMissedDetect algorithm. In Matlab, this would
 *be C(2:numRow,2:numCol)=bsxfun(@minus,C(2:numRow,2:numCol),C(2:numRow,1));
 *The matrix CMod is directly modified.
 */
    {
        double *CModCur=CMod+numRow;
        for(curCol=1;curCol<numCol;curCol++) {
            double *CRow=CMod;
            
            //Nothing is subtracted from the first row.
            for(curRow=1;curRow<numRow;curRow++) {
                CModCur++;
                CRow++;

                (*CModCur)=(*CModCur)-(*CRow);
            }
            //Move to the next column.
            CModCur++;
        }
    }
    
/*Call the function to perform the assignment with only one unconstrained
 *index. In Matlab, this would be
 *[tuplesRet,~,u,v]=assign2DMissedDetect(C(1:numRow,2:numCol),maximize,true);
 *Here, we have to give the proper offset to skip the first column of CMod.
 */
    if(numConstCol==0) {
        for(curRow=0;curRow<numRow-1;curRow++) {
            v[curRow]=0;
        }

    } else {
        isInfeasible=assign2DMissedDetectC(maximize, CMod+numRow, gain, tuples, restBuffer, u, v, numRow, numConstCol);

        if(isInfeasible) {
            (*numTuples)=0;
            (*gain)=-1;
            //Return that the problem is infeasible.
            return true;
        }
    }
    
    //Adjust the dual variables for the transformations that were performed
    //on C. v=v+C(2:numRow,1);
    {
        double *CCur=CMod+1;
        for(curRow=0;curRow<numRow-1;curRow++) {
            v[curRow]+=(*CCur);
            CCur++;
        }
    }
    
    /*Determine the number of rows that were assigned to something other
     *than the missed detection hypothesis.
     *numAssignedRows=sum(tuplesRet(1,:)~=1); */
    {
        ptrdiff_t *curTup=tuples;
        
        numAssignedRows=0;
        for(curCol=0;curCol<numConstCol;curCol++) {
            numAssignedRows+=((*curTup)!=0);

            curTup+=2;
        }
    }

    //The number of rows that have to be assigned to the unconstrained
    //column.
    numUnconstRow=numConstRow-numAssignedRows;
    
    (*numTuples)=numUnconstRow+numConstCol+hasUnconstTuple;
    {
        //Adjust the tuples, because the only rows assigned are not missed
        //detections. tuples(2,1:numConstCol)=tuples(2,1:numConstCol)+1;
        ptrdiff_t *curTup=tuples+1;
        size_t curTupIdx, curIdx;
        
        for(curCol=0;curCol<numConstCol;curCol++) {
            (*curTup)++;
            curTup+=2;
        }
        
    //The special case for the completely unconstrained element. It is
    //assigned if it does not worsen the cost.
        curTup=tuples+2*numConstCol;
        
        if(hasUnconstTuple) {
        //If the gain were taken above, we would add in the cost of the new
        //assignment. gain=gain+C[0];
            //This would be tuples(:,(numConstCol+1))=[1;1]; in Matlab
            curTup[0]=0;
            curTup[1]=0;

            curTotalTuple=numConstCol+1;
            curTup+=2;
        } else {
            curTotalTuple=numConstCol;
        }
        
        //Sort the tuples by the assigned row in ascending order, so that
        //we can tell which rows are not assigned.
        qsort(tuples, curTotalTuple, 2*sizeof(ptrdiff_t), compRows);
        
        //Next, all unassigned rows are assigned to the unconstrained row element
        //and are added to the end of tuples.
        
        //Skip the missed detection rows assigned to columns.
        curTup=tuples;
        curTupIdx=0;
        while(curTupIdx<*numTuples&&*curTup==0) {
            curTup+=2;
            curTupIdx++;
        }
        
        //Next, find which rows are missing and add them.
        curIdx=1;
        while(curIdx<numRow&&curTotalTuple<=*numTuples) {
            if((curTupIdx<numConstCol+hasUnconstTuple)&&((size_t)(*curTup)==curIdx)) {
                 curTup+=2;
                 curTupIdx++;
                 curIdx++;
            } else {
                ptrdiff_t *newTup=tuples+2*curTotalTuple;
                newTup[0]=(ptrdiff_t)curIdx;
                newTup[1]=0;
                 //If the gain were taken above, we would add in the cost
                 //of the new assignment. gain=gain+C(curIdx,1); in Matlab.

                 curTotalTuple++;
                 curIdx++;
            }
        }
        
        /* If a transposed array was used, then undo the flip. */
        if(didFlip==true) {
            size_t temp,i;
            //Not needed unless we flip back u and v, which is not
            //necessary here.
            //double *tempDPtr;
            ptrdiff_t *tempSTPtr;

            temp=numCol;
            numCol=numRow;
            numRow=temp;
            
            //These values were swapped, but they do not need to be swapped
            //back, because they are not going to be used again.
            //temp=numConstRow;
            //numConstRow=numConstCol;
            //numConstCol=temp;

            //tempDPtr=u;
            //u=v;
            //v=tempDPtr;

            tempSTPtr=tuples;
            for(i=0;i<*numTuples;i++) {
               const ptrdiff_t tempP=tempSTPtr[1];
               tempSTPtr[1]=tempSTPtr[0];
               tempSTPtr[0]=tempP;
               tempSTPtr+=2;
            }
        }
        
        //Now, we compute the gain from the original cost matrix.
        (*gain)=0;
        curTup=tuples;
        for(curTupIdx=0;curTupIdx<curTotalTuple;curTupIdx++) {
            (*gain)+=C[(size_t)(*curTup)+(size_t)(*(curTup+1))*numRow];
            curTup+=2;
        }
    }

    return false;
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
