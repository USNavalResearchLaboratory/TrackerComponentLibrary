/**ASSIGN2DMISSEDDETECTC This file contains C language functions to solve
 *              the two-dimensional assignment problem common to target
 *              tracking, where, in a tracking context, one row represents
 *              missed detection costs. The function assign2DMissedDetectC
 *              (commented below) can perform cost maximization or
 *              minimization with cost matrices having positive and/or
 *              negative elements whereas the function
 *              assign2DCMissedDetectBasic assumes that all elements of C
 *              are positive and that only minimization is performed. See
 *              the comments to the Matlab implementation of
 *              assign2DMissedDetect for more details on the algorithm.
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

#ifndef INFINITY
#include<stdint.h>
const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)&infVal)
#endif

size_t assign2DMissedDetectCBufferSize(const size_t numRowsTrue, const size_t numCol) {
/**ASSIGN2DMISSEDDETECTCBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign2DMissedDetectC and assign2DCMissedDetectBasic
 *      algorithms.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t numRow=numRowsTrue-1+numCol; 
    return (numRow+numCol)*sizeof(ptrdiff_t)+(numCol+numRow*2)*sizeof(size_t)+numRow*sizeof(double)+numRow*sizeof(bool);
}

bool assign2DMissedDetectC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  tuples, void *  tempBuffer, double *  u, double *  v, const size_t numRowsTrue, const size_t numCol) {
/**ASSIGN2DMISSEDDETECT Solve the optimization problem 
 *         minimize (or maximize)
 *         sum_{i=0}^{numMeas}sum_{j=1}^{numTar}C_{i,j}*x_{i,j}
 *         subject to
 *         sum_{i=0}^{numMeas}x_{i,j}=1 for all j
 *         sum_{j=0}^{numTar-1}x_{i,j}<=1 for i=1:numMeas
 *         x_{i,j}=0 or 1. 
 *         where the index i=0 has the costs of assigning a target to a
 *         missed detection, so multiple targets can be assigned to a
 *         missed detection hypothesis. Measurements are indexed starting
 *         at 1. Let numCol=numTar and numRowsTrue=numMeas+1. This function
 *         can handle C with positive and negative entries. The values in C
 *         might be modified.
 *
 *%INPUTS: maximize A boolean value indicating whether maximization is
 *                 performed. False indicates that a minimization problem
 *                 is to be solved.
 *               C A pointer to a numRowsTrueXnumCol cost matrix of doubles
 *                 stored by column, where the first row is the cost of a
 *                 missed detection for each target and the subsequent rows
 *                 are the costs of assigning a measurement to a target.
 *                 This matrix might be modified.
 *            gain A pointer to the double variable that will hold the gain
 *                 (the cost of the assignment) returned by this function.
 *                 This is the sum of the values of the assigned elements
 *                 in C.
 *          tuples A pointer to a 2XnumCol matrix of type size_t,
 *                 which will hold the assigned tuples in the optimal
 *                 solution. As with C, the tuples are stored by column.
 *      tempBuffer A pointer to a buffer of memory that is at least
 *                 assign2DMissedDetectCBufferSize(numRowsTrue,numCol) in
 *                 size,
 *            u, v Pointers to arrays of doubles that hold the dual
 *                 variables. u must be at least numCol in size and v at
 *                 least numRowsTrue-1+numCol in size.
 *     numRowsTrue The number of the rows of the cost matrix, as defined
 *                 above.
 *          numCol The number of columns of the cost matrix.
 *
 *OUTPUTS: The outputs are placed in gain, tuples, u, and v. The return
 *         value of this function is 0 if the problem is feasible (and
 *         gain, tuples, u, and v are valid, and it is 1 if the problem is
 *         infeasible and the values in gain, tuples, u, and v have no
 *         meaning. In either case, the values in tempBuffer and C might be
 *         modified.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    double CDelta;
    size_t i;
    const size_t totalNumElsInC=numRowsTrue*numCol;
    //The number of rows in the virtually augmented matrix.
    const size_t numRow=numRowsTrue-1+numCol;

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
    
    if(numCol==0||numRowsTrue==0) {
        //If the matrix is empty, the problem is infeasible.
        *gain=-1;
    } else {
        (*gain)=assign2DCMissedDetectBasic(C,tuples,tempBuffer,u,v,numRowsTrue,numCol);
    }

    if((*gain)==-1) {
        //If the problem is infeasible.
        return true;
    } else if(maximize) {
        (*gain)=-(*gain)+CDelta*numCol;

        for(i=0;i<numCol;i++){
            u[i]=-u[i];
        }
        
        for(i=0;i<numRow;i++){
            v[i]=-v[i];
        }
    } else {
        (*gain)+=CDelta*numCol;
    }

    return false;    
}

double assign2DCMissedDetectBasic(const double *C, ptrdiff_t *  tuples, void *tempBuffer, double *  u, double *  v, const size_t numRowsTrue, const size_t numCol) {
/**ASSIGN2DMISSEDDETECT Solve the optimization problem 
 *         minimize
 *         sum_{i=0}^{numMeas}sum_{j=1}^{numTar}C_{i,j}*x_{i,j}
 *         subject to
 *         sum_{i=0}^{numMeas}x_{i,j}=1 for all j
 *         sum_{j=0}^{numTar-1}x_{i,j}<=1 for i=1:numMeas
 *         x_{i,j}=0 or 1. 
 *         where the index i=0 has the costs of assigning a target to a
 *         missed detection, so multiple targets can be assigned to a
 *         missed detection hypothesis. Measurements are indexed starting
 *         at 1. Let numCol=numTar and numRowsTrue=numMeas+1. All entries
 *         in C must be non-negative.
 *
 *%INPUTS: C A pointer to a numRowsTrueXnumCol cost matrix of all non-
 *           negative doubles stored by column, where the first row is the
 *           cost of a missed detection for each target and the subsequent
 *           rows are the costs of assigning a measurement to a target.
 *           This matrix will not be be modified.
 *    tuples A pointer to a 2XnumCol matrix of type size_t, which will hold
 *           the assigned tuples in the optimal solution. As with C, the
 *           tuples are stored by column.
 * tempBuffer A pointer to a buffer of memory that is at least
 *           assign2DMissedDetectCBufferSize(numRowsTrue,numCol) in size,
 *      u, v Pointers to arrays of doubles that hold the dual variables. u
 *           must be at least numCol in size and v at least
 *           numRowsTrue-1+numCol in size.
 * numRowsTrue The number of the rows of the cost matrix, as defined above.
 *    numCol The number of columns of the cost matrix.
 *
 *OUTPUTS: The outputs are placed in tuples, u, and v. The return variable
 *         is the gain, that is the cost of the assignment. If a gain of
 *         -1 is returned, that indicates that the assignment problem is
 *         not feasible.
 *
 *December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    const size_t numRow=numRowsTrue-1+numCol;
    size_t curRow,curCol,curUnassignedCol;

    //Make sure that the u and v buffers are all zeros.
    memset(u,0,sizeof(double)*numCol);
    memset(v,0,sizeof(double)*numRow);

    /* Assign the memory in tempBuffer to the temporary arrays that are
     * needed. We assume that assigning the boolean variables as the last
     * in the sequence will avoid data alignment issues. 
     */
    ptrdiff_t *  col4row=(ptrdiff_t*)tempBuffer;
    tempBuffer=(void*)((ptrdiff_t*)tempBuffer+numRow);
    
    ptrdiff_t *  row4col=(ptrdiff_t *)tempBuffer;
    tempBuffer=(void*)((ptrdiff_t*)tempBuffer+numCol);
    
    //This holds the INDICES of the scanned columns.
    size_t *  ScannedColIdx=(size_t*)(tempBuffer);
    tempBuffer=(void*)((size_t*)tempBuffer+numCol);
    
    size_t *  pred=(size_t*)tempBuffer;
    tempBuffer=(void*)((size_t*)tempBuffer+numRow);
    
    size_t *  Row2Scan=(size_t*)tempBuffer;
    tempBuffer=(void*)((size_t*)tempBuffer+numRow);
    
    double *  shortestPathCost=(double*)(tempBuffer);
    tempBuffer=(void*)((double*)tempBuffer+numRow);
    
    //This holds 1's and 0's for which rows were scanned.
    bool *  ScannedRows=(bool*)(tempBuffer);
    //ScannedRows is length numRow.
    
    /* These will hold the indices of the assigned things. col4row will be
     * initiaized with -1 values to indicate unassigned rows. The
     * row4col array does not need to be initialized, because columns will
     * be assigned from the minimum index up, so one always knows which
     * columns are unassigned.*/
    for(curRow=0;curRow<numRow;curRow++){
        col4row[curRow]=-1;
    }

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
                
                if(curRow>=numRowsTrue-1) {
                    //If scanning a missed detection cost.
                    if(curRow==numRowsTrue-1+curCol) {
                        reducedCost=delta+C[curCol*numRowsTrue]-u[curCol]-v[curRow];
                    } else {
                        reducedCost=(double)INFINITY;
                    }
                } else {
                    reducedCost=delta+C[curRow+1+curCol*numRowsTrue]-u[curCol]-v[curRow];
                }

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

            if(minVal==(double)INFINITY) {
               /* If the minimum cost column is not finite, then the
                * problem is not feasible.*/
                return -1;
            }

            /* Change the index from the relative column index to the
             * absolute column index.*/
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

            //If we have reached an unassigned row.
            if(col4row[closestRow]==-1) {
                sink=(ptrdiff_t)closestRow;
            } else{
                curCol=(size_t)col4row[closestRow];
            }
        } while(sink==-1);
        
/* Next, update the dual variables.*/
        //Update the first row in the augmenting path.
        u[curUnassignedCol]+=delta;

        //Update the rest of the rows in the augmenting path.
        //curRow starts from 1, not zero, so that it skips curUnassignedRow.
        for(curCol=1;curCol<numColsScanned;curCol++) {
            size_t curScannedIdx=ScannedColIdx[curCol];
            u[curScannedIdx]+=delta-shortestPathCost[row4col[curScannedIdx]];
        }

        //Update the columns in the augmenting path.
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
    
    {
        double gain=0;
    
        for(curCol=0;curCol<numCol;curCol++) {
            const size_t idx=2*curCol;
            size_t rowIdx=(size_t)row4col[curCol];

            if(rowIdx>=numRowsTrue-1) {
                tuples[idx]=0;
                
                //If the assignment is to the virtually augmented part of
                //the cost matrix (a missed detection).
                gain+=C[curCol*numRowsTrue];
            } else {
                tuples[idx]=(ptrdiff_t)(rowIdx+1);
                
                //If the assignment is not a missed detection.
                gain+=C[curCol*numRowsTrue+rowIdx+1];
            }
            tuples[idx+1]=(ptrdiff_t)curCol;
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
