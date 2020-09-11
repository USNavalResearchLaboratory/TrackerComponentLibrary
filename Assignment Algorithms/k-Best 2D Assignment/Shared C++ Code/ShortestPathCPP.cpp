/**ShortestPathCPP  Functions in C++ implementing the shortest augmenting
 *                  path algorithms to get the best or the k-best
 *                  hypotheses for a 2D assignment problem. The code is
 *                  easiest to understand after examining the Matlab
 *                  implementation. 
 *
 *  This file relies on the file ShortestPathCPP.hpp. Much of the
 *  documentation for the functions is found in that header.
 *
 *November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "ShortestPathCPP.hpp"
#include <stdexcept>
#include <queue>
/*This is needed for memcpy*/
#include <cstring>
#include <algorithm>
#include <limits>
//Needed for isfinite
#include <math.h>
//Needed for bsearch and qsort.
#include <cstdlib>

using namespace std;

struct pMurtyHyp {
/* This structure is defined so that pointers to MurtyHyp classes can be
 * placed in the standard C++ priority queue structure.*/
	MurtyHyp *ptr;

    inline bool operator< (const pMurtyHyp &other) const {
		return ptr->gain > other.ptr->gain;
	}
    
	inline pMurtyHyp(MurtyHyp *in) {
		ptr = in;
	}
};

//Prototypes for functions used in this file that are not present in
//the header ShortestPathCPP.hpp.
void calcGain(MurtyHyp *problemSol,const ScratchSpace &workMem,const size_t numRow,const size_t numCol4Gain);
void updateDualAndAugment(MurtyHyp *problemSol,const ScratchSpace& workMem,const size_t curUnassignedCol, const size_t numColsScanned,const size_t numDim,const ptrdiff_t sink,const double delta);
void split(MurtyHyp *parentHyp,priority_queue<pMurtyHyp> &HypQueue, ScratchSpace &workMem, const size_t numVarCol,const size_t numDim);
double makeCostMatrixSafe(double *CMod,const double *COrig,const size_t numEl, const bool maximize);
MurtyHyp *shortestPathUpdateCPP(const MurtyHyp *parentHyp, ScratchSpace &workMem,const size_t curUnassignedCol, size_t numRow2Scan, const size_t numVarCol, const size_t numDim);


inline int compare (const void * a, const void * b) {
/*This comparison function is needed to use the qsort and bsearch functions
 *that are part of the standard C library.*/
  return static_cast<int>(*reinterpret_cast<const ptrdiff_t*>(a) - *reinterpret_cast<const ptrdiff_t*>(b));
}

void calcGain(MurtyHyp *problemSol,const ScratchSpace &workMem,const size_t numRow,const size_t numCol4Gain) {
/*CALCGAIN: Compute the cost of a particular assignment specified by
 *          problemSol.
 *
 *INPUTS: problemSol  A pointer to a MurtyHyp in which the gain result will
 *                    be placed and that holds a solution.
 *        workMem     An instance of the ScratchSpace class that holds a
 *                    pointer to the cost matrix C.
 *
 *OUTPUTS: Nothing is returned; the cost (gain) of the assignment specified
 *         in problemSol is placed in problemSol->gain.
 **/
    
    size_t curCol;
    double gain=0;
    
    for(curCol=0;curCol<numCol4Gain;curCol++){
       gain=gain+workMem.C[curCol*numRow+static_cast<size_t>(problemSol->row4col[curCol])];
    }

    problemSol->gain=gain;
}

void updateDualAndAugment(MurtyHyp *problemSol,const ScratchSpace& workMem,const size_t curUnassignedCol, const size_t numColsScanned,const size_t numDim,const ptrdiff_t sink,const double delta){
/*UPDATEDUALANDAUGMENT Update the dual random variables when computing the
 *               shortest path or when updating a shortest path hypothesis.
 *               This also removes the sink row from those that need to be
 *               assigned.
 */
    
    size_t curRow, curCol;

    //Update the first row in the augmenting path.
    problemSol->u[curUnassignedCol]=problemSol->u[curUnassignedCol]+delta;
    
    //Update the rest of the rows in the augmenting path.
    //curRow starts from 1, not zero, so that it skips curUnassignedRow.
    for(curCol=1;curCol<numColsScanned;curCol++) {
        size_t curScannedIdx=workMem.ScannedColIdx[curCol];
        problemSol->u[curScannedIdx]=problemSol->u[curScannedIdx]+delta-workMem.shortestPathCost[problemSol->row4col[curScannedIdx]];
    }
    
    //Update the columns in the augmenting path.
    for(curRow=0;curRow<numDim;curRow++){
        if(workMem.ScannedRows[curRow]==true){
            problemSol->v[curRow]=problemSol->v[curRow]-delta+workMem.shortestPathCost[curRow];
        }
    }
//Remove the current node from those that must be assigned.
    curRow=static_cast<size_t>(sink);
    do{
        ptrdiff_t h;
        curCol=workMem.pred[curRow];
        problemSol->col4row[curRow]=static_cast<ptrdiff_t>(curCol);
        h=problemSol->row4col[curCol];
        problemSol->row4col[curCol]=static_cast<ptrdiff_t>(curRow);
        curRow=static_cast<size_t>(h);
    } while(curCol!=curUnassignedCol);
}

int shortestPathCPP(MurtyHyp *problemSol,ScratchSpace &workMem,const size_t numRow, const size_t numCol, const size_t numCol4Gain) {
/*SHORTESTPATHCPP A C++ implementation of the basic shortest augmenting
 *                path 2D assignment algorithm.
 **/
    size_t curRow, curCol, curUnassignedCol;
    
    /* These will hold the indices of the assigned things. row4col will be
     * initiaized with -1 values to indicate unassigned columns. The
     * col4row array does not need to be initialized, because rows will be
     * assigned from the minimum index up, so one always knows which rows
     * are unassigned.*/
    fill_n(problemSol->col4row,numRow,-1);
    
    /* These will hold the dual variable values. They are all initialized
     * to zero.*/
    fill_n(problemSol->u,numCol,0);
    fill_n(problemSol->v,numRow,0);
    problemSol->activeCol=0;
    fill_n(problemSol->forbiddenActiveRows,numRow,false);
 
    for(curUnassignedCol=0;curUnassignedCol<numCol;curUnassignedCol++){
        size_t numRow2Scan,numColsScanned;
        ptrdiff_t sink;
        double delta;
/* First, find the shortest augmenting path starting at
 * curUnassignedCol.*/
        
        /* Mark everything as not yet scanned. A 1 will be placed in each
         * row entry as it is scanned.*/
        numColsScanned=0;
        fill_n(workMem.ScannedRows,numRow,false);
        /* Initially, the cost of the shortest path to each column is not
         * known and will be made infinite.*/
        fill_n(workMem.shortestPathCost,numRow,numeric_limits<double>::infinity());
        
        /*All rows need to be scanned.*/
        for(curRow=0;curRow<numRow;curRow++){
            workMem.Row2Scan[curRow]=static_cast<ptrdiff_t>(curRow);
        }

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
            workMem.ScannedColIdx[numColsScanned]=curCol;
            numColsScanned++;
            
            /*Scan all of the columns that have not already been scanned.*/
            minVal=numeric_limits<double>::infinity();
            for(curRowScan=0;curRowScan<numRow2Scan;curRowScan++) {
                double reducedCost;
                
                curRow=static_cast<size_t>(workMem.Row2Scan[curRowScan]);
                reducedCost=delta+workMem.C[curRow+curCol*numRow]-problemSol->u[curCol]-problemSol->v[curRow];
                
                if(reducedCost<workMem.shortestPathCost[curRow]){
                    workMem.pred[curRow]=curCol;
                    workMem.shortestPathCost[curRow]=reducedCost;
                }

                //Find the minimum unassigned column that was scanned.
                if(workMem.shortestPathCost[curRow]<minVal){
                    minVal=workMem.shortestPathCost[curRow];
                    closestRowScan=curRowScan;
                }
            }
            
            if(minVal==numeric_limits<double>::infinity()) {
               /* If the minimum cost column is not finite, then the
                * problem is not feasible.*/
                problemSol->gain=-1;
                return 1;
            }
            
            /* Change the index from the relative column index to the
             * absolute column index.*/
            closestRow=static_cast<size_t>(workMem.Row2Scan[closestRowScan]);
        
            /* Add the closest column to the list of scanned columns and
             * delete it from the list of columns to scan by shifting all
             * of the items after it over by one.
             */
            workMem.ScannedRows[closestRow]=true;
            
            memmove(workMem.Row2Scan+closestRowScan,workMem.Row2Scan+closestRowScan+1,(numRow2Scan-closestRowScan)*sizeof(ptrdiff_t));
            numRow2Scan--;//One fewer row to scan.           
            
            delta=workMem.shortestPathCost[closestRow];
            
            //If we have reached an unassigned row.
            if(problemSol->col4row[closestRow]==-1) {
                sink=static_cast<ptrdiff_t>(closestRow);
            } else{
                curCol=static_cast<size_t>(problemSol->col4row[closestRow]);
            }            
        } while(sink==-1);
        
/* Next, update the dual variables.*/
        
        updateDualAndAugment(problemSol,workMem,curUnassignedCol, numColsScanned,numRow,sink,delta);
    }
    
    //Determine the gain to return
    calcGain(problemSol,workMem,numRow,numCol4Gain);    
    problemSol->forbiddenActiveRows[problemSol->row4col[0]]=true;
    return 0;
}

MurtyHyp *shortestPathUpdateCPP(const MurtyHyp *parentHyp, ScratchSpace &workMem,const size_t curUnassignedCol, size_t numRow2Scan, const size_t numVarCol, const size_t numDim) {
/*SHORTESTPATHUPDATECPP
 *
 * This is a realization of the update inheriting dual variables as
 * described in 
 * M. L. Miller, H. S. Stone, and I. J. Cox, "Optimizing Murty's ranked
 * assignment method," IEEE Transactions on Aerospace and Electronic
 * Systems, vol. 33, no. 3, pp. 851-862, Jul. 1997.
 * to efficiently implement Murty's algorithm. This function is only used
 * in the k-best 2D assignment algorithm, not in the regular 2D assignment
 * algorithm. The return value is the new solution. The dual variables from
 * the parent hypothesis parentHyp are inherited and modified. This
 * function is called by the split function. If the problem is infeasible,
 * then the gain in problemSol is set to -1.
 *
 */
    
    /* Col2Scan, numCol2Scan, and and forbiddenActiveCols in the scratch
     * space are assumed to have already been set by the parent hypothesis.
     **/
    size_t curRow,curCol,numColsScanned;
    ptrdiff_t sink;
    double delta;
    MurtyHyp *problemSol;

    problemSol= new MurtyHyp(numDim,numDim);

    //Copy the appropriate things that are to be inherited.
    problemSol->activeCol=curUnassignedCol;
    
    memcpy(problemSol->row4col,parentHyp->row4col,numDim*sizeof(ptrdiff_t));
    memcpy(problemSol->col4row,parentHyp->col4row,numDim*sizeof(ptrdiff_t));
    memcpy(problemSol->u,parentHyp->u,numDim*sizeof(double));
    memcpy(problemSol->v,parentHyp->v,numDim*sizeof(double));
    memcpy(problemSol->forbiddenActiveRows,workMem.forbiddenActiveRows,numDim*sizeof(bool));

    //Remove the association of the current row/ column.
    problemSol->col4row[problemSol->row4col[curUnassignedCol]]=-1;
    problemSol->row4col[curUnassignedCol]=-1;

/* Find the shortest augmenting path starting at
 * curUnassignedRow. It is the only unassigned row.*/

    /* Mark everything as not yet scanned. A 1 will be placed in each
     * column entry as it is scanned.*/
    numColsScanned=0;
    fill_n(workMem.ScannedRows,numDim,false);
    /* Initially, the cost of the shortest path to each column is not
     * known and will be made infinite.*/
    fill_n(workMem.shortestPathCost,numDim,numeric_limits<double>::infinity());

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
        /*Mark the current row as having been visited.*/
        workMem.ScannedColIdx[numColsScanned++]=curCol;

        /*Scan all of the columns that have not already been scanned.*/
        minVal=numeric_limits<double>::infinity();
        for(curRowScan=0;curRowScan<numRow2Scan;curRowScan++) {
            curRow=static_cast<size_t>(workMem.Row2Scan[curRowScan]);
                
            if(curCol!=curUnassignedCol||workMem.forbiddenActiveRows[curRow]==false) {
                double reducedCost;
                
                reducedCost=delta+workMem.C[curRow+curCol*numDim]-problemSol->u[curCol]-problemSol->v[curRow];
                if(reducedCost<workMem.shortestPathCost[curRow]){
                    workMem.pred[curRow]=curCol;
                    workMem.shortestPathCost[curRow]=reducedCost;
                }

                //Find the minimum unassigned column that was scanned.
                if(workMem.shortestPathCost[curRow]<minVal){
                    minVal=workMem.shortestPathCost[curRow];
                    closestRowScan=curRowScan;
                }
            }
        }
        
        if(minVal==numeric_limits<double>::infinity()) {
           /* If the minimum cost column is not finite, then the
            * problem is not feasible.*/
            
            problemSol->gain=-1;            
            return problemSol;
        }
        /* Change the index from the relative column index to the
         * absolute column index.*/
        closestRow=static_cast<size_t>(workMem.Row2Scan[closestRowScan]);

        /* Add the closest column to the list of scanned columns and
         * delete it from the list of columns to scan by shifting all
         * of the items after it over by one.
         */
        workMem.ScannedRows[closestRow]=true;

        memmove(workMem.Row2Scan+closestRowScan,workMem.Row2Scan+closestRowScan+1,(numRow2Scan-closestRowScan)*sizeof(ptrdiff_t));
        numRow2Scan--;//One fewer column to scan.

        delta=workMem.shortestPathCost[closestRow];

        //If we have reached an unassigned row.
        if(problemSol->col4row[closestRow]==-1) {
            sink=static_cast<ptrdiff_t>(closestRow);
        } else{
            curCol=static_cast<size_t>(problemSol->col4row[closestRow]);
        }            
    } while(sink==-1);
    
    updateDualAndAugment(problemSol,workMem,curUnassignedCol, numColsScanned,numDim,sink,delta);

    //Determine the gain to return
    calcGain(problemSol,workMem,numDim,numVarCol);
    problemSol->forbiddenActiveRows[problemSol->row4col[curUnassignedCol]]=true;
    return problemSol;
}

void split(MurtyHyp *parentHyp,priority_queue<pMurtyHyp> &HypQueue, ScratchSpace &workMem, const size_t numVarCol,const size_t numDim) {
/*SPLIT
 *
 * This is a realization of a function to split hypothesis for the k-best
 * hypothesis algorithm by adding more constraints as described in
 * K. G. Murty, "An algorithm for ranking all the assignments in order of
 * increasing cost," Operations Research, vol. 16, no. 3, pp. 682-687,
 * May-Jun. 1968.
 * This function then calls shortestPathUpdateCPP to get the best 2D
 * assignments for the hypotheses with varying numbers of constraints.
 *
 */    
    size_t curCol,numRow2Scan,activeCol;
    char *row2RemovePtr;
    MurtyHyp *newHyp;

    activeCol=parentHyp->activeCol;
    
    /* This function splits this hypotheses into a number of
     * sub-hypotheses and solves the association problem for each of
     * them, putting the results into a queue.
     **/
    /* All of the columns that are initially available are those that
     * have not already been assigned in a superproblem.*/
    numRow2Scan=0;
    for(curCol=activeCol;curCol<numDim;curCol++) {
        workMem.Row2ScanParent[numRow2Scan]=parentHyp->row4col[curCol];
        numRow2Scan++;
    }

    //Sort the indices of the columns, so that they can be quickly removed for each of the subproblems.
    qsort(workMem.Row2ScanParent,numRow2Scan,sizeof(ptrdiff_t),compare);
    
    copy(workMem.Row2ScanParent,workMem.Row2ScanParent+numRow2Scan,workMem.Row2Scan);
    //For the first split, all of the extra constraints imposed upon the active row are present.
    copy(parentHyp->forbiddenActiveRows,parentHyp->forbiddenActiveRows+numDim,workMem.forbiddenActiveRows);
    
    /*This function does not modify numCol2ScanParent */    
    newHyp=shortestPathUpdateCPP(parentHyp,workMem,activeCol,numRow2Scan,numVarCol,numDim);
    
    /*If it is not a missed detection, add it to the queue*/
    if(newHyp->gain==-1) {
        delete newHyp;
    }
    else {
        HypQueue.push(pMurtyHyp(newHyp));
    }    

    /* Remove the current assignment from the list of columns that can 
     * be scanned. First, we must find the index of the column that is 
     * to be removed.*/
    row2RemovePtr=reinterpret_cast<char*>(bsearch(&parentHyp->row4col[activeCol],workMem.Row2ScanParent,numRow2Scan,sizeof(ptrdiff_t),compare));
    memmove(row2RemovePtr,row2RemovePtr+sizeof(ptrdiff_t),static_cast<size_t>(reinterpret_cast<char*>(workMem.Row2ScanParent+numRow2Scan)-row2RemovePtr));
    numRow2Scan--;//One fewer column to scan.

    fill_n(workMem.forbiddenActiveRows,numDim,false);
    for(curCol=activeCol+1;curCol<numVarCol;curCol++){
        copy(workMem.Row2ScanParent,workMem.Row2ScanParent+numRow2Scan,workMem.Row2Scan);
        /* Since these are rows after the active row, only the current
         * assignment is invalid.
         */
        workMem.forbiddenActiveRows[parentHyp->row4col[curCol]]=true;
        
        newHyp=shortestPathUpdateCPP(parentHyp,workMem,curCol,numRow2Scan,numVarCol,numDim);
        
        /*If it is not a missed detection, add it to the queue*/
        if(newHyp->gain==-1) {delete newHyp;}
        else {HypQueue.push(pMurtyHyp(newHyp));}

        /*Remove the current assignment from the list of columns to be scanned.*/
        row2RemovePtr=reinterpret_cast<char*>(bsearch(parentHyp->row4col+curCol,workMem.Row2ScanParent,numRow2Scan,sizeof(ptrdiff_t),compare));
        memmove(row2RemovePtr,row2RemovePtr+sizeof(ptrdiff_t),static_cast<size_t>(reinterpret_cast<char*>(workMem.Row2ScanParent+numRow2Scan)-row2RemovePtr));
        numRow2Scan--;//One fewer column to scan.

        //Unforbid the current column.
        workMem.forbiddenActiveRows[parentHyp->row4col[curCol]]=false;
    }
}

double makeCostMatrixSafe(double *CMod,const double *COrig,const size_t numEl, const bool maximize) {
/*MAKECOSTMATRIXSAFE
 * This function adjusts the cost matrix so that a shortest augmenting path
 * algorithm can be used. Specifically, it offsets the elements so that the
 * optimal solution can be obtained as a minimization problem of a cost
 * matrix with all positive elements. The results are placed in CMod and
 * the amount by which the elements of the matrix matrix were shifted
 * (after possibly being negatived) is returned.
 *
 *INPUTS: CMod  A pointer to a cost matrix that can hold the result.
 *        COrig A pointer to the original cost matrix that is to be
 *              adjusted. CMod and COrig can be the same matrix.
 *        numEl The number of elements in the matrix COrig.
 *     maximize True is the adjustment is for a cost matrix that will be
 *              used in a maximization problem, false otherwise.
 *
 **/
double CDelta;
size_t i;

    if(maximize==false) {
        CDelta = *min_element(COrig, COrig + numEl);

        for(i=0;i<numEl;i++) {
            CMod[i]=COrig[i]-CDelta;
        }
    } else {
        CDelta = *max_element(COrig, COrig + numEl);

        for(i=0;i<numEl;i++) {
            CMod[i]=-COrig[i]+CDelta;
        }
    }
    
    return CDelta;
}

size_t kBest2D(const size_t k,const size_t numRow,const size_t numCol,const bool maximize,const double *C, ScratchSpace &workMem,ptrdiff_t *col4rowBest,ptrdiff_t *row4colBest,double *gainBest) {
    size_t curSweep;
    MurtyHyp *curHyp=new MurtyHyp(numRow,numRow);
    priority_queue<pMurtyHyp> HypQueue;
    double CDelta;
    
    /* The cost matrix must have all non-negative elements for the
     * assignment algorithm to work. This forces all of the elements to be
     * positive. The delta is added back in when computing the gain in the
     * end. workMem.C must hold enough space for a numRowXnumRow matrix.
     * since numrow>=numCol, the extra columns will be set to zero.*/
    CDelta=makeCostMatrixSafe(workMem.C,C,numRow*numCol,maximize);
    CDelta=CDelta*numCol;
    //The rest of the entries in the working memory matrix are zero. 
    fill_n(workMem.C+numRow*numCol,numRow*(numRow-numCol),0);

    //First, solve the full problem for the best hypothesis.
    if(shortestPathCPP(curHyp,workMem,numRow,numRow,numCol)) {
    /*If the problem is infeasible, then identify it as such and return.
     */
        delete curHyp;
        return 0;
    }

    copy(curHyp->col4row,curHyp->col4row+numRow,col4rowBest);        
    copy(curHyp->row4col,curHyp->row4col+numCol,row4colBest);
    gainBest[0]=curHyp->gain;
    /* Adjust for shifting that was done to make everything positive.*/
    if(maximize==false) {
        gainBest[0]=gainBest[0]+CDelta;
    } else {
        gainBest[0]=-gainBest[0]+CDelta;
    }
    
    HypQueue.push(pMurtyHyp(curHyp));
    //Enter the loop to generate all of the other hypotheses from this one.
    for(curSweep=1;curSweep<k;curSweep++) {
        curHyp=HypQueue.top().ptr;
        HypQueue.pop();
        
        //The current hypothesis must now be split.
        split(curHyp,HypQueue, workMem, numCol,numRow);
        
        //The hypothesis queue now contains all of the ordered, new
        //hypotheses and the old current hypothesis is no longer needed.
        delete curHyp;

        //Save the result in the output. If the queue is empty, then all
        //possible hypotheses have been generated and we can return.
        if(HypQueue.empty()==false) {
            curHyp=HypQueue.top().ptr;
            copy(curHyp->col4row,curHyp->col4row+numRow,col4rowBest+curSweep*numRow);
            copy(curHyp->row4col,curHyp->row4col+numCol,row4colBest+curSweep*numCol);
            gainBest[curSweep]=curHyp->gain;
        /* Adjust for shifting that was done to make everything positive.*/
            if(maximize==false) {
                gainBest[curSweep]=gainBest[curSweep]+CDelta;
            } else {
                gainBest[curSweep]=-gainBest[curSweep]+CDelta;
            }
        } else {
            break;
        }
    }
    
    /* We are done and all of the old junk on the queue can be
     * deleted.*/
    while(HypQueue.empty()==0){
        delete HypQueue.top().ptr;
        HypQueue.pop();
    }

    return curSweep;
}

int assign2D(const size_t numRow,const size_t numCol,const bool maximize,const double *C, ScratchSpace &workMem,MurtyHyp *problemSol) {
/*ASSIGN2D Perform 2D assignment after making adjusting the cost matrix to
 *         be safe and transforming the optimization problem into a
 *         minimization problem.
 **/
    double CDelta;
    
    /* The cost matrix must have all non-negative elements for the
     * assignment algorithm to work. This forces all of the elements to be
     * positive. The delta is added back in when computing the gain in the
     * end. workMem.C must hold enough space for a numRowXnumCol matrix.*/
    CDelta=makeCostMatrixSafe(workMem.C,C,numRow*numCol,maximize);
    CDelta=CDelta*numCol;
    
    if(shortestPathCPP(problemSol,workMem,numRow,numCol,numCol)) {
    /*If the problem is infeasible, then identify it as such and return.
     */
        return 0;
    }
    
    if(maximize==false) {
        problemSol->gain=problemSol->gain+CDelta;
    } else {
        problemSol->gain=-problemSol->gain+CDelta;
    }
    
    return 1;
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
