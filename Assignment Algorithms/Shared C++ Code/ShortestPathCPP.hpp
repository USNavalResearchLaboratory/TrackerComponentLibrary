/*SHORTESTPATHCPP A header file for functions and data structures 
 *                implementing the shortest augmenting path algorithms to
 *                get the best or the k-best hypotheses for a 2D assignment
 *                problem. 
 *
 *This file needs to be compiled with the file ShortestPathCPP.cpp
 *
 *November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef SPALGS
#define SPALGS
#include <stddef.h>

/**The MurtyHyp class is used to hold a solution to the 2D assignment
 * algorithm as well as additional information that is useful when
 * implementing Murty's k-Best 2D assignment algorithm.
 **/
class MurtyHyp {
private:
    char *buffer;
public:
    ptrdiff_t *col4row;
    ptrdiff_t *row4col;
    double gain;
    double *u;
    double *v;
    /*activeCol and forbiddenActiveRows are used in the k-best 2D 
     *assignment algorithm, but not in the regular 2D assignment
     *algorithm.*/
    size_t activeCol;
    bool *forbiddenActiveRows;
    
    MurtyHyp(){
        buffer=NULL;
    }
    
    MurtyHyp(const size_t numRow, const size_t numCol) {
       char *basePtr;
    /*To minimize the number of calls to memory allocation and deallocation
     * routines, a big chunk of memory is allocated at once and pointers
     * to parts of it for the different variables are saved.*/
       buffer=new char[sizeof(ptrdiff_t)*numRow+sizeof(ptrdiff_t)*numCol+sizeof(double)*numCol+sizeof(double)*numRow+sizeof(bool)*numRow];
       basePtr=buffer;
       col4row=reinterpret_cast<ptrdiff_t*>(basePtr);
       basePtr+=numRow*sizeof(ptrdiff_t);
       row4col=reinterpret_cast<ptrdiff_t*>(basePtr);
       basePtr+=numCol*sizeof(ptrdiff_t);
       u=reinterpret_cast<double*>(basePtr);
       basePtr+=sizeof(double)*numCol;
       v=reinterpret_cast<double*>(basePtr);
       basePtr+=sizeof(double)*numRow;
       forbiddenActiveRows=reinterpret_cast<bool*>(basePtr);
    }
    
    ~MurtyHyp(){
        if(buffer!=NULL){
            delete[] buffer;
        }
    }
}; 

/* The ScratchSpace class holds the parameters and scratch space for the 2D
 * assignment algorithm so as to reduce the number of memory allocation and
 * deallocation operations that are necessary when the 2D assignment 
 * algorithm must be called repeatedly. Some of the entries are only used
 * in the k-best implementation of the algorithm.
 */
class ScratchSpace {
public:
    char *buffer;
    double *C;
    size_t *ScannedColIdx;
    bool *ScannedRows;
    size_t *pred;
    double *shortestPathCost;
    ptrdiff_t *Row2ScanParent;
    ptrdiff_t *Row2Scan;
    bool* forbiddenActiveRows;
    
    //The constructor
    ScratchSpace(){
        buffer=NULL;
    }
    
    ScratchSpace(const size_t numRow,const size_t numCol){
        this->init(numRow,numCol);
    }
    
    void init(const size_t numRow,const size_t numCol){
        char *basePtr;
    /*To minimize the number of calls to memory allocation and deallocation
     * routines, a big chunk of memory is allocated at once and pointers
     * to parts of it for the different variables are saved.*/
        buffer=new char[numCol*sizeof(size_t)+2*numRow*sizeof(ptrdiff_t)+numRow*sizeof(size_t)+numRow*(1+numCol)*sizeof(double)+2*numRow*sizeof(bool)];
        basePtr=buffer;
        ScannedColIdx=reinterpret_cast<size_t*>(basePtr);
        basePtr+=sizeof(size_t)*numCol;
        Row2ScanParent=reinterpret_cast<ptrdiff_t*>(basePtr);
        basePtr+=sizeof(ptrdiff_t)*numRow;
        Row2Scan=reinterpret_cast<ptrdiff_t*>(basePtr);
        basePtr+=sizeof(ptrdiff_t)*numRow;
        pred=reinterpret_cast<size_t*>(basePtr);
        basePtr+=sizeof(size_t)*numRow;
        shortestPathCost=reinterpret_cast<double*>(basePtr);
        basePtr+=sizeof(double)*numRow;
        ScannedRows=reinterpret_cast<bool*>(basePtr);
        basePtr+=sizeof(bool)*numRow;
        forbiddenActiveRows=reinterpret_cast<bool*>(basePtr);
        basePtr+=sizeof(bool)*numRow;

        C=reinterpret_cast<double*>(basePtr);
    }
    
    ~ScratchSpace(){
        if(buffer!=NULL) {
            delete[] buffer;
        }
    }
};

int assign2D(const size_t numRow,
             const size_t numCol,
             const bool maximize,
             const double *C,
             ScratchSpace &workMem,
             MurtyHyp *problemSol);
/*ASSIGN2D Perform 2D assignment using a shortest augmenting path algorithm
 *         that scans by row. This function transforms the maximization
 *         problems into equivalent minimization problems having cost
 *         matrices of all positive elements. 
 *
 *INPUTS:numRow The number of rows in the cost matrix.
 *       numCol The number of columns in the cost matrix. Note that
 *              numRow>=numCol.
 *     maximize True if the optimization is a maximization
 *          C   The cost matrix.
 *      workMem An instance of the ScratchSpace class that was initialized
 *              with workMem.init(numRow,numCol);
 *   ProblemSol An instance of MurtyHyp created using
 *              MurtyHyp(numRow,numCol) in which the solution to the
 *              assignment problem is placed.
 *
 *OUTPUTS: The results are placed in MurtyHyp. The return value is 0 if no
 *         optimal solution with finite cost exists. It is one otherwise.
 *         That is, the return value is the number of solutions found.
 *
 * The algorithm of
 * D. F. Crouse, "Advances in displaying uncertain estimates of multiple
 * targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, and
 * Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
 * is used.
 *
 **/

int shortestPathCPP(MurtyHyp *problemSol,
                   ScratchSpace &workMem,
                   const size_t numRow,
                   const size_t numCol,
                   const size_t numCol4Gain);
/*SHORTESTPATHCPP
 *
 * The shortest augmenting path algorithm for 2D assignment. This algorithm
 * assumes that one is performing a minimization, that numRow>=numCol, and
 * that all of the elements in the cost matrix are non-negative. If these
 * conditions do not hold, then the function assign2D should be used. That
 * function prepares the input and then calls this function.
 * Unlike assign2D, this function has an additional input called
 * numCol4Gain, which is useful when this function is used to start Murty's
 * algorithm for the k-best hypotheses. numCol4Gain is the number of
 * columns used when computing the gain. Normally, this will be the same as
 * numCol, but if the matrix was padded for use in a k-best assignment
 * algorithm, then this will be less than the gain.
 * 
 * The result is placed in workMem. the return value is 1 if the problem is
 * infeasible; otherwise it is zero. If the problem is infeasible, then the
 * gain in problemSol is set to -1.
 *
 **/


size_t kBest2D(const size_t k,
               const size_t numRow,
               const size_t numCol,
               const bool maximize,
               const double *C,
               ScratchSpace &workMem,
               ptrdiff_t *col4rowBest,
               ptrdiff_t *row4colBest,
               double *gainBest);
/*KBEST2D         Finds the k-Best 2D assignments using a shortest
 *                augmenting path algorithm that scans by row.
 *
 *INPUTS:   k   The number of solutions to find. If fewer than k solutions
 *              exist, then the maximum number of solutions will be
 *              returned.
 *       numRow The number of rows in the cost matrix.
 *       numCol The number of columns in the cost matrix. Note that
 *              numRow>=numCol.
 *     maximize True if the optimization is a maximization
 *          C   The cost matrix.
 *      workMem An instance of the ScratchSpace class that was initialized
 *              with workMem.init(numRow,numRow);
 *  col4rowBest An array with space for numRow*k elements to hold the
 *              assignments of rows to columns for each of the hypotheses.
 *  row4colBest An array with space for numCol*k elements to hold the
 *              assignments of rows to columns for each of the hypotheses.
 *     gainBest A length-k array to hold the gain (cost) of each
 *              assignment.
 *
 *OUTPUTS: The results are placed in col4rowBest, row4colBest and gainBest.
 *         The function returns the number of solutions found. That will
 *         always be less than or equal to k.
 *
 * This is an implementation of Murty's method, which is described in 
 * K. G. Murty, "An algorithm for ranking all the assignments in order of
 * increasing cost," Operations Research, vol. 16, no. 3, pp. 682-687,
 * May-Jun. 1968.
 * The algorithm relies on the existence of a 2D assignment algorithm. The
 * 2D assignment algorithm of
 * D. F. Crouse, "Advances in displaying uncertain estimates of multiple
 * targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, and
 * Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
 * is used. Additionally, the dual variable inheritance methods described
 * in
 * M. L. Miller, H. S. Stone, and I. J. Cox, "Optimizing Murty's ranked
 * assignment method," IEEE Transactions on Aerospace and Electronic
 * Systems, vol. 33, no. 3, pp. 851-862, Jul. 1997.
 * is used to reduce the computational complexity of the technique.
 *
 **/

template <class T> void increment(T &x){
    x++;
}
/* increment                A function that increments its input. This is
 *                          useful when converting C indices to Matlab
 *                          indices.
 */

#endif

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
