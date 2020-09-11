/**KBEST2DASSIGNMENT Find the k lowest (or highest) cost 2D assignments for
 *                   the two-dimensional assignment problem with a
 *                   rectangular cost matrix C.
 *
 *INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
 *          where the largest finite element minus the smallest element is
 *          a finite quantity (does not overflow) when performing
 *          minimization and where the smallest finite element minus the
 *          largest element is finite when performing maximization.
 *          Forbidden assignments can be given costs of +Inf for
 *          minimization and -Inf for maximization.
 *        k The number >=1 of hypotheses to generate. If k is less than the
 *          total number of unique hypotheses, then all possible hypotheses
 *          will be returned.
 * maximize If true, the minimization problem is transformed into a
 *          maximization problem. The default if this parameter is omitted
 *          is false.
 *
 *OUTPUTS: col4rowBest A numRowXk vector where the entry in each element
 *                      is an assignment of the element in that row to a
 *                      column. 0 entries signify unassigned rows.
 *          row4colbest A numColXk vector where the entry in each element
 *                      is an assignment of the element in that column to a
 *                      row. 0 entries signify unassigned columns.
 *          gainBest    A kX1 vector containing the sum of the values of
 *                      the assigned elements in C for all of the
 *                      hypotheses.
 *DEPENDENCIES: mex.h
 *              <algorithm>
 *              MexValidation.h
 *              ShortestPathCPP.hpp
 *              ShortestPathCPP.cpp
 *
 *This is an implementation of Murty's method, which is described in [1].
 *The algorithm relies on the existence of a 2D assignment algorithm. The
 *2D assignment algorithm of [2] is used. Additionally, the dual variable
 *inheritance methods described in [3] is used to reduce the computational
 *complexity of the technique.
 *
 *Murty's algorithm runs 2D assignment algorithms a number of times with an
 *increasing number of constraints. Much of the assignment code is in the
 *handle subclass MurtyData. Instances of MurtyData are stored in an
 *ordered list implemented using the BinaryHeap class.
 *
 * The algorithm can be compiled for use in Matlab  using the command
 * mex('-v','-largeArrayDims','kBest2DAssign.cpp','ShortestPathCPP.cpp')
 *
 * The algorithm relies on the files
 * MexValidation.h
 * ShortestPathCPP.hpp
 * ShortestPath.cpp
 * This file is primarily a MATLAB wrapper for the algorithm that
 * validates the input and formats it for use in the function kBest2DCPP,
 * which is implemented with the single-hypothesis 2D assignment algorithm
 * in ShortestPath.cpp.
 *
 * The algorithm is run in Matlab using the command format
 * [col4row,row4col,gain]=kBest2DAssign(C,k,maximize)
 *
 *REFERENCES:
 *[1] K. G. Murty, "An algorithm for ranking all the assignments in order 
 *    of increasing cost," Operations Research, vol. 16, no. 3, pp. 682-
 *    687, May-Jun. 1968.
 *[2] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
 *    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no.
 *    4, pp. 1679-1696, Aug. 2016.
 *[3] M. L. Miller, H. S. Stone, and J. Cox, Ingemar, "Optimizing Murty's
 *    ranked assignment method," IEEE Transactions on Aerospace and
 *    Electronic Systems, vol. 33, no. 3, pp. 851-862, Jul. 1997.
 *
 *November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
#include <algorithm>
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "ShortestPathCPP.hpp"

using namespace std;

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t numRow,numCol,k,numFound;
    mxArray *CMat, *col4rowMATLAB, *row4colMATLAB, *gainMATLAB;//These will hold the values to be returned.
    ptrdiff_t *col4rowBest, *row4colBest;
    double *gainBest;
    ScratchSpace workMem;//Scratch space needed for the assignment algorithm.
    bool didFlip=false;
    bool maximize=false;
    
    if(nrhs<2){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    k=getSizeTFromMatlab(prhs[1]);
    if(k<=0){
        mexErrMsgTxt("Invalid number of hypotheses requested.");
        return;
    }

    if(nrhs==3) {
        maximize=getBoolFromMatlab(prhs[2]);
    }
    
    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Too many outputs.");
        return;
    }
    
    /*Verify the validity of the assignment matrix.*/
    checkRealDoubleArray(prhs[0]);
    /* Get the dimensions of the input data and the pointer to the matrix.
     * It is assumed that the matrix is not so large in M or N as to cause
     * an overflow when using a SIGNED integer data type.*/
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);

    /* Transpose the matrix, if necessary, so that the number of rows is
     * >= the number of columns.*/
    if(numRow>=numCol) {
        CMat=mxDuplicateArray(prhs[0]);
    } else {        
        swap(numRow,numCol);
        
        //This is freed using mxDestroyArray
        CMat=mxCreateDoubleMatrix(numCol,numRow,mxREAL);
        mexCallMATLAB(1, &CMat, 1,  (mxArray **)&prhs[0], "transpose");
        didFlip=true;
    }
    
    //Allocate scratch space.
    workMem.init(numRow,numRow);

    //Allocate space for the return variables
    col4rowMATLAB=allocSignedSizeMatInMatlab(numRow,k);
    row4colMATLAB=allocSignedSizeMatInMatlab(numCol,k);

    gainMATLAB = mxCreateNumericMatrix(k,1,mxDOUBLE_CLASS,mxREAL);
    if(sizeof(ptrdiff_t)==4) {//32 bit
    	col4rowBest=(ptrdiff_t*)mxGetInt32s(col4rowMATLAB);
        row4colBest=(ptrdiff_t*)mxGetInt32s(row4colMATLAB);
    } else {//64 bit
        col4rowBest=(ptrdiff_t*)mxGetInt64s(col4rowMATLAB);
        row4colBest=(ptrdiff_t*)mxGetInt64s(row4colMATLAB);
    }
    gainBest=mxGetDoubles(gainMATLAB);

    /*The assignment algorithm returns a nonzero value if no valid
     * solutions exist.*/
    numFound=kBest2D(k,numRow,numCol,maximize, mxGetDoubles(CMat), workMem, col4rowBest,row4colBest,gainBest);
    mxDestroyArray(CMat);
    
    if(numFound==0){
        mwSize dims[2] = {0,0};
        const mwSize numDims=2;
        mxSetDimensions(col4rowMATLAB, dims, numDims);
        mxSetDimensions(row4colMATLAB, dims, numDims);
        
        dims[0]=1;
        dims[1]=1;
        mxSetDimensions(gainMATLAB, dims, numDims);
        gainBest[0]=-1;
        
        //Set the outputs
        mxDestroyArray(col4rowMATLAB);
        mxDestroyArray(row4colMATLAB);
        mxDestroyArray(gainMATLAB);
        
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        if(nlhs>1){
            plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);
            if(nlhs>2){
                plhs[2]=mxCreateDoubleScalar(-1.0);
            }
        }
        return;
    }

    {
        size_t dims[2];
        const size_t numDims=2;
        size_t i;
        /* Reduce the return matrices to the number of hypotheses that
         * were actually found, increment them to be regular Matlab
         * indices, zero all the row4colBest that were just padding. */

        switch(nlhs) {
        case 3:
            dims[0]=numFound;
            dims[1]=1;
        case 2:
            mxSetDimensions(gainMATLAB, dims, numDims); 
            
            dims[0]=numRow;
            dims[1]=numFound;
            /*Convert the indices to Matlab indices*/
            for(i=0;i<dims[0]*dims[1];i++) {
                if(col4rowBest[i]<(ptrdiff_t)numCol){
                    col4rowBest[i]+=1;
                } else {
                  col4rowBest[i]=0;  
                }
            }
            mxSetDimensions(col4rowMATLAB, dims, numDims);
        default:
            dims[0]=numCol;
            dims[1]=numFound;
            /*Convert the indices to Matlab indices*/
            for_each(row4colBest, row4colBest+dims[0]*dims[1], increment<ptrdiff_t>);

            mxSetDimensions(row4colMATLAB, dims, numDims);
        }
    }

    /* If a transposed array was used */
    if(didFlip==true) {        
        swap(numRow,numCol);
        swap(row4colMATLAB,col4rowMATLAB);
    }

    //Let Matlab know that these are the return variables.
    switch(nlhs) {
        case 3:
            plhs[2]=gainMATLAB;
        case 2:
            plhs[1]=row4colMATLAB;
        default:
            plhs[0]=col4rowMATLAB;
    }
    return;
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
