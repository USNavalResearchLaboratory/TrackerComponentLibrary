/**ASSIGN2DFULL A C (for Matlab) implementation of a function to perform
*      the type of 2D assignment that often arises as a subproblem when
*      solving an S-dimensional assignment problem for target tracking,
*      where one row and one column hold missed detection costs.
*      See the comments to the Matlab implementation for more details.
*
*INPUTS: C A (numRow)X(numCol) cost matrix where the first row is the
*          cost of not assigning each column (from number 2) and the first
*          column is the cost of not assigning each row (from number 2) and
*          C(1,1) is an assignment cost that is unconstrained and only
*          assigned if it improves the cost function (given above). The
*          matrix cannot contain any NaNs and if maximizing, cannot
*          contain any Inf values and if minimizing cannot contain any -Inf
*          values.
* maximize If true, the minimization problem is transformed into a
*          maximization problem. The default if this parameter is omitted
*          or an empty matrix is passed is false.
* algorithm A value selecting which algorithm will be used. The basis for
*          the algorithms are discussed in more detail below. Possible
*          values are:
*          0 (The default if omitted or an empty matrix is passed) Use an
*            algorithm that effectively solves an optimization problem on a
*            (numRow-1)+(numCol-1) augmented matrix, though the matrix is
*            not explicitly formed.
*          1 Use an algorithm that removes the first column and modifies
*            the cost matrix so that the optimization problem can be
*            solved using the assign2DMissedDetect function, which
*            effectively solves a problem on a
*            (numRow+numCol-1)X(numCol-1) augmented matrix. Unlike
*            algorithm 0, this algorithm requires that there not be
*            non-finite terms in both the first row and the first column of
*            C. Otherwise, the problem will be labeled as infeasible.
*
*OUTPUTS: tuples A 2XnumTuples set of assignment values. The returned data
*               type is double, to be consistent with Matlab's defaults,
*               even though the values are all unsigned integers and an
*               integer data type might normally be chosen. This is ordered
*               [Row Index;Column Index]. the number of tuples ranges from
*               max(numRow-1,numCol-1), for the case where as many things
*               as possible are not assigned to the uncontrained values and
*               C(1,1) is not assigned, to numRow+numCol-1 for the case
*               where all rows and columns are assigned to the
*               unconstrained hypothesis and C(1,1) is provided. If the
*               problem is infeasible, this is an empty matrix.
*          gain This is the value of the cost. This is the sum of the
*               values in C corresponding to the tuples. If
*               the problem is infeasible, this is -1.
*           u,v These are the dual variables for the constrained columns
*               and for constrained rows of C. These will satisfy the
*               complementary slackness condition of the problem as
*               transformed. See the example below.
*
*The algorithm transforms the problem into a minimization problem and
*solves it using a modified Jonker-Volgenant algorithm as described in the
*comments to the Matlab implementation.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm)
*
*January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//The header for the assignment algorithms in C.
#include "assignAlgs2D.h"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t numRow, numCol;
    bool isInfeasible, maximize=false;
    mxArray *gainMATLAB, *uMATLAB, *vMATLAB;
    double *CIn, *gainPtr, *u, *v;
    int algorithm=0;
    void *tempBuffer;
    ptrdiff_t *tuples;
    size_t numTuples;

    if(nrhs>3||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }

    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        maximize=getBoolFromMatlab(prhs[1]);
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        algorithm=getIntFromMatlab(prhs[2]);
    }
    
    /*Verify the validity of the assignment matrix.*/
    checkRealDoubleArray(prhs[0]);
    
    //If an empty matrix is passed, return the appropriate values.
    if(mxIsEmpty(prhs[0])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//tuples
        if(nlhs>1) {
            gainMATLAB=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            gainPtr=mxGetDoubles(gainMATLAB);
            *gainPtr=0;
            plhs[1]=gainMATLAB;//Gain
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleMatrix(0,0,mxREAL);//u
                if(nlhs>3) {
                    plhs[3]=mxCreateDoubleMatrix(0,0,mxREAL);//v
                }
            }
        }
        return;
    }

    /* Get the dimensions of the input data. It is assumed that the matrix
     * is not so large in either dimension as to cause an overflow when
     * using a SIGNED integer data type.*/
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    CIn=mxGetDoubles(prhs[0]);

    switch(algorithm){
        case 0:
        {
            double *C;
            const size_t augmentedSize=(numRow-1)+(numCol-1);
            
            {
                const size_t numBytes=numRow*numCol*sizeof(double);
                C=(double*)mxMalloc(numBytes);
                //Duplicate the input array.
                memcpy(C,CIn,numBytes);
            }

            //Allocate space for temporary variables required by the
            //assign2DFullC function plus the space for the tuples, which
            //require at most 2*maxNumTuples size_t elements.
            {
                const size_t maxNumTuples=(augmentedSize+1);
                tempBuffer=mxMalloc(assign2DFullCBufferSize(numRow,numCol)+2*maxNumTuples*sizeof(ptrdiff_t));
            }

            //Allocate the return values.
            gainMATLAB=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            uMATLAB=mxCreateNumericMatrix(augmentedSize,1,mxDOUBLE_CLASS,mxREAL);
            vMATLAB=mxCreateNumericMatrix(augmentedSize,1,mxDOUBLE_CLASS,mxREAL);

            gainPtr=mxGetDoubles(gainMATLAB);
            u=mxGetDoubles(uMATLAB);
            v=mxGetDoubles(vMATLAB);

            {
                //The first part of the buffer holds the tuples, the second
                //part is the buffer to pass to assign2DFullC.
                void* bufferStart=(void*)((char*)(tempBuffer)+2*(augmentedSize+1)*sizeof(size_t));
                tuples=(ptrdiff_t*)(tempBuffer);

                isInfeasible=assign2DFullC(maximize,C,gainPtr,tuples,&numTuples,bufferStart,u,v,numRow,numCol);
            }
            mxFree(C);
            break;
        }
        case 1:
        {                
            //Allocate space for temporary variables required by the
            //assign2DFullCAlt function plus the space for the tuples,
            //which require a maximum of 2*maxNumTuples size_t elements.
            void* bufferStart;
            //The number of rows needed in the augmented matrix inside the
            //assign2DMissedDetectC function, given that in
            //assign2DFullCAlt we remove a column before calling the
            //assign2DMissedDetectC function.
            const size_t numRowAug=(numRow-1)+(numCol-1);
            const size_t maxNumTuples=numRowAug+1;
            
            tempBuffer=mxMalloc(assign2DFullCAltBufferSize(numRow,+numCol)+2*maxNumTuples*sizeof(ptrdiff_t));

            //Allocate the return values.
            gainMATLAB=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);

            if((numRow-1)==0&&(numCol-1)==0) {
                vMATLAB=mxCreateNumericMatrix(0,1,mxDOUBLE_CLASS,mxREAL);
                uMATLAB=mxCreateNumericMatrix(0,1,mxDOUBLE_CLASS,mxREAL);
            } else {
                size_t maxMemVal=(numRow-1)+(numCol-1);

                if((numCol-1)>maxMemVal) {
                    maxMemVal=numCol-1;
                }

                uMATLAB=mxCreateNumericMatrix(maxMemVal,1,mxDOUBLE_CLASS,mxREAL);
                vMATLAB=mxCreateNumericMatrix(maxMemVal,1,mxDOUBLE_CLASS,mxREAL);
            }

            gainPtr=mxGetDoubles(gainMATLAB);
            u=mxGetDoubles(uMATLAB);
            v=mxGetDoubles(vMATLAB);

            //The first part of the buffer holds the tuples, the second
            //part is the buffer to pass to assign2DFullCAlt.
            tuples=(ptrdiff_t*)(tempBuffer);
            bufferStart=(void*)((size_t*)(tempBuffer)+2*maxNumTuples);
            isInfeasible=assign2DFullCAlt(maximize,CIn,gainPtr,tuples,&numTuples,bufferStart,u,v,numRow,numCol);
            
            break;
        }
        default:
            mexErrMsgTxt("Invalid algorithm specified.");
            return;
    }

    if(isInfeasible) {
        mxFree(tempBuffer);
        mxDestroyArray(gainMATLAB);
        mxDestroyArray(uMATLAB);
        mxDestroyArray(vMATLAB);
        plhs[0]=mxCreateDoubleMatrix(2,0,mxREAL);//tuples
        if(nlhs>1) {
            gainMATLAB=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            gainPtr=mxGetDoubles(gainMATLAB);
            *gainPtr=-1;
            plhs[1]=gainMATLAB;//Gain
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleMatrix(0,0,mxREAL);//u
                if(nlhs>3) {
                    plhs[3]=mxCreateDoubleMatrix(0,0,mxREAL);//v
                }
            }
        }

        return;
    } else {
        mxArray *tuplesMATLAB;
        size_t i;
        
        //The tuples in C++ are indexed from zero. Add 1 to all of them to
        //get the Matlab indexation.
        for(i=0;i<2*numTuples;i++) {
            tuples[i]++;
        }
        
        //Extract the tuples from the temporary buffer and convert them to
        //doubles for return to Matlab.
        tuplesMATLAB=ptrDiffTMat2MatlabDoubles(tuples,2,numTuples);
        //Delete the temporary buffer, which held the tuples as a matrix of
        //size_t elements.
        mxFree(tempBuffer);

        plhs[0]=tuplesMATLAB;
        if(nlhs>1) {
            plhs[1]=gainMATLAB;
            if(nlhs>2) {
                //Only keep the dual variables for real columns and rows.
                mxSetM(uMATLAB,numCol-1);
                plhs[2]=uMATLAB;
                if(nlhs>3) {
                    mxSetM(vMATLAB,numRow-1);
                    plhs[3]=vMATLAB;
                } else {
                    mxDestroyArray(vMATLAB);
                }
            } else {
                mxDestroyArray(uMATLAB);
                mxDestroyArray(vMATLAB);
            }
        } else {
            mxDestroyArray(gainMATLAB);
            mxDestroyArray(uMATLAB);
            mxDestroyArray(vMATLAB);
        }
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
