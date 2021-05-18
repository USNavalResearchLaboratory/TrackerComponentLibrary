/**ASSIGN2DMISSEDDETECT A C-code (for Matlab) implementation of an
*                algorithm to solve the two-dimensional assignment problem
*                common to target tracking, where, in a tracking context,
*                one row (or column) represents missed detection costs.
*                See the comments to the Matlab implementation for more
*                details.
*
*INPUTS: C A (numMeas+1)XnumTar (or if targetByRow=false, a
*          numTarX(numMeas+1) cost matrix) where the first row (or column
*          if targetByRow=false) is the cost of a missed detection for each
*          target and the subsequent rows (columns) are the costs of
*          assigning a measurement to a target. The matrix cannot contain
*          any NaNs and if maximizing, cannot containing any Inf values and
*          if minimizing cannot contain any -Inf values.
* maximize If true, the minimization problem is transformed into a
*          maximization problem. The default if this parameter is omitted
*          or an empty matrix is passed is false.
* missedDetectByRow If this is true, the first row of C is the set of
*          missed detection costs. Otherwise, the first column of C is the
*          set of missed detection costs. The default if this parameter is
*          omitted or an empty matrix is passed is true.
*
*OUTPUTS: tuples A 2XnumTar set of assignment values. If
*               targetsByRow=false (or is omitted), then this is ordered
*               [Measurement Index;Target Index]; Otherwise, the ordering
*               is [Target Index; Measurement Index]. If the problem is
*               infeasible, this is an empty matrix.
*          gain This is the value of the cost. This is the sum of the
*               values in C corresponding to the tuples. If the problem is
*               infeasible, this is -1.
*           u,v These are the dual variables for the constrained columns
*               and for constrained rows of C. These will satisfy the
*               complementary slackness condition of a transformed
*               problem.
*
*The algorithm transforms the problem into a minimization rpoblems and
*solve it using a modified Jonker-Volgenant algorithm as described in the
*comments to the Matlab implementation.
*
*The algorithm is described in more detail in the comments to the Matlab
*implementation.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
*[tuples,gain,u,v]=assign2DMissedDetect(C,maximize,missedDetectByRow)
*
*December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"

//This makes sure that the bool type is defined and provides functions for
//handling Matlab data.
#include "MexValidation.h"

//The header for the assignment algorithms in C.
#include "assignAlgs2D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    bool maximize=false;
    bool missedDetectByRow=true;
    mxArray *gainMATLAB, *uMATLAB, *vMATLAB;
    double *gainPtr;
    ptrdiff_t *tuples;
    size_t numRow, numCol, numRowsTrue;
    double *CIn, *C;
    void *tempBuffer;
    bool isInfeasible;

    if(nrhs<1) {
        mexErrMsgTxt("Not enough inputs.");
    }

    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        maximize=getBoolFromMatlab(prhs[1]);
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        missedDetectByRow=getBoolFromMatlab(prhs[2]);
    }
 
    /*Verify the validity of the assignment matrix.*/
    checkRealDoubleArray(prhs[0]);

    /* Get the dimensions of the input data. It is assumed that the matrix
     * is not so large in either dimension as to cause an overflow when
     * using a SIGNED integer data type.*/
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    
    /*Transpose the matrix, if necessary, so that we actually do the
     *missed detection analysis by column. The matrix must be duplicated
     *either way, because it must be modified.*/
    CIn=mxGetDoubles(prhs[0]);
    C=(double*)mxMalloc(numRow*numCol*sizeof(double));
    if(missedDetectByRow) {
        //Duplicate the input array.
        memcpy(C,CIn,numRow*numCol*sizeof(double));
    } else {
        size_t temp, curRow, curCol;

        //Copy the elements in prhs[0] into C in a transposed order.
        for(curCol=0;curCol<numCol;curCol++) {
            const size_t colOffset=numRow*curCol;
            
            for(curRow=0;curRow<numRow;curRow++) {
                const size_t rowOffset=numCol*curRow;
                C[curCol+rowOffset]=CIn[curRow+colOffset];
            }
        }

        temp=numRow;
        numRow=numCol;
        numCol=temp;
    }
    numRowsTrue=numRow;
    //The number of rows of the virtually augmented matrix.
    numRow=numRowsTrue-1+numCol;
    
    /* These will hold the dual variable values. They are all initialized
     * to zero.*/
    gainMATLAB=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    uMATLAB=mxCreateNumericMatrix(numCol,1,mxDOUBLE_CLASS,mxREAL);
    vMATLAB=mxCreateNumericMatrix(numRow,1,mxDOUBLE_CLASS,mxREAL);
    gainPtr=mxGetDoubles(gainMATLAB);
    
    //Allocate the temporary space used by the assign2DMissedDetectC
    //function. We also allocate an additional 2*numCol*sizeof(ptrdiff_t)
    //bytes to hold the tuples output, which will be converted into doubles
    //in other memory prior to returning it to Matlab, because floating
    //point doubles are the default format of values in Matlab.
    tempBuffer=mxMalloc(assign2DMissedDetectCBufferSize(numRowsTrue,numCol)+2*numCol*sizeof(ptrdiff_t));
    tuples=(ptrdiff_t*)tempBuffer;

    {
        void *bufferStart=(void*)(tuples+2*numCol);

        isInfeasible=assign2DMissedDetectC(maximize, C, gainPtr, tuples, bufferStart, mxGetDoubles(uMATLAB), mxGetDoubles(vMATLAB), numRowsTrue, numCol);
    }

    //If there is no feasible solution.
    if(isInfeasible) {
        mxFree(tempBuffer);
        
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        if(nlhs>1){
            plhs[2]=mxCreateDoubleScalar(-1);

            if(nlhs>2) {
                plhs[3]=uMATLAB;

                if(nlhs>3){
                    plhs[4]=vMATLAB;
                }
            }
        }
        return;
    } else {
        //Get rid of the entries in v that are related to the
        //"unconstrained" column.
        mxSetM(vMATLAB,numRowsTrue-1);
        
        /*Convert the C indices in tuples into indices for MATLAB. Also,
         *if missedDetectByRow==false, then the ordering of the indices in
         *tuples needs to be flipped.
         */
        if(missedDetectByRow==false) {
            size_t curCol;
            mxArray *temp;
            
            for(curCol=0;curCol<numCol;curCol++) {
                const size_t idx=2*curCol;
                ptrdiff_t tempP;
                
                tempP=tuples[idx];
                tuples[idx]=tuples[idx+1]+1;
                tuples[idx+1]=tempP+1;
            }

            temp=uMATLAB;
            uMATLAB=vMATLAB;
            vMATLAB=temp;
        } else {
            size_t curCol;

            for(curCol=0;curCol<numCol;curCol++) {
                const size_t idx=2*curCol;

                tuples[idx]++;
                tuples[idx+1]++;
            }
        }

        plhs[0]=ptrDiffTMat2MatlabDoubles(tuples,2,numCol);
        mxFree(tempBuffer);

        if(nlhs>1) {
            plhs[1]=gainMATLAB;
            if(nlhs>2) {
                plhs[2]=uMATLAB;
                if(nlhs>3) {
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
