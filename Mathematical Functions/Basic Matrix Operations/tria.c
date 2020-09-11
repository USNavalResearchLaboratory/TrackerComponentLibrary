/*TRIA Square root matrix triangularization. Given a rectangular square
*     root matrix, obtain a lower-triangular square root matrix that is
*     square.
*
*INPUTS: A A numRowXnumCol matrix that is generally not square.
*
*OUTPUTS: S A lower-triangular matrix such that S*S'=A*A'. If
*           numCol>=numRow, then S is a square numRowXnumRow matrix.
*           Otherwise, S is a numRowXnumCol matrix.
*
*This is the tria function needed for various steps in the cubature Kalman
*filter and the square root Kalman filter. It is described in [1]. It has
*been slightly modified from the paper so that the diagonal elements remain
*positive.
*
*The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
* S=tria(A)
*
*REFERENCES:
*[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
*    bistatic measurements," IEEE Aerospace and Electronic Systems
*    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
*
*April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
#include "MexValidation.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    ptrdiff_t numRows,numCols,lwork,info;
    double *AOrig, *A, *R, *tau, *work;
    size_t tauLength,curRow,curCol;
    mxArray *RMat;
    
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    checkRealDoubleArray(prhs[0]);
    
    if(mxIsEmpty(prhs[0])) {
        plhs[0]=mxDuplicateArray(prhs[0]);
        return;
    }
    
    numRows=(ptrdiff_t)mxGetM(prhs[0]);
    numCols=(ptrdiff_t)mxGetN(prhs[0]);

    //We must duplicate A, because the function dgeqrf modifies it for the
    //return value. However, the first step is to also take a transpose of
    //A. Thus, we shall insert it into the new array in the appropriate
    //order.
    A=mxMalloc((size_t)numRows*(size_t)numCols*sizeof(double));
    AOrig=mxGetDoubles(prhs[0]);
    
    for(curRow=0;curRow<numRows;curRow++) {
        for(curCol=0;curCol<numCols;curCol++) {
            A[curCol+numCols*curRow]=AOrig[curRow+numRows*curCol];
        }
    }
    
    //Swap numRows and numCols.
    {
       ptrdiff_t temp;
       
       temp=numRows;
       numRows=numCols;
       numCols=temp;
    }

    if(numRows<numCols) {
        tauLength=numRows;
    } else {
        tauLength=numCols;
    }
    
    tau=(double*)mxMalloc(tauLength*sizeof(double));

    //The optimal work array
    {
    double workArraySize=0;
    lwork=-1;

    //Run the function with lwork=-1 to query the optimal array size
    dgeqrf(
        &numRows,
        &numCols,
        A,
        &numRows,
        tau,
        &workArraySize,
        &lwork,
        &info
        );
    
    if(info!=0) {
        mexErrMsgTxt("An error occurred determining the optimal work array size.");
    }
    
    lwork=(ptrdiff_t)workArraySize;
    }

    //Run the function
    work=(double*)mxMalloc(lwork*sizeof(double));
    dgeqrf(
        &numRows,
        &numCols,
        A,
        &numRows,
        tau,
        work,
        &lwork,
        &info
        );
    
    mxFree(work);
    mxFree(tau);
    
    if(info!=0) {
        mexErrMsgTxt("An error occurred performing a QR decomposition.");
    }

    //A is now essentially the output of [~,A]=qr(A,0) in Matlab.
    //However, the zero elements have not been set and it has not been
    //shrunk. We want to only keep min(numRows,numCols) rows of A and we
    //want to zero the part that should be zero. That equals tauLength.
    //First, we zero the parts that should be zero.
    for(curRow=1;curRow<tauLength;curRow++) {
        for(curCol=0;curCol<curRow;curCol++) {
           A[curRow+numRows*curCol]=0;
        }
    }
    //As opposed to trying to do the transpose in-place and then resizing
    //the matrix to be returned to Matlab, we just allocate another matrix
    //to return.
    RMat=mxCreateDoubleMatrix(numCols,tauLength,mxREAL);
    R=mxGetDoubles(RMat);
    for(curRow=0;curRow<tauLength;curRow++) {
        for(curCol=0;curCol<numCols;curCol++) {
            R[curCol+numCols*curRow]=A[curRow+numRows*curCol];
        }
    }
    
    mxFree(A);
    //Next, we have to make sure that the diagonal elements are positive.
    //For any diagonal element that is not positive, we flip the sign of
    //the entire column.
    for(curCol=0;curCol<tauLength;curCol++) {
        if(R[curCol+numCols*curCol]<0) {
            for(curRow=curCol;curRow<numCols;curRow++) {
                R[curRow+numCols*curCol]*=-1.0;
            }
        }
    }

    plhs[0]=RMat;
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
