/**ASSIGN2D A C-code (for Matlab) implementation of the shortest path
 *          assignment algorithm to solve the two-dimensional assignment
 *          problem with a rectangular cost matrix C. This implementation
 *          scans the cost matrix by row rather than by column. See the
 *          comments to the Matlab implementation for more details.
 *
 *INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
 *         where the largest finite element minus the smallest element is a
 *         finite quantity (does not overflow) when performing minimization
 *         and where the smallest finite element minus the largest
 *         element is finite when performing maximization. Forbidden
 *         assignments can be given costs of +Inf for minimization and -Inf
 *         for maximization.
 * maximize If true, the minimization problem is transformed into a
 *          maximization problem. The default if this parameter is omitted
 *          is false.
 *
 *OUTPUTS: col4row A numRowX1 vector where the entry in each element is an
 *                assignment of the element in that row to a column. 0
 *                entries signify unassigned rows. If the problem is
 *                infeasible, this is an empty matrix.
 *        row4col A numColX1 vector where the entry in each element is an
 *                assignment of the element in that column to a row. 0
 *                entries signify unassigned columns. If the problem is
 *                infeasible, this is an empty matrix.
 *           gain The sum of the values of the assigned elements in C. If
 *                the problem is infeasible, this is -1.
 *            u,v The dual variable for the columns and for the rows. Note
 *                that if C contains any negative entries for minimization
 *                or positive entries for maximization, then these are from
 *                version of C transformed to be all positive
 *                (minimization) or negative (maximization).
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
 * [col4row,row4col,gain,u,v]=assign2D(C,maximize)
 *
 *October 2013, notably modified December 2017 David F. Crouse, Naval
 *Research Laboratory, Washington D.C.
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
    double gain, *CIn;
    size_t numRow, numCol, i;
    ptrdiff_t *col4row, *row4col;
    mxArray *uMATLAB, *vMATLAB;
    bool maximize=false;
    bool isInfeasible;
    void *tempBuffer;
    
    if(nrhs<1) {
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>2) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nrhs==2&&!mxIsEmpty(prhs[1])){
        maximize=getBoolFromMatlab(prhs[1]);
    }

    if(nlhs>5) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    /*Verify the validity of the assignment matrix.*/
    checkRealDoubleArray(prhs[0]);

    /* Get the dimensions of the input data. It is assumed that the matrix
     * is not so large in either dimension as to cause an overflow when
     * using a SIGNED integer data type.*/
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);

    /* Transpose the matrix, if necessary, so that the number of rows is
     * >= the number of columns. The matrix must be duplicated either way,
     * because it will be modified.*/
    CIn=mxGetDoubles(prhs[0]);
  
    /* These will hold the dual variable values. They are all initialized
     * to zero.*/
    uMATLAB=mxCreateNumericMatrix(numCol,1,mxDOUBLE_CLASS,mxREAL);
    vMATLAB=mxCreateNumericMatrix(numRow,1,mxDOUBLE_CLASS,mxREAL);

    //Allocate the temporary space needed by the assign2DC function. We 
    //also assign an additional (numRow+numCol)*sizeof(ptrdiff_t) values to
    //hold the col4row and row4col outputs of the function. These outputs
    //will be converted into doubles in other memory prior to returning
    //them to Matlab, because floating point doubles are the default
    //format of values in Matlab.
    tempBuffer=mxMalloc(assign2DSimpBufferSize(numRow,numCol)+(numRow+numCol)*sizeof(ptrdiff_t));

    col4row=(ptrdiff_t*)tempBuffer;
    row4col=col4row+numRow;
    
    {
        void* bufferStart=(void*)(row4col+numCol);
        
        isInfeasible=assign2DSimp(maximize,CIn,&gain, col4row, row4col,bufferStart,mxGetDoubles(uMATLAB), mxGetDoubles(vMATLAB), numRow, numCol);
    }
    
    //If there is no feasible solution.
    if(isInfeasible){
        mxFree(tempBuffer);
        
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
    } else {
        /* Convert the C indices into indices for MATLAB.*/
        for(i=0;i<numRow;i++){
            col4row[i]++;
        }
        
        for(i=0;i<numCol;i++){
            row4col[i]++;
        }

        plhs[0]=ptrDiffTMat2MatlabDoubles(col4row,numRow,1);
        if(nlhs>1) {
            plhs[1]=ptrDiffTMat2MatlabDoubles(row4col,numCol,1);
            
            for(i=0;i<numCol;i++){
                row4col[i]=row4col[i]+1;
            }

            if(nlhs>2) {
                plhs[2]=mxCreateDoubleScalar(gain);
                
                if(nlhs>3) {
                    plhs[3]=uMATLAB;
                    
                    if(nlhs>4){
                        plhs[4]=vMATLAB;
                    }
                }
            }
        }
        
        mxFree(tempBuffer);
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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
*/

