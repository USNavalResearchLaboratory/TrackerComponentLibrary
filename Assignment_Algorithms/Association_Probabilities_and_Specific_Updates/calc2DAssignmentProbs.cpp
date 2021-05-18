/**CALC2DASSIGNMENTPROBS  A C++ implementation of an algorithm that,
*                given a matrix of all-positive likelihood ratios,
*                determines the probability that each row is assigned to
*                each target. Whereas the assign2D algorithm provides the
*                optimal assignment of all rows to columns, this provides
*                the probability of each assignment. The diagAugment
*                parameter allows for a more efficient algorithm to be used
*                in the case that for a numRow X numCol matrix where
*                numCol>=numRow, the last numRow columns are a weighted
*                identity matrix. Such a matrix strcture arises when
*                performing target-measurement data association with missed
*                detection hypotheses. In such an instance, the rows are
*                targets, the first numCol-numRow columns correspond to
*                measurements and the final numRow columns have missed
*                detection likelihoods.
*
*INPUTS: A A matrix of positive likelihoods or likelihood ratios (NOT log-
*          likelihood ratios). If diagAugment=true, then A is a
*          numTarX(numMeas+numTar) matrix of all-positive likelihoods or
*          likelihood ratios for assigning the target specified by the row
*          to the measurement specified by the column. Columns > numMeas
*          hold missed-detection likelihoods/ likelihood ratios. Thus, off-
*          diagonal terms for columns > numMeas should be set to 0 and the
*          diagonal terms set to the costs of a missed detection for each
*          given target.
* diagAugment A boolean variable indicating whether the probabilties
*          should be found for a general assignment problem or assuming A
*          ends with a diagonal submatrix of missed detection
*          probabilities. The default if omitted is false (the general
*          problem). Setting diagAugment to true changes the shape of the
*          output. The default if omitted is false. See the description of
*          the output beta for how this affects the output.
*
*OUTPUTS: beta If diagAugment is omitted or false, then beta has the same
*              dimensionality as A and hold the probability of assigning
*              each row to each column in the traditional 2D assignment
*              problem (Each row must be assigned to at least one column,
*              each column can be assigned to at most one row, or vice
*              versa if there are more rows than columns. If diagAugment
*              is true, then beta is a numTar X (numMeas+1) matrix of
*              probabilities of assigning the target given by the row to
*              the measurement given by the column. The final column is a
*              set of missed detection probabilities. Thus, in this case,
*              beta has fewer columns than A, because the final numRow
*              columns are collapsed into one column on return.
*
*Details of the algorithm are given in the comments to the Matlab
*implementation.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*beta=calc2DAssignmentProbs(A);
*or
*beta=calc2DAssignmentProbs(A,diagAugment);
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "MexValidation.h"
/*This header is required by Matlab*/
#include "mex.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool diagAugment=false;
    size_t numRow, numCol,numBetaCols;
    mxArray *betaMatlab;
    double *A, *beta;

    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>2) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    
    //If an empty matrix is passed, then return an empty matrix.
    if(numRow==0||numCol==0) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

    checkRealDoubleArray(prhs[0]);

    if(nrhs>1) {
        diagAugment=getBoolFromMatlab(prhs[1]);
    }
    
    //Get the matrix.
    A=mxGetDoubles(prhs[0]);
    
    //If we are solving the general probability problem.
    if(diagAugment==false) {
        //If we are here, then we just want general assignment probabilities,
        //not specialized to target tracking applications.
        size_t *buffer,*rows2Keep,*cols2Keep;
        size_t *buff4PermFunc;//Buffer for the permanent function.
        const size_t numRowsKept=numRow-1;
        const size_t numColsKept=numCol-1;
        size_t curRow,curCol;
        
        //Allocate the return values.
        numBetaCols=numCol;
        betaMatlab=mxCreateDoubleMatrix(numRow,numBetaCols,mxREAL);
        beta=mxGetDoubles(betaMatlab);
        
        //Allocate the space for the indices of the rows and columns that
        //are kept in the submatrix to be passed to the permCPPSkip
        //function.
        buffer=new size_t[numRowsKept+2*numColsKept];
        rows2Keep=buffer;
        cols2Keep=rows2Keep+numRowsKept;
        buff4PermFunc=cols2Keep+numColsKept;
        
        for(curRow=0;curRow<numRow;curRow++) {
             size_t i;
             //Fill in the indices of the rows to keep.
             for(i=curRow;i<numRowsKept;i++) {
                rows2Keep[i]=i+1;
             }

             for(curCol=0;curCol<numCol;curCol++){
                 double ati=A[curRow+curCol*numRow];
                 //The A matrix removing the row and column corresponding
                 //to the selected association.
                 //Fill in the indices of the columns to keep.
                 for(i=curCol;i<numColsKept;i++) {
                    cols2Keep[i]=i+1;
                 }

                 beta[curRow+curCol*numRow]=ati*permCPPSkip(A,numRow,rows2Keep, cols2Keep, numRowsKept, numColsKept, buff4PermFunc);
                 cols2Keep[curCol]=curCol;
             }
             rows2Keep[curRow]=curRow;
        }
        
        //Delete the buffers
        delete buffer;
    }else {
        //This is the case where we are solving for target-measurement
        //assignment probabilities with missed detections.
        size_t *buffer,*rows2Keep,*cols2Keep;
        size_t *buff4PermFunc;//Buffer for the permanent function.
        const size_t numRowsKept=numRow-1;
        const size_t numColsKept=numCol-1;
        size_t numTar,numMeas;
        size_t curTar,curMeas;
        
        if(numCol<numRow) {
            mexErrMsgTxt("The number of columns cannot be less than the number of rows when diagAugment is true.");
        }

        numTar=numRow;
        numMeas=numCol-numTar;
        numBetaCols=numMeas+1;

        //Allocate the return values.
        betaMatlab=mxCreateDoubleMatrix(numRow,numBetaCols,mxREAL);
        beta=mxGetDoubles(betaMatlab);

        //Allocate the space for the indices of the rows and columns
        //that are kept in the submatrix to be passed to the
        //permCPPSkip function.
        buffer=new size_t[numRow+2*numCol];
        rows2Keep=buffer;
        cols2Keep=rows2Keep+numRowsKept;
        buff4PermFunc=cols2Keep+numColsKept;

        for(curTar=0;curTar<numTar;curTar++) {
            size_t i;
            double ati;
            //Fill in the indices of the rows to keep.
            for(i=curTar;i<numRowsKept;i++) {
                rows2Keep[i]=i+1;
            }
            
            for(curMeas=0;curMeas<numMeas;curMeas++) {
                //The measurement hypotheses
                ati=A[curTar+curMeas*numRow];
                //The A matrix removing the row and column corresponding
                //to the selected association.
                //Fill in the indices of the columns to keep.
                for(i=curMeas;i<numColsKept;i++) {
                    cols2Keep[i]=i+1;
                }

                beta[curTar+curMeas*numRow]=ati*permCPPSkip(A,numRow,rows2Keep, cols2Keep, numRowsKept, numColsKept, buff4PermFunc);
                cols2Keep[curMeas]=curMeas;
            }
            //The missed detection hypothesis
            curMeas=numMeas+curTar;
            ati=A[curTar+curMeas*numRow];
            //Fill in the indices of the columns to keep.
            for(i=numMeas;i<curMeas;i++) {
                cols2Keep[i]=i;
            }
            for(i=curMeas;i<numColsKept;i++) {
                cols2Keep[i]=i+1;
            }
            
            beta[curTar+numMeas*numRow]=ati*permCPPSkip(A,numRow,rows2Keep, cols2Keep, numRowsKept, numColsKept, buff4PermFunc);
            rows2Keep[curTar]=curTar;
        }
        //Delete the buffers
        delete buffer;
    }
    
    //Normalize the return values.
    //It is faster to normalize the betas this way then to compute the
    //normalization constant by finding the permanent of the entire A
    //matrix.
    {
        double sumVal;
        size_t curRow,curCol;
        
        for(curRow=0;curRow<numRow;curRow++) {
            sumVal=0;
            //Compute the sum across the columns
            for(curCol=0;curCol<numBetaCols;curCol++) {
                sumVal+=beta[curRow+curCol*numRow];
            }
            //Normalize the value across the columns.
            for(curCol=0;curCol<numBetaCols;curCol++) {
                size_t idx=curRow+curCol*numRow;
                beta[idx]/=sumVal;
            }
        }
    }
 
    //Set the return values.
    plhs[0]=betaMatlab;
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
