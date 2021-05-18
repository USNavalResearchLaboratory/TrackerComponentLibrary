/**BINSEARCHDOUBLES A C-code implementation for performing  a binary search
*             for the value key in the vector vec, which is a vector of
*             double-precision floating point numbers that has been sorted
*             such that the elements are in increasing order. The optional
*             choice parameter sets what is returned if key is not in vec.
*             This differs from the binSearch function in that it has been
*             implemented in C for speed and thus assumes that the inputs
*             are doubles, whereas binSearch will work with a variety of
*             data types.
*
*INPUTS: vec A vector with elements sorted in increasing order.
*        key The value that one wishes to find in the vector vec.
*    choice  An optional parameter that determines what is returned if key
*            is not found.
*            0 means return the closest value.
*            1 means return the next lower value if there is one, otherwise
*              return the lowest value in vec.
*            2 means return the next higher value if there is one,
*              otherwise return the highest value in vec.
*
*OUTPUTS: val Either key, if found, or a value nearby as determined by the
*             parameter closest.
*         idx The index of val in the vector vec.
*
* This is just a basic binary search. The search space is cut in half each
* time. In some cases, such as are elaborated in [1] the Fibonacci search
* can be faster. However, for most problems, there is little difference
* between a binary and Fibonacci search and in [2], it is shown that the
* binary search has a degree of optimality over the more complicated
* Fibonacci search in many instances.
*
* If one does not require as fine a choice of the output for the case when
* a key is not found, then the bsearch function in stdlib could just be
* used.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [val, idx]=binSearchDoubles(vec,key,choice);
*
*REFERENCES:
* [1] S. Nishihara and H. Nishino, "Binary search revisited: Another
*     advantage of Fibonacci search," IEEE Transactions on Computers, vol.
*     C-36, no. 9, pp. 1132-1135, Sep. 1987.
* [2] K. J. Overholt, "Optimal binary search methods," BIT Numerical
*     Mathematics, vol. 13, no. 1, pp. 84-91, 1973.
*
*December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
/*This is needed for fabs*/
#include <math.h>
/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "mathFuncsC.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int choice;
    size_t numRow, numCol,numPoints,foundIdx;
    double *vec;
    double key;
    mxArray *retVal;
    
    if(nrhs<2){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Too many outputs.");
        return;
    }
    
    if(nrhs==2){
        choice=0;
    } else {
        choice=getIntFromMatlab(prhs[2]);
        
        if(choice <0 || choice >2){
            mexErrMsgTxt("Invalid maximization parameter given.");
            return;
        }
    }
    
    /*Verify the validity of the vector to sort.*/
    checkRealDoubleArray(prhs[0]);
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    if(numRow>1&&numCol>1){
        mexErrMsgTxt("Too many dimensions in the search vector.");
        return;
    } else {
        if(numRow>numCol) {
            numPoints=numRow;
        } else {
            numPoints=numCol;
        }
    }

    vec=mxGetDoubles(prhs[0]);
    key=getDoubleFromMatlab(prhs[1]);
    
    foundIdx=binSearchC(numPoints, vec, key, choice);
    
    retVal=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetDoubles(retVal)=vec[foundIdx];
    plhs[0]=retVal;

    if(nlhs>1){
        mxArray *retIdx;

        retIdx=allocUnsignedSizeMatInMatlab(1,1);
        if(sizeof(size_t)==4) {//32 bit
            *(size_t *)mxGetUint32s(retIdx)=foundIdx+1;
        } else {//64 bit
            *(size_t *)mxGetUint64s(retIdx)=foundIdx+1;
        }
        
        plhs[1]=retIdx;
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
