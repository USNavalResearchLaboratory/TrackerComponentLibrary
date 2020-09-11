/**FINDFIRSTMAX Given an array of sorted vector of values in increasing
*               order, which might be full of duplicates, find the first
*               occurrence of the maximum value. An implementation in C++.
*
*INPUTS: arr A vector of values, with possible repeats, sorted in
*            increasing order.
*
*OUTPUTS: idx The index of the first occurrence of the maximum value (last)
*             element in arr.    
*
* The easiest way to find the first occurrence of the maximum value is to
* scan downward in the array from the end until a value different from the
* maximum value is found. However, that approach has a bad worst-case
* accuracy if all of the elements in the array have the same value. Thus,
* this method first checks for the easiest solutions, and then uses
* something akin to a binary search technique to find the first occurrence
* of the maximum value.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* idx=findFirstMax(arr);
*
*December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "mathFuncs.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numRow, numCol,numPoints, foundIdx;
    double *arr;
    mxArray *retIdx;
    
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nrhs>1) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    /*Verify the validity of the vector to search.*/
    checkRealDoubleArray(prhs[0]);
    numRow = mxGetM(prhs[0]);
    numCol = mxGetN(prhs[0]);
    if(numRow>1&&numCol>1){
        mexErrMsgTxt("Too many dimensions in the array provided.");
        return;
    } else {
        if(numRow>numCol) {
            numPoints=numRow;
        } else {
            numPoints=numCol;
        }
    }
    
    arr=mxGetDoubles(prhs[0]);
    foundIdx= findFirstMaxCPP(arr,numPoints);
    
    //Set the return value
    retIdx=allocUnsignedSizeMatInMatlab(1,1);
    if(sizeof(size_t)==4) {//32 bit
        *reinterpret_cast<size_t *>(mxGetUint32s(retIdx))=foundIdx+1;
    } else {//64 bit
        *reinterpret_cast<size_t *>(mxGetUint64s(retIdx))=foundIdx+1;
    }

    plhs[0]=retIdx;
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
