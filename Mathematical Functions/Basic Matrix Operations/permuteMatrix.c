/**PERMUTEMATRIX Given an n1Xn2X...nS matrix, rearrange the dimensions of C
 *     to be in the order specified by order. All elements of order must be
 *     unique, real, positive, integer values from 1 to S. This functions
 *     in the same manner as Matlab's permute function. The permutation
 *     performed by this function is not done in-place. This type of
 *     permutation is a type of tensor matrix transpose.
 *
 *INPUTS: C An n1Xn2X...XnS matrix.
 *   order The length S order vector.
 *
 *OUTPUTS: CPerm C with its dimensions permuted according to order.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *CPerm=permuteMatrix(C,order)
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This is for input validation*/
#include "MexValidation.h"

#include "basicMatOps.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *orderMatlab;
    void *CPermReal, *CPermImag=NULL;
    void *CReal, *CImag;
    mxArray *CPermMATLAB;
    const size_t *nVals;
    size_t *order, *nValsNew, *nValsExtended=NULL;
    size_t S, numOrderEls, numElsInC;
    size_t i;
    void *tempBuffer, *curBufferSpot;
    mxClassID CClass;
    size_t classDataSize;

    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    S=mxGetNumberOfDimensions(prhs[0]);
    numElsInC=mxGetNumberOfElements(prhs[0]);
    numOrderEls=mxGetNumberOfElements(prhs[1]);

    if(mxIsEmpty(prhs[1])&&!mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("order cannot be empty with a non-empty C.");
        return;
    }
    
    if(numOrderEls<S){
        mexErrMsgTxt("The ordering vector must be >=ndims(C).");
        return;
    }
    
    if(mxIsEmpty(prhs[0])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//CPerm=[];
        return;
    }
    
    checkRealDoubleArray(prhs[1]);
    
    CClass=mxGetClassID(prhs[0]);
    switch(CClass) {
        case mxLOGICAL_CLASS:
             classDataSize=sizeof(mxLogical);
             break;
        case mxCHAR_CLASS:
            classDataSize=sizeof(char);
            break;
        case mxDOUBLE_CLASS:
            classDataSize=sizeof(double);
            break;
        case mxSINGLE_CLASS:
            classDataSize=sizeof(float);
            break;
        case mxINT8_CLASS:
            classDataSize=sizeof(int8_t);
            break;
        case mxUINT8_CLASS:
            classDataSize=sizeof(uint8_t);
            break;
        case mxINT16_CLASS:
            classDataSize=sizeof(int16_t);
            break;
        case mxUINT16_CLASS:
            classDataSize=sizeof(uint16_t);
            break;
        case mxINT32_CLASS:
            classDataSize=sizeof(int32_t);
            break;
        case mxUINT32_CLASS:
            classDataSize=sizeof(uint32_t);
            break;
        case mxINT64_CLASS:
            classDataSize=sizeof(int64_t);
            break;
        case mxUINT64_CLASS:
            classDataSize=sizeof(uint64_t);
            break;
        case mxUNKNOWN_CLASS:
        case mxFUNCTION_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("C is of an unsupported class.");
            return;
    }
    
    CReal=mxGetData(prhs[0]);
    CImag=mxGetImagData(prhs[0]);

    orderMatlab=mxGetPr(prhs[1]);
    for(i=0;i<S;i++) {
        if(orderMatlab[i]<1||orderMatlab[i]>numOrderEls) {
            mexErrMsgTxt("The values in order must be from 1 to numOrderEls.");
            return;
        }
    }

    //If numOrderEls>S, then we will create a new nVals and pad it with 1s
    //at the end.
    nVals=mxGetDimensions(prhs[0]);
    
    //Allocate temporary space plus the space needed for order and
    //nValsNew Also, if  numOrderEls>S, then space must be allocated for a
    //new nVals vector.
    if(numOrderEls>S) {
        const size_t SOrig=S;
        
        S=numOrderEls;
        tempBuffer=mxMalloc(numElsInC*sizeof(double)+permuteMatrixCBufferSize(S)+3*S*sizeof(size_t));
        order=(size_t*)tempBuffer;
        nValsNew=order+S;
        nValsExtended=nValsNew+S;
        curBufferSpot=(void*)(nValsExtended+S);
        
        //Copy nVals into nValsExtended and add the 1s at the end.
        memcpy(nValsExtended,nVals,SOrig*sizeof(size_t));
        
        for(i=SOrig;i<S;i++) {
            nValsExtended[i]=1;
        }
        nVals=nValsExtended;
    }
    else {
        tempBuffer=mxMalloc(numElsInC*sizeof(double)+permuteMatrixCBufferSize(S)+2*S*sizeof(size_t));
        order=(size_t*)tempBuffer;
        nValsNew=order+S;
        curBufferSpot=(void*)(nValsNew+S);
    }
            
    for(i=0;i<S;i++) {//Copy order into a format usable by the function.
        size_t j;
        order[i]=(size_t)(orderMatlab[i]-1);
        
        for(j=0;j<i;j++) {
            if(order[j]==order[i]) {
                mexErrMsgTxt("The values in order must be unique.");
                return;
            }
        }
    }

    //Allocate the return matrix. We initially give it the original
    //dimensions. the dimensions will be set to the permuted dimensions
    //after calling permuteMatrixC.
    if(CImag==NULL) {//If it is real.
        CPermMATLAB=mxCreateNumericArray(S,nVals,CClass,mxREAL);
        CPermReal=mxGetData(CPermMATLAB);
    } else {//If it is imaginary.
        CPermMATLAB=mxCreateNumericArray(S,nVals,CClass,mxCOMPLEX);
        CPermReal=mxGetData(CPermMATLAB);
        CPermImag=mxGetImagData(CPermMATLAB);
    }

    if(CReal==NULL) {
        mexErrMsgTxt("The real portion of the data is missing.");
        return;
    }
    
    permuteMatrixC(S,nVals, nValsNew, CPermReal, CReal, classDataSize, curBufferSpot, order);
    if(CImag!=NULL) {
        permuteMatrixC(S,nVals, nValsNew, CPermImag, CImag, classDataSize, curBufferSpot, order);
    }

    if(mxSetDimensions(CPermMATLAB,nValsNew,S)) {
        mexErrMsgTxt("An error occurred setting return array dimension vector.");
        return;
    }
    
    mxFree(tempBuffer);
    
    plhs[0]=CPermMATLAB;    
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
