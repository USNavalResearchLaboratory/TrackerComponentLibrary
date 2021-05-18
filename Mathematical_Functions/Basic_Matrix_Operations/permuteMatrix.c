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
#include "matrix.h"
/*This is for input validation*/
#include "MexValidation.h"
#include "basicMatOps.h"
#include <stdbool.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *orderMatlab;
    void *CPerm;
    void *CData;
    mxArray *CPermMATLAB;
    const size_t *nVals;
    size_t *order, *nValsNew, *nValsExtended=NULL;
    size_t S, numOrderEls;
    size_t i;
    void *tempBuffer, *curBufferSpot;
    mxClassID CClass;
    size_t classDataSize;
    bool isComplex;

    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    S=mxGetNumberOfDimensions(prhs[0]);
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
    orderMatlab=mxGetDoubles(prhs[1]);
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
    //nValsNew. Also, if  numOrderEls>S, then space must be allocated for a
    //new nVals vector.
    if(numOrderEls>S) {
        const size_t SOrig=S;
        
        S=numOrderEls;
        tempBuffer=mxMalloc(permuteMatrixCBufferSize(S)+3*S*sizeof(size_t));
        order=(size_t*)tempBuffer;//Length S
        nValsNew=order+S;//Length S
        nValsExtended=nValsNew+S;//Length S
        curBufferSpot=(void*)(nValsExtended+S);//Length permuteMatrixCBufferSize(S).
        
        //Copy nVals into nValsExtended and add the 1s at the end.
        memcpy(nValsExtended,nVals,SOrig*sizeof(size_t));
        
        for(i=SOrig;i<S;i++) {
            nValsExtended[i]=1;
        }
        nVals=nValsExtended;
    }
    else {
        tempBuffer=mxMalloc(permuteMatrixCBufferSize(S)+2*S*sizeof(size_t));
        order=(size_t*)tempBuffer;//Length S
        nValsNew=order+S;//Length S
        curBufferSpot=(void*)(nValsNew+S);// Length permuteMatrixCBufferSize(S).
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

    //Determine the data size and allocate the return matrix.
    isComplex=mxIsComplex(prhs[0]);
    CClass=mxGetClassID(prhs[0]);
    
    //Allocate the return matrix. We initially give it the original
    //dimensions. The dimensions will be set to the permuted dimensions
    //after calling permuteMatrixC.
    if(isComplex) {
        CPermMATLAB=mxCreateNumericArray(S,nVals,CClass,mxCOMPLEX);    
    } else {
        CPermMATLAB=mxCreateNumericArray(S,nVals,CClass,mxREAL);
    }

    switch(CClass) {
        case mxLOGICAL_CLASS:
             if(isComplex) {
                mexErrMsgTxt("Complex logical variables are not supported.");
                return;
             }
             classDataSize=sizeof(mxLogical);
             CData=(void*)mxGetLogicals(prhs[0]);
             CPerm=(void*)mxGetLogicals(CPermMATLAB);
             break;
        case mxCHAR_CLASS:
             if(isComplex) {
                mexErrMsgTxt("Complex character variables are not supported.");
                return;
             }
            classDataSize=sizeof(mxChar);
            CData=(void*)mxGetChars(prhs[0]);
            CPerm=(void*)mxGetChars(CPermMATLAB);
            break;
        case mxDOUBLE_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexDouble);
                CData=(void*)mxGetComplexDoubles(prhs[0]);
                CPerm=(void*)mxGetComplexDoubles(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxDouble);
                CData=(void*)mxGetDoubles(prhs[0]);
                CPerm=(void*)mxGetDoubles(CPermMATLAB);
            }
            break;
        case mxSINGLE_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexSingle);
                CData=(void*)mxGetComplexSingles(prhs[0]);
                CPerm=(void*)mxGetComplexSingles(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxSingle);
                CData=(void*)mxGetSingles(prhs[0]);
                CPerm=(void*)mxGetSingles(CPermMATLAB);
            }
            break;
        case mxINT8_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexInt8);
                CData=(void*)mxGetComplexInt8s(prhs[0]);
                CPerm=(void*)mxGetComplexInt8s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxInt8);
                CData=(void*)mxGetInt8s(prhs[0]);
                CPerm=(void*)mxGetInt8s(CPermMATLAB);
            }
            break;
        case mxUINT8_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexUint8);
                CData=(void*)mxGetComplexUint8s(prhs[0]);
                CPerm=(void*)mxGetComplexUint8s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxUint8);
                CData=(void*)mxGetUint8s(prhs[0]);
                CPerm=(void*)mxGetUint8s(CPermMATLAB);
            }
            break;
        case mxINT16_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexInt16);
                CData=(void*)mxGetComplexInt16s(prhs[0]);
                CPerm=(void*)mxGetComplexInt16s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxInt16);
                CData=(void*)mxGetInt16s(prhs[0]);
                CPerm=(void*)mxGetInt16s(CPermMATLAB);
            }
            break;
        case mxUINT16_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexUint16);
                CData=(void*)mxGetComplexUint16s(prhs[0]);
                CPerm=(void*)mxGetComplexUint16s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxUint16);
                CData=(void*)mxGetUint16s(prhs[0]);
                CPerm=(void*)mxGetUint16s(CPermMATLAB);
            }
            break;
        case mxINT32_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexInt32);
                CData=(void*)mxGetComplexInt32s(prhs[0]);
                CPerm=(void*)mxGetComplexInt32s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxInt32);
                CData=(void*)mxGetInt32s(prhs[0]);
                CPerm=(void*)mxGetInt32s(CPermMATLAB);
            }
            break;
        case mxUINT32_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexUint32);
                CData=(void*)mxGetComplexUint32s(prhs[0]);
                CPerm=(void*)mxGetComplexUint32s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxUint32);
                CData=(void*)mxGetUint32s(prhs[0]);
                CPerm=(void*)mxGetUint32s(CPermMATLAB);
            }
            break;
        case mxINT64_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexInt64);
                CData=(void*)mxGetComplexInt64s(prhs[0]);
                CPerm=(void*)mxGetComplexInt64s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxInt64);
                CData=(void*)mxGetInt64s(prhs[0]);
                CPerm=(void*)mxGetInt64s(CPermMATLAB);
            }
            break;
        case mxUINT64_CLASS:
            if(isComplex) {
                classDataSize=sizeof(mxComplexUint64);
                CData=(void*)mxGetComplexUint64s(prhs[0]);
                CPerm=(void*)mxGetComplexUint64s(CPermMATLAB);
            } else {
                classDataSize=sizeof(mxUint64);
                CData=(void*)mxGetUint64s(prhs[0]);
                CPerm=(void*)mxGetUint64s(CPermMATLAB);
            }
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

    permuteMatrixC(S,nVals, nValsNew, CPerm, CData, classDataSize, curBufferSpot, order);

    if(mxSetDimensions(CPermMATLAB,nValsNew,S)) {
        mexErrMsgTxt("An error occurred setting return array dimension vector.");
        return;
    }
    
    mxFree(tempBuffer);
    plhs[0]=CPermMATLAB;    
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
