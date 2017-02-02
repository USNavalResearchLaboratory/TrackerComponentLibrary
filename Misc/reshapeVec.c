/**reshapeVec Reshape a matrix given a vector of the sizes of the
 *            dimensions. This is essentially the same as the reshape
 *            function in matlab, except the reshape function only takes
 *            the dimensions as multiple inputs, not as a vector. This
 *            function is useful in implementing functions where the total
 *            number of indices of a matrix is not known a priori.
 *
 *INPUTS:   A      A matrix whose dimensions are to be changed.
 *          dimVec A numDimX1 or 1XnumDim vector of the dimensions that are
 *                 to be set for the matrix A. A scalar is considered 1X1
 *                 in Matlab, thus normally numDim=2 is the minimum. Here,
 *                 if dimVec contains only one element, a trailing one is
 *                 automatically added to dimVec.
 *
 *OUTPUTS: A  A copy of the input matrix A with its dimensions changed to
 *            those of dimVec
 *
 *Note that this function, like Matlab's reshape function, automatically
 *removes any trailing singleton dimensions (beyond the second dimension).
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *A=reshapeVec(A,dimVec);
 *
 *August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize numEls, numIdx,extraIdx=0;
    mxArray *AOut;
    mxClassID dimVecType;
    mwSize *dims;
    size_t i;
    
    if(nrhs!=2){
        mexErrMsgTxt("Wrong number of inputs");
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
    }

    //Determine the total number of elements in the input matrix
    numEls=mxGetNumberOfElements(prhs[0]);

    //The number of indices for the reshaped matrix.
    numIdx=mxGetNumberOfElements(prhs[1]);
    if(numIdx<1) {
        mexErrMsgTxt("A minimum of two indices must be specified in dimVec.");
    } else if(numIdx==1) {
        //Matlab requires at least two indices, so pad an extra singleton index.
        extraIdx=1;
    } else {
        extraIdx=0;
    }
    
    //Allocate space for the vector of indices to be passed to
    dims=mxMalloc((numIdx+extraIdx)*sizeof(mwSize));
    
    //Copy the passed vector of indices into a vector that can be passed to
    //the mxSetDimensions function.
    dimVecType=mxGetClassID(prhs[1]);
    switch(dimVecType) {
        case mxDOUBLE_CLASS:
        {
            double *dimsVec=(double*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *dimsVec=(float*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *dimsVec=(int8_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *dimsVec=(uint8_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *dimsVec=(int16_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *dimsVec=(uint16_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *dimsVec=(int32_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *dimsVec=(uint32_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *dimsVec=(int64_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *dimsVec=(uint64_T*)mxGetData(prhs[1]);
            for(i=0;i<numIdx;i++) {
               dims[i]=(mwSize)dimsVec[i];
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mxFree(dims);//Free allocated memory
            mexErrMsgTxt("The dimensions vector is an unsupported type.");
    }
    //If an extra signleton index is added.
    if(extraIdx!=0) {
        dims[numIdx]=1;
        numIdx++;
    }
    
    //Check to make sure that the total number of elements is not changing.
    {
        mwSize newNumEls=dims[0];
        
        for(i=1;i<numIdx;i++) {
            newNumEls*=dims[i];
        }
        
        if(newNumEls!=numEls) {
            mxFree(dims);//Free allocated memory
            mexErrMsgTxt("The number of elements in matrix A cannot change."); 
        }
    }
    
    //Duplicate the input matrix as we cannot modify in place.
    AOut=mxDuplicateArray(prhs[0]);
    //Set the number of dimensions.
    if(mxSetDimensions(AOut,dims,numIdx)!=0) {
        mxFree(dims);
        mxDestroyArray(AOut);
        mexErrMsgTxt("An error occurred setting the dimensions.");
    } else {
        mxFree(dims);
        plhs[0]=AOut;
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
