/**GETNEXTGRAYCODE A C++ implementation of a function to return the next
*                 gray code value in the sequence given the current gray
*                 code value. The first code to start the sequence is all
*                 logical false (zeros), which is returned if this function
*                 is called with an empty matrix for the code and nCard=n
*                 (the length of the code). A gray code is a  binary code
*                 such that only one entry changes from 0 to 1 or back each
*                 step. This function can be used to get all subsets of an
*                 n-set.
*
*INPUTS: code An nX1 or 1Xn vector consisting of zeros and ones
*             representing the current gray code value. The first code in
*             the sequence is all zeros. If getNexGrayCode is called with
*             an empty matrix for code and nCard=n, then the returned code
*             will be the first in the sequence.
*       nCard If code is empty, then this is the dimensionality of the code
*             sequence desired. Otherwise, this is the number of ones in
*             code. When getting a next code, nCard can speed things up,
*             but if omitted, it is just found as sum(code).
*
*OUTPUTS: code The next length codeLen gray code value in the sequence. If
*              the final gray code in the sequence was passed, then an
*              empty matrix is returned.
*        nCard The number of ones in the returned code, or n if the last
*              code was passed.
*       isLast True if the returned code is the last in the series and
*              passing it to getNexGrayCode would return an empty matrix.
*            j The index of the entry in code that was changed by this
*              function call. If this is the first function call and code
*              was just created, then j is an empty matrix.
*
*The algorithm is based on NEXSUB in Chapter 1 of [1]. Gray codes are
*also discussed in Chapter 7.2.1.1 of [2].
*
*The function called checks against the length of code so that if an
*invalid value of nCard is passed, it will not read/ write past the end of
*the array, even though it will return a code that is not the next in the
*sequence and value of nCard that has nothing to do with the code returned.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[code,nCard,isLast,j]=getNextGrayCode(code,nCard);
*
*REFERENCES:
*[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
*    and Calculators, 2nd ed. New York: Academic press, 1978.
*[2] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 2:
*    Generating all Tuples and Permutations, Upper Saddle River, NJ:
*    Addison-Wesley, 2009.
*
*October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "MexValidation.h"
/*This header is required by Matlab*/
#include "mex.h"
#include "combinatorialFuns.hpp"
//For the max command.
#include <algorithm>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t n=0;
    size_t j=0;
    size_t nCard;
    bool isLast=false;
    bool lastPassed=false;
    mxArray *codeArray=NULL;
    
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>2){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Too many outputs.");
    }
    
    //If an empty code matrix is passed, then the second argument is
    //required, and we will return the first gray code in the sequence.
    if(mxIsEmpty(prhs[0])) {
        mwSize dims[2];
        if(nrhs<2) {
            mexErrMsgTxt("The second argument is required when an empty code matrix is passed.");
        }
        
        dims[0]=getSizeTFromMatlab(prhs[1]);
        dims[1]=1;
        //Allocate the array; this also initializes all of the elements to
        //0.
        codeArray=mxCreateLogicalArray(2, dims);
        
        //This is the code
        plhs[0]=codeArray;
        
        if(nlhs>1) {
            //This is nCard
            plhs[1]=mxCreateDoubleScalar(0.0);
            
            if(nlhs>2) {
                //This is isLast
                plhs[2]=mxCreateLogicalScalar(false);
                
                if(nlhs>3) {
                    //This is j.
                    plhs[3]=mxCreateDoubleMatrix(0, 0, mxREAL);
                }
            }
        }
        return;
    }
    
    //The code array cannot be complex.
    if(mxIsComplex(prhs[0])!=false) {
        mexErrMsgTxt("The code array cannot be complex.");
    }
    
    //Copy the code value so that the input matrix is not modified on
    //return.
    {
        size_t n1,n2;
        n1=mxGetM(prhs[0]);
        n2=mxGetN(prhs[0]);

        if((n1==1&&n2>=1)||(n2==1&&n1>=1)) {
            n=std::max(n1,n2);
            
            codeArray=mxDuplicateArray(prhs[0]);
        } else {
            mexErrMsgTxt("The code vector has the wrong dimensionality.");
        }
    }

    if(nrhs>1) {
        nCard=getSizeTFromMatlab(prhs[1]);
    } else {
        mxArray *lhs[1];
        //Sum up the ones in codeArray to get nCard if it is not provided.
        mexCallMATLAB(1,lhs,1, &codeArray, "sum");
        nCard=getSizeTFromMatlab(lhs[0]);
        mxDestroyArray(lhs[0]);
    }
    
    //The type of the data in the code array is whatever the user passed to
    //the function. The function has to be called with the correct template 
    //value for the type of the code.
    lastPassed=false;
    switch(mxGetClassID(codeArray)){
        case mxCHAR_CLASS:
        {
            mxChar *code=mxGetChars(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            mxLogical* code=mxGetLogicals(codeArray); 
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double* code=mxGetDoubles(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float* code=mxGetSingles(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T* code=mxGetInt8s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T* code=mxGetUint8s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T* code=mxGetInt16s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T* code=mxGetUint16s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T* code=mxGetInt32s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T* code=mxGetUint32s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T* code=mxGetInt64s(codeArray);
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T* code=mxGetUint64s(codeArray);   
            
            if(nCard==static_cast<size_t>(code[n-1])&&nCard!=0) {
                lastPassed=true;
            } else{
                isLast=getNextGrayCodeCPP(n, code, nCard, j);
            }
            break;
        }
        case mxUNKNOWN_CLASS:
        case mxCELL_CLASS:
        case mxSTRUCT_CLASS:
        case mxVOID_CLASS:
        case mxFUNCTION_CLASS:
        case mxOPAQUE_CLASS:
        case mxOBJECT_CLASS:
        default:
            mexErrMsgTxt("The code vector is of an unknown data type.");
    }
    
    //If the final gray code was passed, then just return empty matrices
    //and put n in nCard. That way, if called again, the function will
    //start from the beginning.
    if(lastPassed==true) {
        mxDestroyArray(codeArray);
        plhs[0]=mxCreateDoubleMatrix(0, 0, mxREAL);
        if(nlhs>1) {
            mxArray *nCardMat=allocUnsignedSizeMatInMatlab(1, 1);
            if(sizeof(size_t)==4) {//32 bit
                *reinterpret_cast<size_t*>(mxGetUint32s(nCardMat))=n;
            } else {//64 bit
                *reinterpret_cast<size_t*>(mxGetUint64s(nCardMat))=n;
            }
 
            plhs[1]=nCardMat;
            if(nlhs>2) {
                //This is isLast
                plhs[2]=mxCreateDoubleMatrix(0, 0, mxREAL);;

                if(nlhs>3) {
                    //This is j
                    plhs[3]=mxCreateDoubleMatrix(0, 0, mxREAL);;
                }
            }
        }
        
        return;   
    }
    
    //Set the return values for the case when the last code was not passed.
    plhs[0]=codeArray;
    if(nlhs>1) {
        mxArray *nCardMat=allocUnsignedSizeMatInMatlab(1, 1);
        if(sizeof(size_t)==4) {//32 bit
            *reinterpret_cast<size_t*>(mxGetUint32s(nCardMat))=nCard;
        } else {//64 bit
            *reinterpret_cast<size_t*>(mxGetUint64s(nCardMat))=nCard;
        }

        plhs[1]=nCardMat;
        if(nlhs>2) {
            //This is isLast
            plhs[2]=mxCreateLogicalScalar(isLast);

            if(nlhs>3) {
                mxArray *jMat=allocUnsignedSizeMatInMatlab(1, 1);
                //Increment j to be an index for Matlab.
                j++;
                
                if(sizeof(size_t)==4) {//32 bit
                    *reinterpret_cast<size_t*>(mxGetUint32s(jMat))=j;
                } else {//64 bit
                    *reinterpret_cast<size_t*>(mxGetUint64s(jMat))=j;
                }
                
                plhs[3]=jMat;
            }
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
