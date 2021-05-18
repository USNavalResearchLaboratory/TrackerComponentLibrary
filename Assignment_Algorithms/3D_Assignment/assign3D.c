/**ASSIGN3D A C-code (for Matlab) implementation of a dual-primal
 *          Lagrangian relaxation algorithm for approximating the solution
 *          to the operations research axial 3D assignment problem. See the
 *          comments to the Matlab implementation for more details.
 *
 *It is assumed that all matrices are sufficiently small that it does not
 *matter whether lengths are held in size_t or ptrdiff_t variables.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * [tuples,fStar,qStar,u,exitCode]=assign3D(C,maximize,subgradMethod,subgradParams,maxIter,AbsTol,RelTol);
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//Prototypes for the actual assignment functions.
#include "assignAlgs3D.h"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    const size_t *nDims;
    double AbsTol=1e-10;
    double RelTol=0.05;
    double *C, *CCopy;
    size_t maxIter=20;
    bool maximize=false;
    int subgradMethod=0;
    //These are the parameters for the subgradient algorithm. Their
    //specific definitions vary depending on the subgradient algorithm that
    //has been selected.
    double param1, param2=0;
    size_t param3=0;
    mxArray *uMATLAB;
    ptrdiff_t *tuples;
    double *u;
    void *tempSpace;
    double fStar, qStar;
    ptrdiff_t exitCode;
    
    if(nrhs>7||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>5) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    ////////
    //Check the matrix C
    ///////
    if(mxIsComplex(prhs[0])) {
        mexErrMsgTxt("C must be real.");
        return;
    }
    
    if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("C must be of type double.");
        return;
    }

    if(mxIsEmpty(prhs[0])) {//The empty matrix special case.
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//tuples=[];
        
        if(nlhs>1) {
            plhs[1]=mxCreateDoubleScalar(0);//fStar=0
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleScalar(0);//qStar=0
            
                if(nlhs>3) {
                    plhs[3]=mxCreateDoubleMatrix(0,0,mxREAL);//u=[]
                    
                    if(nlhs>4) {
                        plhs[4]=mxCreateDoubleScalar(0);//exitCode=0
                    }
                }
            }
        }
        return;
    }

    C=mxGetDoubles(prhs[0]);
    
    //The scalar special case.
    if(mxIsScalar(prhs[0])) {
        mxArray *tuplesMATLAB=mxCreateDoubleMatrix(3,1,mxREAL);
        double *tuplesD=mxGetDoubles(tuplesMATLAB);
        tuplesD[0]=1;
        tuplesD[1]=1;
        tuplesD[2]=1;
        
        plhs[0]=tuplesMATLAB;
        if(nlhs>1) {
            plhs[1]=mxCreateDoubleScalar(C[0]);//fStar
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleScalar(C[0]);//qStar
                
                if(nlhs>3) {
                    plhs[3]=mxCreateDoubleScalar(0);//u=0
                    
                    if(nlhs>4) {
                        plhs[4]=mxCreateDoubleScalar(0);//exitCode=0;
                    }
                }
            }
        }
        
        return;
    }
    
    if(mxGetNumberOfDimensions(prhs[0])!=3) {
        mexErrMsgTxt("C must be 3D.");
        return;
    }

    nDims=mxGetDimensions(prhs[0]);

    if(!(nDims[0]<=nDims[1]&&nDims[1]<=nDims[2])) {
        mexErrMsgTxt("It is required that size(C,1)<=size(C,2)<=size(C,3)");
        return;
    }

    ////////
    //Check the other inputs
    ///////
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        maximize=getBoolFromMatlab(prhs[1]);
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        subgradMethod=getIntFromMatlab(prhs[2]);
        
        if(subgradMethod<0||subgradMethod>4) {
            mexErrMsgTxt("Invalid subgradMethod specified.");
            return;
        }
    }
    
    //Set the default parameters for the selected subgradient method.
    switch(subgradMethod) {
        case 0:
            param1=1;//gammaParam
            break;
        case 1:
            param1=2.8;//M
            param2=0.06;//r
            break;
        case 2:
            param1=0.3;//a
            param2=1.5;//b
            break;
        case 3:
        default://subgradMethod==4
            param1=2;//M=2
            param2=2.220446049250313e-16;//normBound=eps();
            param3=(size_t)nDims[2];//NR=n3;
    }

    //Get any provided subgradient parameters.
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        mxArray *curField;
        if(!mxIsStruct(prhs[3])) {
            mexErrMsgTxt("subgradMethod must be a structure.");
            return;
        }
        
        switch(subgradMethod) {
            case 0:
                curField=mxGetField(prhs[3],0,"gamma");
                if(curField!=NULL) {
                    param1=getDoubleFromMatlab(curField);
                }

                break;
            case 1:
                curField=mxGetField(prhs[3],0,"M");
                if(curField!=NULL) {
                    param1=getDoubleFromMatlab(curField);
                }
                
                curField=mxGetField(prhs[3],0,"r");
                if(curField!=NULL) {
                    param2=getDoubleFromMatlab(curField);
                }

                break;
            case 2:
                curField=mxGetField(prhs[3],0,"a");
                if(curField!=NULL) {
                    param1=getDoubleFromMatlab(curField);
                }
                
                curField=mxGetField(prhs[3],0,"b");
                if(curField!=NULL) {
                    param2=getDoubleFromMatlab(curField);
                }

                break;
            case 3:
            default:
                curField=mxGetField(prhs[3],0,"M");
                if(curField!=NULL) {
                    param1=getDoubleFromMatlab(curField);
                }
                
                curField=mxGetField(prhs[3],0,"normBound");
                if(curField!=NULL) {
                    param2=getDoubleFromMatlab(curField);
                }
                
                curField=mxGetField(prhs[3],0,"NR");
                if(curField!=NULL) {
                    param3=getSizeTFromMatlab(curField);
                }
        }
    }
        
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        maxIter=getSizeTFromMatlab(prhs[4]);
        
        if(maxIter<1) {
            mexErrMsgTxt("maxIter must be >=1.");
            return;
        }
    }
    
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        AbsTol=getDoubleFromMatlab(prhs[5]);
        
        if(AbsTol<0) {
            mexErrMsgTxt("AbsTol should be non-negative.");
            return;
        }
    }
    
    if(nrhs>6&&!mxIsEmpty(prhs[6])) {
        RelTol=getDoubleFromMatlab(prhs[6]);
        
        if(RelTol<0) {
            mexErrMsgTxt("RelTol should be non-negative.");
            return;
        }
    }
    
    ////////
    //Allocate space for temporary parameters and return variables.
    ///////
    //Duplicate the input array. If maximization is performed, the
    //original one will be overwritten.
    {
        const size_t CSize=(size_t)(nDims[0]*nDims[1]*nDims[2])*sizeof(double);
        CCopy=(double*)mxMalloc(CSize);
    
        memcpy(CCopy,C,CSize);
    }

    uMATLAB=mxCreateNumericMatrix((size_t)nDims[2],1,mxDOUBLE_CLASS,mxREAL);
    u=mxGetDoubles(uMATLAB);
    
    //Determine the size of the buffer needed for the function. We
    //allocate space for an additional 3*n1 ptrdiff_t values to hold the
    //tuples. The tuples returned by this function will be those converted
    //to doubles, because doubles are the default format in Matlab.
    tempSpace=mxMalloc(assign3DCBufferSize(nDims,subgradMethod)+3*(size_t)nDims[0]*sizeof(ptrdiff_t));
    tuples=(ptrdiff_t*)tempSpace;
    
    {
        void *bufferStart=(void*)(tuples+3*nDims[0]);
        exitCode=assign3DC(tuples,&fStar,&qStar,u,bufferStart,nDims,CCopy,maximize,subgradMethod,maxIter,AbsTol,RelTol,param1,param2,param3);
        mxFree(CCopy);
    }
    
    if(exitCode==-3||exitCode==-1) {
        // If no valid assignment was found.   
        mxFree(tempSpace);
        mxDestroyArray(uMATLAB);
        
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//tuples=[];
        
        if(nlhs>1) {
            plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);//fStar=[]'
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleMatrix(0,0,mxREAL);//qStar=[]
            
                if(nlhs>3) {
                    plhs[3]=mxCreateDoubleMatrix(0,0,mxREAL);//u=[]
                    
                    if(nlhs>4) {
                        plhs[4]=mxCreateDoubleScalar(exitCode);//exitCode
                    }
                }
            }
        }
        return;
    } else {
        const size_t numEls=3*nDims[0];
        size_t i;
        
        /* Convert the C indices into indices for MATLAB.*/
        for(i=0;i<numEls;i++) {
            tuples[i]++;
        }
        
        plhs[0]=ptrDiffTMat2MatlabDoubles(tuples,3,(size_t)nDims[0]);
        mxFree(tempSpace);
        if(nlhs>1) {
            plhs[1]=mxCreateDoubleScalar(fStar);
            
            if(nlhs>2) {
                plhs[2]=mxCreateDoubleScalar(qStar);
                        
                if(nlhs>3) {
                    plhs[3]=uMATLAB;
                    
                    if(nlhs>4) {
                        plhs[4]=mxCreateDoubleScalar(exitCode);
                    }
                } else {
                    mxDestroyArray(uMATLAB);
                }
            } else {
                mxDestroyArray(uMATLAB);
            }
        } else {
            mxDestroyArray(uMATLAB);
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
