/**SPHERHARMONICEVAL A mex file implementation of the function
 *                   spherHarmonicEval. See the comments to the Matlab
 *                   implementation for more details.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The function is called in Matlab using the format:
 *[V,gradV,HessianV]=spherHarmonicSetEval(C,S,point,a,c,systemType,spherDerivs,scalFactor,algorithm);
 *or using 
 *[V]=spherHarmonicSetEval(C,S,point,a,c,systemType,spherDerivs,scalFactor,algorithm);
 *if one only wants the potential. The function executes faster if only the
 *potential and not the gradient and Hessian need be computed.
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathFuncs.hpp"
//Needed for sqrt
#include <cmath>
//For the case that complex values are passed.
#include <complex>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t numSets;
    double scalFactor;
    std::complex <double> a,c;
    bool spherDerivs;
    size_t algorithm;
    double *point, *pointCopy=NULL;
    CountingClusterSetVecCPP<double> CReal;
    CountingClusterSetVecCPP<double> SReal;
    CountingClusterSetVecCPP<double> CImag;
    CountingClusterSetVecCPP<double> SImag;
    size_t numPoints, pointDim;
    size_t i, M, totalNumEls;
    mxArray *VMATLAB;
    mxArray *gradVMATLAB;
    mxArray *HessianVMATLAB;
    double *VReal;
    double *VImag;
    double *gradVReal=NULL;
    double *gradVImag=NULL;
    double *HessianVReal=NULL;
    double *HessianVImag=NULL;
    bool useComplexAlg=false;
    mxComplexity complexVal;
    size_t systemType;
    //Buffer is only used to hold temporary space in case the complex
    //algorithm is run but not all (real/ complex) parts are available. For
    //example, if S is real and C is complex.
    double *buffer=NULL;
    
    if(nrhs<3||nrhs>9) {
        mexErrMsgTxt("Wrong number of inputs.");
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Invalid number of outputs.");
    }

    //Check the validity of the ClusterSets
    checkDoubleArray(prhs[0]);
    checkDoubleArray(prhs[1]);
    numSets=mxGetN(prhs[0]);
    
    //The number of elements in the ClusterSets
    totalNumEls=mxGetM(prhs[0]);
    
    //Make sure that both of the coefficient sets have the same number of
    //elements and are not empty.
    if(mxIsEmpty(prhs[0])||mxIsEmpty(prhs[1])||totalNumEls!=mxGetM(prhs[1])||numSets!=mxGetN(prhs[1])) {
        mexErrMsgTxt("Invalid data passed.");
    }
    
    CReal.clusterEls=mxGetPr(prhs[0]);
    CImag.clusterEls=mxGetPi(prhs[0]);
    useComplexAlg=useComplexAlg|(CImag.clusterEls!=NULL);
    
    SReal.clusterEls=mxGetPr(prhs[1]);
    SImag.clusterEls=mxGetPi(prhs[1]);
    useComplexAlg=useComplexAlg|(SImag.clusterEls!=NULL);

    //Since the total number of points in a CountingClusterSetCPP
    //is (M+1)*(M+2)/2, where M is the number of clusters -1, we can easily
    //verify that a valid number of points was passed.
    //Some versions of Windows SDK have problems with passing an integer
    //type to the sqrt function, so this explicit typecasting should take
    //care of the problem.
    M=(-3+static_cast<size_t>(sqrt(static_cast<double>(1+8*totalNumEls))))/2;        
    if((M+1)*(M+2)/2!=totalNumEls) {
        mexErrMsgTxt("S and C contain an inconsistent number of elements.");
    }
    
    //Otherwise, fill in the values for the cluster sizes and the ofset
    //array.
    CReal.numClust=M+1;
    CImag.numClust=M+1;
    SReal.numClust=M+1;
    SImag.numClust=M+1;
    CReal.totalNumEl=totalNumEls;
    CImag.totalNumEl=totalNumEls;
    SReal.totalNumEl=totalNumEls;
    SImag.totalNumEl=totalNumEls;
    CReal.numSets=numSets;
    CImag.numSets=numSets;
    SReal.numSets=numSets;
    SImag.numSets=numSets;

    //Get the other parameters.
    pointDim=mxGetM(prhs[2]);
    if((pointDim!=3&&pointDim!=2)||mxIsEmpty(prhs[2])) {
        mexErrMsgTxt("The points have an incorrect dimensionality.");
    }

    checkRealDoubleArray(prhs[2]);
    numPoints=mxGetN(prhs[2]);
    point=reinterpret_cast<double*>(mxGetData(prhs[2]));
    if(pointDim==2) {
        size_t curPoint;
        double *curCopyPoint;
        double *curOrigPoint;
        //Add a range component if one was not provided.
        pointCopy=new double[3*numPoints];

        curCopyPoint=pointCopy;
        curOrigPoint=point;
        for(curPoint=1;curPoint<numPoints;curPoint++) {
            *(curCopyPoint)=1;
            *(curCopyPoint+1)=*(curOrigPoint);
            *(curCopyPoint+2)=*(curOrigPoint+1);

            curCopyPoint+=3;
            curOrigPoint+=2;
        }

        point=pointCopy;
    }

    if(nrhs<4||mxIsEmpty(prhs[3])) {
        if(pointDim==2) {
            a=1;
        } else {
            a=getScalarMatlabClassConst("Constants", "EGM2008SemiMajorAxis");
        }
    } else {
        a=getComplexDoubleFromMatlab(prhs[3]);
        useComplexAlg=useComplexAlg|(imag(a)!=0);
    }
    
    if(nrhs<5||mxIsEmpty(prhs[4])) {
        if(pointDim==2) {
            c=1;
        } else {
            c=getScalarMatlabClassConst("Constants", "EGM2008GM");
        }
    } else {
        c=getComplexDoubleFromMatlab(prhs[4]);
        useComplexAlg=useComplexAlg|(imag(c)!=0);
    }
    
    if(nrhs<6||mxIsEmpty(prhs[5])) {
        systemType=0;
    } else {
        systemType=getSizeTFromMatlab(prhs[5]);
    }

    if(systemType!=0&&systemType!=2) {
        mexErrMsgTxt("An unsupported systemType was specified.");
    }
    
    if(nrhs<7||mxIsEmpty(prhs[6])) {
        spherDerivs=false;
    } else {
        spherDerivs=getBoolFromMatlab(prhs[6]);
    }
    
    if(nrhs<8||mxIsEmpty(prhs[7])) {
        scalFactor=1e-280;
    } else {
        scalFactor=getDoubleFromMatlab(prhs[7]);
    }
    
    if(nrhs<9||mxIsEmpty(prhs[8])) {
        algorithm=0;
    } else {
        algorithm=getSizeTFromMatlab(prhs[8]);
    }
    
    if(algorithm>2) {
        mexErrMsgTxt("Unknown algorithm option specified.");
    }
    
    //Allocate space for the return values
    if(useComplexAlg) {
        size_t numNulls;
        double *buffCur;
        complexVal=mxCOMPLEX;
        
        //If the real or complex component of C or S is missing, then
        //allocate a matrix ful of zeros.
        numNulls=(CReal.clusterEls==NULL)+(CImag.clusterEls==NULL)+(SReal.clusterEls==NULL)+(SImag.clusterEls==NULL);
        buffer=reinterpret_cast<double*>(mxCalloc(numNulls*totalNumEls*numSets,sizeof(double)));
        
        buffCur=buffer;
        if(CReal.clusterEls==NULL) {
            CReal.clusterEls=buffCur;
            buffCur+=totalNumEls*numSets;
        }
        
        if(CImag.clusterEls==NULL) {
            CImag.clusterEls=buffCur;
            buffCur+=totalNumEls*numSets;
        }
        
        if(SReal.clusterEls==NULL) {
            SReal.clusterEls=buffCur;
            buffCur+=totalNumEls*numSets;
        }
        
        if(SImag.clusterEls==NULL) {
            SImag.clusterEls=buffCur;
        }
    } else {  
        complexVal=mxREAL;
    }
    
    VMATLAB=mxCreateDoubleMatrix(numSets,numPoints,complexVal);
    VReal=mxGetPr(VMATLAB);
    VImag=mxGetPi(VMATLAB);//If it is purely real, this is just NULL.

    if(nlhs>1) {
        mwSize dimsV[3];
        dimsV[0]=3;
        dimsV[1]=numSets;
        dimsV[2]=numPoints;
        
        gradVMATLAB=mxCreateNumericArray(3,dimsV,mxDOUBLE_CLASS,complexVal);
        gradVReal=mxGetPr(gradVMATLAB);
        gradVImag=mxGetPi(gradVMATLAB);

        if(nlhs>2) {
            mwSize dims[4];
            dims[0]=3;
            dims[1]=3;
            dims[2]=numSets;
            dims[3]=numPoints;

            HessianVMATLAB=mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,complexVal);
            HessianVReal=mxGetPr(HessianVMATLAB);
            HessianVImag=mxGetPi(HessianVMATLAB);
        }
    }

    if(useComplexAlg) {
        spherHarmonicSetEvalCPPComplex(VReal, VImag, gradVReal, gradVImag, HessianVReal, HessianVImag, CReal, CImag, SReal, SImag, point, numPoints, a, c, systemType, spherDerivs,scalFactor,algorithm);
    } else {
        spherHarmonicSetEvalCPPReal(VReal,gradVReal,HessianVReal,CReal,SReal,point,numPoints,real(a),real(c),systemType, spherDerivs,scalFactor,algorithm);
    }

    plhs[0]=VMATLAB;
    if(nlhs>1) {
        plhs[1]=gradVMATLAB;
        if(nlhs>2) {
            plhs[2]=HessianVMATLAB;
        }
    }
    
    if(pointCopy!=NULL) {
        mxFree(pointCopy);
    }
    
    if(buffer!=NULL) {
        mxFree(buffer);
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
