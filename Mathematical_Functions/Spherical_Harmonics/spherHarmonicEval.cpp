/**SPHERHARMONICEVAL A mex file implementation of the function
 *                  spherHarmonicEval. See the comments to the Matlab
 *                  implementation for more details.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The function is called in Matlab using the format:
 *[V,gradV,HessianV]=spherHarmonicEval(C,S,point,a,c,systemType,spherDerivs,scalFactor,algorithm);
 *or using 
 *[V]=spherHarmonicEval(C,S,point,a,c,systemType,spherDerivs,scalFactor,algorithm);
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

using namespace std;

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double scalFactor;
    //Copies are only made if real and complex inputs are mixed. In such an
    //instance, we convert everything to complex. If copies are made, they
    //will need to be freed prior to returning.
    mxArray *CCopy=NULL;
    mxArray *SCopy=NULL;
    std::complex <double> a,c;
    bool spherDerivs;
    size_t algorithm;
    double *point, *pointCopy=NULL;
    size_t numPoints, pointDim;
    size_t M, totalNumEls;
    mxArray *VMATLAB;
    mxArray *gradVMATLAB;
    mxArray *HessianVMATLAB;
    bool useComplexAlg=false;
    mxComplexity complexVal;
    size_t systemType;
    
    if(nrhs<3||nrhs>9) {
        mexErrMsgTxt("Wrong number of inputs.");
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Invalid number of outputs.");
    }

    //Check the validity of the coefficients.
    checkDoubleArray(prhs[0]);
    checkDoubleArray(prhs[1]);
    
    //Make sure that both of the coefficient sets have the same number of
    //elements and are not empty.
    if(mxIsEmpty(prhs[0])||mxIsEmpty(prhs[1])||mxGetM(prhs[0])!=mxGetM(prhs[1])||mxGetN(prhs[0])!=mxGetN(prhs[1])) {
        mexErrMsgTxt("Invalid data passed.");
    }
    
    //The number of elements in the ClusterSets
    totalNumEls=mxGetM(prhs[0])*mxGetN(prhs[0]);

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
        
    //Get the other parameters before getting C and S.
    pointDim=mxGetM(prhs[2]);
    if((pointDim!=3&&pointDim!=2)||mxIsEmpty(prhs[2])) {
        mexErrMsgTxt("The points have an incorrect dimensionality.");
    }

    checkRealDoubleArray(prhs[2]);
    numPoints=mxGetN(prhs[2]);
    point=mxGetDoubles(prhs[2]);
    if(pointDim==2) {
        size_t curPoint;
        double *curCopyPoint;
        double *curOrigPoint;
        //Add a range component if one was not provided.
        pointCopy=new double[3*numPoints];

        curCopyPoint=pointCopy;
        curOrigPoint=point;
        for(curPoint=0;curPoint<numPoints;curPoint++) {
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
            a=1.0;
        } else {
            a=getScalarMatlabClassConst("Constants", "EGM2008SemiMajorAxis");
        }
    } else {
        if(mxIsComplex(prhs[3])) {
            a=getComplexDoubleFromMatlab(prhs[3]);
            useComplexAlg=true;
        } else {
            a=static_cast<double>(getDoubleFromMatlab(prhs[3]));
        }
    }
    
    if(nrhs<5||mxIsEmpty(prhs[4])) {
        if(pointDim==2) {
            c=1;
        } else {
            c=getScalarMatlabClassConst("Constants", "EGM2008GM");
        }
    } else {
        if(mxIsComplex(prhs[4])) {
            c=getComplexDoubleFromMatlab(prhs[4]);
            useComplexAlg=true;
        } else {
            c=getDoubleFromMatlab(prhs[4]);
        }
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
    
    useComplexAlg=useComplexAlg||mxIsComplex(prhs[0])||mxIsComplex(prhs[1]);
    
    if(useComplexAlg) {
        complexVal=mxCOMPLEX;
    } else {
        complexVal=mxREAL;
    }
    
    //Allocate space for the return values.
    VMATLAB=mxCreateDoubleMatrix(numPoints,1,complexVal);
    if(nlhs>1) {
        gradVMATLAB=mxCreateDoubleMatrix(3,numPoints,complexVal);
        
        if(nlhs>2) {
            mwSize dims[3];

            dims[0]=3;
            dims[1]=3;
            dims[2]=numPoints;

            HessianVMATLAB=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,complexVal);
        }
    }    
    
    if(useComplexAlg) {
        complex<double> *V;
        complex<double> *gradV=NULL;
        complex<double> *HessianV=NULL;
        CountingClusterSetCPP<complex<double>> C;
        CountingClusterSetCPP<complex<double>> S;
        
        C.numClust=M+1;
        S.numClust=M+1;
        C.totalNumEl=totalNumEls;
        S.totalNumEl=totalNumEls;
        
        if(!mxIsComplex(prhs[0])) {
            CCopy=mxDuplicateArray(prhs[0]);
            mxMakeArrayComplex(CCopy);
            C.clusterEls=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(CCopy));
        } else {
            C.clusterEls=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(prhs[0]));
        }
        
        if(!mxIsComplex(prhs[1])) {
            SCopy=mxDuplicateArray(prhs[1]);
            mxMakeArrayComplex(SCopy);
            S.clusterEls=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(SCopy));
        } else {
            S.clusterEls=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(prhs[1]));
        }
        
        V=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(VMATLAB));
        if(nlhs>1) {
            gradV=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(gradVMATLAB));
            if(nlhs>2) {
                HessianV=reinterpret_cast<complex<double>*>(mxGetComplexDoubles(HessianVMATLAB));
            }
        }

        spherHarmonicEvalCPPComplex(V,gradV,HessianV,C,S,point,numPoints,a,c,systemType,spherDerivs,scalFactor,algorithm);
    } else {//All data is real.
        double *V;
        double *gradV=NULL;
        double *HessianV=NULL;
        CountingClusterSetCPP<double> C;
        CountingClusterSetCPP<double> S;
        
        C.numClust=M+1;
        S.numClust=M+1;
        C.totalNumEl=totalNumEls;
        S.totalNumEl=totalNumEls;
        C.clusterEls=mxGetDoubles(prhs[0]);
        S.clusterEls=mxGetDoubles(prhs[1]);
        
        V=mxGetDoubles(VMATLAB);
        if(nlhs>1) {
            gradV=mxGetDoubles(gradVMATLAB);
            if(nlhs>2) {
                HessianV=mxGetDoubles(HessianVMATLAB);
            }
        }
  
       spherHarmonicEvalCPPReal(V,gradV,HessianV,C,S,point,numPoints,real(a),real(c),systemType,spherDerivs,scalFactor,algorithm);  
    }

    plhs[0]=VMATLAB;
    if(nlhs>1) {
        plhs[1]=gradVMATLAB;
        if(nlhs>2) {
            plhs[2]=HessianVMATLAB;
        }
    }

    if(pointCopy!=NULL) {
        delete pointCopy;
    }
    
    if(CCopy!=NULL) {
        mxDestroyArray(CCopy);
    }
    
    if(SCopy!=NULL) {
        mxDestroyArray(SCopy);
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
