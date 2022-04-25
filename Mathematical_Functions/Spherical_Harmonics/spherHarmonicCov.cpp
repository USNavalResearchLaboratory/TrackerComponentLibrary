/**SPHERHARMONICCOV A mex file implementation of the function
 *                  spherHarmonicCov. See the comments to the Matlab
 *                  implementation for more details.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The function is called in Matlab using the format:
 *[sigma2,Sigma,didOverflow]=spherHarmonicCov(CStdDev,SStdDev,point,a,c,scalFactor);
 *or using 
 *[sigma2]=spherHarmonicCov(CStdDev,SStdDev,point,a,c,scalFactor);
 *if one only the variance of the potential is desired. The function
 *executes faster if only the variance of the potential and not the
 *covariance matrix of the gradient need be computed.
 *
 *April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include"matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathFuncs.hpp"
//For sqrt
#include <cmath>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double a,c,scalFactor;
    double *point, *pointCopy=NULL;
    CountingClusterSetCPP<double> CStdDev;
    CountingClusterSetCPP<double> SStdDev;
    size_t numPoints, pointDim;
    size_t M, totalNumEls;
    mxArray *sigma2MATLAB;
    //This variable is only used if nlhs>1. It is set to zero here to
    //suppress a warning if compiled using -Wconditional-uninitialized.
    mxArray *SigmaMATLAB=NULL;
    double *sigma2,*Sigma;
    bool didOverflow;
    
    if(nrhs<3||nrhs>6) {
        mexErrMsgTxt("Wrong number of inputs.");
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Invalid number of outputs.");
    }
    
    //Check the validity of the ClusterSets
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    //Make sure that both of the coefficient sets have the same number of
    //elements and are not empty.
    if(mxIsEmpty(prhs[0])||mxIsEmpty(prhs[1])||mxGetM(prhs[0])!=mxGetM(prhs[1])||mxGetN(prhs[0])!=mxGetN(prhs[1])) {
        mexErrMsgTxt("Invalid data passed.");
    }

    CStdDev.clusterEls=mxGetDoubles(prhs[0]);
    SStdDev.clusterEls=mxGetDoubles(prhs[1]);

    //The number of elements in the CountingClusterSetCPP
    totalNumEls=mxGetM(prhs[0])*mxGetN(prhs[0]);
    
    //Since the total number of points in a CountingClusterSetCPP
    //is (M+1)*(M+2)/2, where M is the number of clusters -1, we can easily
    //verify that a valid number of points was passed.
    //Some versions of Windows SDK have problems with passing an integer
    //type to the sqrt function, so this explicit typecasting should take
    //care of the problem.
    M=(-3+static_cast<size_t>(sqrt(static_cast<double>(1+8*totalNumEls))))/2;
    
    if(M<3) {
        mexErrMsgTxt("The coefficients must be provided to at least degree 3. To use a lower degree, one can insert zero coefficients.");
    }
    
    if((M+1)*(M+2)/2!=totalNumEls) {
        mexErrMsgTxt("The set of harmonic coefficients contain an inconsistent number of elements.");
    }
    
    //Otherwise, fill in the values for the cluster sizes and the ofset
    //array.
    CStdDev.numClust=M+1;
    SStdDev.numClust=M+1;

    CStdDev.totalNumEl=totalNumEls;
    SStdDev.totalNumEl=totalNumEls;
    
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
            a=1;
        } else {
            a=getScalarMatlabClassConst("Constants", "EGM2008SemiMajorAxis");
        }
    } else {
        a=getDoubleFromMatlab(prhs[3]);
    }
    
    if(nrhs<5||mxIsEmpty(prhs[4])) {
        if(pointDim==2) {
            c=1;
        } else {
            c=getScalarMatlabClassConst("Constants", "EGM2008GM");
        }
    } else {
        c=getDoubleFromMatlab(prhs[4]);
    }
    
    if(nrhs<6||mxIsEmpty(prhs[5])) {
        scalFactor=3.280212943147993e-142;
    } else {
        scalFactor=getDoubleFromMatlab(prhs[5]);
    }
    
    //Allocate space for the return values
    sigma2MATLAB=mxCreateDoubleMatrix(numPoints, 1,mxREAL);
    sigma2=mxGetDoubles(sigma2MATLAB);
    
    if(nlhs>1) {
        mwSize dims[3];
        //A 3X3X numPoints array is needed for the covariance matrices.
        dims[0]=3;
        dims[1]=3;
        dims[2]=numPoints;
        
        SigmaMATLAB=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        Sigma=mxGetDoubles(SigmaMATLAB);
    } else {
        Sigma=NULL;
    }
    didOverflow=spherHarmonicCovCPP(sigma2,Sigma,CStdDev,SStdDev,point,numPoints,a,c,scalFactor);

    plhs[0]=sigma2MATLAB;
    
    if(nlhs>1) {
        plhs[1]=SigmaMATLAB;
        if(nlhs>2) {
            plhs[2]=boolMat2Matlab(&didOverflow,1,1);
        }
    }
    
    if(pointCopy!=NULL) {
        mxFree(pointCopy);
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
