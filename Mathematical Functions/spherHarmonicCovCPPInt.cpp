/**SPHERHARMONICCOVCPPINT A mex file interface to the C++ implementation of
 *                     the algorithm for computing variances and
 *                     covariances associated with estimates derived by
 *                     synthesizing spherical harmonic coefficients.
 *                     Generally, the Matlab function spherHarmonicEval
 *                     should be called instead of this one, as this
 *                     function does no input checking and running the
 *                     function with invalid inputs will crash Matlab.
 *
 *This function exists, because one can not access the properties of a
 *Matlab class in C++ without Matlab making a copy of the data. However, 
 *when large coefficient sets are used, this wastes memory and is slow.
 *Thus, rather than creating a mex file with the name spherHarmonicCov
 *that will get run in place of the .m file, a condition was inserted
 *into spherHarmonicCov.m to check for the existence of a compiled version
 *of this function and the individual elements of the ClusterSet classes
 *for the coefficients are passed.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The function is called in Matlab using the format:
 *[sigma2,Sigma,didOverflow]=spherHarmonicCovCPPInt(CStdDev,SStdDev,point,a,c,scalFactor);
 *or using 
 *[sigma2]=spherHarmonicCovCPPInt(CStdDev,SStdDev,point,a,c,scalFactor);
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
    double *point;
    ClusterSetCPP<double> CStdDev;
    ClusterSetCPP<double> SStdDev;
    size_t numPoints;
    size_t i, M, totalNumEls;
    mxArray *sigma2MATLAB;
    //This variable is only used if nlhs>1. It is set to zero here to
    //suppress a warning if compiled using -Wconditional-uninitialized.
    mxArray *SigmaMATLAB=NULL;
    double *sigma2,*Sigma;
    bool didOverflow;
    size_t *theOffsetArray, *theClusterSizes;
    
    if(nrhs!=6) {
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
        mexErrMsgTxt("Invalid ClusterSet data passed.");
    }

    CStdDev.clusterEls=reinterpret_cast<double*>(mxGetData(prhs[0]));
    SStdDev.clusterEls=reinterpret_cast<double*>(mxGetData(prhs[1]));

    //The number of elements in the ClusterSets
    totalNumEls=mxGetM(prhs[0])*mxGetN(prhs[0]);
    
    //We do not want to trust that the offset arrays and the cluster arrays
    //are valid. Also, the data types of the passed values can vary (e.g.,
    //might be doubles). Thus, we this function does not take the offset
    //and number of clusters data from a ClusterSet class. Rather, it
    //recreates the values based on the length of the data. This avoids
    //possibly reading past the end of an array.
    
    //Since the total number of points in a ClusterSet
    //is (M+1)*(M+2)/2, where M is the number of clusters -1, we can easily
    //verify that a valid number of points was passed.
    //Some versions of Windows SDK have problems with passing an integer
    //type to the sqrt function, so this explicit typecasting should take
    //care of the problem.
    M=(-3+static_cast<size_t>(sqrt(static_cast<double>(1+8*totalNumEls))))/2;        
    if((M+1)*(M+2)/2!=totalNumEls) {
        mexErrMsgTxt("The ClusterSets contain an inconsistent number of elements.");
    }
    
    //Otherwise, fill in the values for the cluster sizes and the ofset
    //array.
    CStdDev.numClust=M+1;
    SStdDev.numClust=M+1;
    
    theOffsetArray=new size_t[M+1];
    theClusterSizes=new size_t[M+1];
    theClusterSizes[0]=1;
    theOffsetArray[0]=0;
    for(i=1;i<=M;i++) {
       theOffsetArray[i]=theOffsetArray[i-1]+theClusterSizes[i-1];
       theClusterSizes[i]=i+1; 
    }

    CStdDev.totalNumEl=totalNumEls;
    SStdDev.totalNumEl=totalNumEls;
    CStdDev.offsetArray=theOffsetArray;
    SStdDev.offsetArray=theOffsetArray;
    CStdDev.clusterSizes=theClusterSizes;
    SStdDev.clusterSizes=theClusterSizes;

    //Get the other parameters.
    checkRealDoubleArray(prhs[2]);
    point=reinterpret_cast<double*>(mxGetData(prhs[2]));
    numPoints=mxGetN(prhs[2]);
    a=getDoubleFromMatlab(prhs[3]);
    c=getDoubleFromMatlab(prhs[4]);
    scalFactor=getDoubleFromMatlab(prhs[5]);
    
    //Allocate space for the return values
    sigma2MATLAB=mxCreateDoubleMatrix(numPoints, 1,mxREAL);
    sigma2=reinterpret_cast<double*>(mxGetData(sigma2MATLAB));
    
    if(nlhs>1) {
        mwSize dims[3];
        //A 3X3X numPoints array is needed for the covariance matrices.
        dims[0]=3;
        dims[1]=3;
        dims[2]=numPoints;
        
        SigmaMATLAB=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        Sigma=reinterpret_cast<double*>(mxGetData(SigmaMATLAB));
    } else {
        Sigma=NULL;
    }
    didOverflow=spherHarmonicCovCPP(sigma2,Sigma,CStdDev,SStdDev,point,numPoints,a,c,scalFactor);
            
    //Free the allocated space for the offet array and cluster sizes
    delete[] theOffsetArray;
    delete[] theClusterSizes;
            
    plhs[0]=sigma2MATLAB;
    
    if(nlhs>1) {
        plhs[1]=SigmaMATLAB;
        if(nlhs>2) {
            plhs[2]=boolMat2Matlab(&didOverflow,1,1);
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
