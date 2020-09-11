/*NALEGENDRECOSRAT A C++ implementation to evaluate
 *                \bar{P}_{nm}(cos(theta))/u^m for all n from 0 to M and
 *                for each n for all m from 0 to n, where
 *                u=abs(sin(theta)). \bar{P}_{nm}(x) is the fully
 *                normalized associated Legendre function of x of degree n
 *                and order m. Also evaluate
 *                D{\bar{P}_{nm}(cos(theta))}/u^m and
 *                D2{\bar{P}_{nm}(cos(theta))}/u^m where D{} is the first
 *                derivative operator with respect to theta and D2{} is the
 *                second derivative operator with respect to theta. All of
 *                the values can be scaled by a factor of scalFac, if
 *                desired, to help prevent overflows with high degrees and
 *                orders.
 *
 *INPUTS: theta An angle in radians.
 *      maxDeg The maximum degree and order of the output. This must be
 *             >=0.
 *  scalFactor A scale factor to help prevent overflow of the results. In
 *             [1], discussed below, a value of 10^(-280) is used.
 *
 *OUTPUTS: PBarUVals An instance of the CountingClusterSet class such that
 *                  PBarUVals(n+1,m+1)=scalFac*\bar{P}_{nm}(cos(theta))/u^m.
 * dPBarUValsdTheta An instance of the CountingClusterSet class such that
 *                  dPBarUValsdTheta(n+1,m+1)=scalFac*D{\bar{P}_{nm}(cos(theta))}/u^m
 * d2PBarUValsdTheta2 An instance of the CountingClusterSet class such that
 *                  d2PBarUValsdTheta2(n+1,m+1)=scalFac*D2{\bar{P}_{nm}(cos(theta))}/u^m
 *
 *Comments are given in the Matlab implementation. This
 *implementation tends to be approximately 11,000 times faster than the
 *Matlab implementation.
 *
 *Note that this function requires a large amount of memory, since Matlab's
 *mxSetProperty function copies the data rather than allowing pointers to
 *be passed. Thus, more than double the normal amount of space is required
 *for the return values.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[PBarUVals,dPBarUValsdTheta,d2PBarUValsdTheta2]=NALegendreCosRat(theta,M,scalFactor);
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double theta, scalFactor;
    CountingClusterSetCPP<double> PBarUVals;
    CountingClusterSetCPP<double> dPBarUValsdTheta;//The first derivatives
    CountingClusterSetCPP<double> d2PBarUValsdTheta2;//The second derivatives
    size_t M, numPBarU;
    mxArray *CSRetVal;
    mxArray *clusterElsMATLAB;
    mxArray *clusterEls1stDerivMATLAB, *numClustMATLAB;
    
    if(nrhs!=3){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return; 
    }

    theta=getDoubleFromMatlab(prhs[0]);
    //M>=0. This command will throw an error if the value is negative.
    M=getSizeTFromMatlab(prhs[1]);
    scalFactor=getDoubleFromMatlab(prhs[2]);
    
    numPBarU=(M+1)*(M+2)/2;
    
    //Allocate space for the results.
    clusterElsMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
    {
        double temp=static_cast<double>(M+1);
        numClustMATLAB=doubleMat2Matlab(&temp,1,1);
    }

    PBarUVals.numClust=M+1;
    PBarUVals.totalNumEl=numPBarU;
    PBarUVals.clusterEls=mxGetDoubles(clusterElsMATLAB);

    NALegendreCosRatCPP(PBarUVals,theta,scalFactor);
    
    //Set the first return value
    mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
    mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
    mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);
    
    plhs[0]=CSRetVal;

    if(nlhs>1) {//Compute the first derivatives, if they are desired.
        clusterEls1stDerivMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
        
        dPBarUValsdTheta.numClust=M+1;
        dPBarUValsdTheta.totalNumEl=numPBarU;
        dPBarUValsdTheta.clusterEls=mxGetDoubles(clusterEls1stDerivMATLAB);

        NALegendreCosRatDerivCPP(dPBarUValsdTheta, PBarUVals, theta);
        
        //Set the second return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls1stDerivMATLAB);
        mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);

        plhs[1]=CSRetVal;
    }
    
    if(nlhs>2) {//Compute the second derivatives if they are desired.
        mxArray *clusterEls2ndDerivMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
        
        d2PBarUValsdTheta2.numClust=M+1;
        d2PBarUValsdTheta2.totalNumEl=numPBarU;
        d2PBarUValsdTheta2.clusterEls=mxGetDoubles(clusterEls2ndDerivMATLAB);
        
        NALegendreCosRatDeriv2CPP(d2PBarUValsdTheta2, dPBarUValsdTheta, PBarUVals,theta);
        
        //Set the third return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls2ndDerivMATLAB);
        mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);

        plhs[2]=CSRetVal;
        mxDestroyArray(clusterEls2ndDerivMATLAB);
        mxDestroyArray(clusterEls1stDerivMATLAB);
    } else if(nlhs>1) {
        mxDestroyArray(clusterEls1stDerivMATLAB);
    }
    
    //Free the buffers. The mxSetProperty command copied the data.
    mxDestroyArray(clusterElsMATLAB);
    mxDestroyArray(numClustMATLAB);
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
