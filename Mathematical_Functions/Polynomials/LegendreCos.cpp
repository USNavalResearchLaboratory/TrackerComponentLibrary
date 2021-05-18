/**LEGENDRECOS A C++ implementation to evaluate the fully associated
 *          Legendre functions of cos(theta) of degree n and order m for
 *          all n from 0 to M and for each n, m goes from 0 to n. Such a
 *          function is typically written as \bar{P}_{nm}(cos(theta)). Also
 *          evaluate D{\bar{P}_{nm}(cos(theta))} and
 *          D2{\bar{P}_{nm}(cos(theta))}, where D{} is the first derivative
 *          with respect to theta and D2{} is the second derivative with
 *          respect to theta. All of the values can be scaled by a factor
 *          of scalFac, if desired, to help prevent overflows with high
 *          degrees and orders. 
 *
 *INPUTS: theta An angle in radians.
 *      maxDeg The maximum degree and order of the output. This must be
 *             >=0.
 *  scalFactor A scale factor to help prevent overflow of the results. If
 *             this parameter is omitted or an empty matrix is passed, the
 *             default of 1 is used.
 *
 *OUTPUTS: PBarVals An instance of the ClusterSet class such that
 *                  PBarVals(n+1,m+1)=scalFac*\bar{P}_{nm}(cos(theta)).
 *                  To extract all coefficients as a vector just call
 *                  PBarVals(:).
 *  dPBarValsdTheta An instance of the ClusterSet class such that
 *                  dPBarValsdTheta(n+1,m+1)=scalFac*D{\bar{P}_{nm}(cos(theta))}
 * d2PBarValsdTheta2 An instance of the ClusterSet class such that
 *                  d2PBarValsdTheta2(n+1,m+1)=scalFac*D2{\bar{P}_{nm}(cos(theta))}
 *
 *Comments are given in the Matlab implementation.
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
 *[PBarVals,dPBarValsdTheta,d2PBarValsdTheta2]=NALegendreCos(theta,M,scalFactor);
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathFuncs.hpp"

#include <math.h>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double scalFactor,theta,u,um;
    CountingClusterSetCPP<double> PBarUVals;
    CountingClusterSetCPP<double> dPBarUValsdTheta;
    CountingClusterSetCPP<double> PBarVals;
    CountingClusterSetCPP<double> dPBarValsdTheta;//The first derivatives
    CountingClusterSetCPP<double> d2PBarValsdTheta2;//The second derivatives
    size_t M, numPBar,n,m,curEl;
    mxArray *CSRetVal;
    mxArray *clusterElsMATLAB, *numClustMATLAB;
    
    if(nrhs<2||nrhs>3){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return; 
    }

    theta=getDoubleFromMatlab(prhs[0]);
    M=getSizeTFromMatlab(prhs[1]);
    
    if(nrhs<3||mxIsEmpty(prhs[2])) {
        scalFactor=1;
    } else {
        scalFactor=getDoubleFromMatlab(prhs[2]);
    }
    
    numPBar=(M+1)*(M+2)/2;
    //Allocate space for the results.
    clusterElsMATLAB=mxCreateDoubleMatrix(numPBar,1,mxREAL);
    {
        double temp=static_cast<double>(M+1);
        numClustMATLAB=doubleMat2Matlab(&temp,1,1);
    }

    PBarUVals.numClust=M+1;
    PBarUVals.totalNumEl=numPBar;
    PBarUVals.clusterEls=reinterpret_cast<double*>(mxMalloc(sizeof(double)*numPBar));
    
    PBarVals.numClust=M+1;
    PBarVals.totalNumEl=numPBar;
    PBarVals.clusterEls=mxGetDoubles(clusterElsMATLAB);
    
    u=fabs(sin(theta));

    NALegendreCosRatCPP(PBarUVals,theta,scalFactor);

    PBarVals.clusterEls[0]=PBarUVals.clusterEls[0];
    curEl=1;
    for(n=1;n<=M;n++) {
        um=u;
        PBarVals.clusterEls[curEl]=PBarUVals.clusterEls[curEl];
        curEl++;
        for(m=1;m<=n;m++) {
            PBarVals.clusterEls[curEl]=PBarUVals.clusterEls[curEl]*um;
            curEl++;
            um*=u;
        }
    }
    
    //Set the first return value
    mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
    mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
    mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);
    mxDestroyArray(clusterElsMATLAB);
    
    plhs[0]=CSRetVal;
    
    if(nlhs>1) {//Compute the first derivatives, if they are desired.
        clusterElsMATLAB=mxCreateDoubleMatrix(numPBar,1,mxREAL);
        
        dPBarUValsdTheta.numClust=M+1;
        dPBarUValsdTheta.totalNumEl=numPBar;
        dPBarUValsdTheta.clusterEls=reinterpret_cast<double*>(mxMalloc(sizeof(double)*numPBar));
        
        dPBarValsdTheta.numClust=M+1;
        dPBarValsdTheta.totalNumEl=numPBar;
        dPBarValsdTheta.clusterEls=mxGetDoubles(clusterElsMATLAB);

        NALegendreCosRatDerivCPP(dPBarUValsdTheta, PBarUVals, theta);
        
        dPBarValsdTheta.clusterEls[0]=dPBarUValsdTheta.clusterEls[0];
        curEl=1;
        for(n=1;n<=M;n++) {
            um=u;
            dPBarValsdTheta.clusterEls[curEl]=dPBarUValsdTheta.clusterEls[curEl];
            curEl++;
            for(m=1;m<=n;m++) {
                dPBarValsdTheta.clusterEls[curEl]=dPBarUValsdTheta.clusterEls[curEl]*um;
                curEl++;
                um*=u;
            }
        }
        
        //Set the second return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
        mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);
        mxDestroyArray(clusterElsMATLAB);

        plhs[1]=CSRetVal;
    }
    
    if(nlhs>2) {//Compute the second derivatives if they are desired.
        clusterElsMATLAB=mxCreateDoubleMatrix(numPBar,1,mxREAL);
        
        d2PBarValsdTheta2.numClust=M+1;
        d2PBarValsdTheta2.totalNumEl=numPBar;
        d2PBarValsdTheta2.clusterEls=mxGetDoubles(clusterElsMATLAB);
        
        NALegendreCosRatDeriv2CPP(d2PBarValsdTheta2, dPBarUValsdTheta, PBarUVals,theta);
        
        curEl=1;
        for(n=1;n<=M;n++) {
            um=u;
            curEl++;
            for(m=1;m<=n;m++) {
                d2PBarValsdTheta2.clusterEls[curEl]*=um;
                curEl++;
                um*=u;
            }
        }

        //Set the third return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "CountingClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
        mxSetProperty(CSRetVal,0,"numClust",numClustMATLAB);
        mxDestroyArray(clusterElsMATLAB);
        
        plhs[2]=CSRetVal;
        mxFree(dPBarUValsdTheta.clusterEls);
        
    } else if(nlhs>1) {
        mxFree(dPBarUValsdTheta.clusterEls);
    }
    
    //Free the buffers. The mxSetProperty command copied the data.
    mxFree(PBarUVals.clusterEls);
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
