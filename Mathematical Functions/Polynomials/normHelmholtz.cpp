/**NORMHELMHOLTZ A C++ implementation of a function that computes the fully
*              normalized derived Legendre functions (also known as fully
*              normalized Helmholtz polynomials) of degree n=0 to M and
*              for each n from order m=0 to n evaluated at the point u.
*              That is, \bar{H}^m_n(u). This also finds the first and
*              second derivatives of the normalized Helmholtz polynomials
*              with respect to the parameter u, D{\bar{H}^m_n(u)} and
*              D2{\bar{H}^m_n(u)}, where D{} and D2{} are respectively the
*              first and second derivative operators. All of the values
*              can be scaled by a factor of scalFac, if desired, to help
*              prevent overflows with high degrees and orders. The
*              Helmholtz polynomials are used in Pine's method for
*              spherical harmonic evalulation.
*
*INPUTS: u The value at which the fully normalized Helmholtz polynomials
*          should be evaluated.
*   maxDeg The maximum degree and order of the output. This should be >=3.
* scalFactor A scale factor to help prevent overflow of the results. A
*          value of 10^(-280) (for example) might be useful at high
*          degrees.
*
*OUTPUTS: HBar An instance of the ClusterSet class such that
*              HBar(n+1,m+1)=scalFac*\bar{H}^m_n(u).
*      dHBardu An instance of the ClusterSet class such that
*              dHBardu(n+1,m+1)=scalFac*D{\bar{H}^m_n(u)}.
*   d2HBardu2  An instance of the ClusterSet class such that
*              dHBardu(n+1,m+1)=scalFac*D{\bar{H}^m_n(u)}.
*
*The fully normalized derived Legendre functions (Helmholtz polynomials)
*are described in [1] and the algorithm of that paper is implemented here
*with the addition of a scale factor, if desired.
*
*The Helmholtz polynomial of degree n and order m is defined to be
*H^m_n(u)=1/(n!2^n)D_{n+m}{(u^2-1)^n}
*where the notation D_y{x} represents the yth derivative of x. A fully
*normalized Helmholtz polynomial is defined as
*\bar{H}^m_n(u)=H^m_n(u)*sqrt(k*(2*n+1)*factorial(n-m)/factorial(n+m))
*where k=1 if m=0 and k=2 otherwise. This normalization is consistent with
*the normalization that the National Geospatial Intelligence Agency (NGA)
*uses with the coefficients for the EGM96 and EGM2008 gravitational models.
*
*This implementation tends to be approximately 600 times faster than the
*Matlab implementation.
*
*Note that this function requires a large amount of memory, since Matlab's
*mxSetProperty function copies the data rather than allowing pointers to
*be passed. Thus, more than double the normal amount of space is required
*for the return values.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[HBar,dHBardu,d2HBardu2]=normHelmholtz(u,M,scalFactor);
*
*REFERENCES:
*[1] E. Fantino and S. Casotto, "Methods of harmonic synthesis for global
*    geopotential models and their first-, second- and third-order
*    gradients," Journal of Geodesy, vol. 83, no. 7, pp. 595-619, Jul.
*    2009.
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
    double u, scalFactor;
    ClusterSetCPP<double> HBar;
    ClusterSetCPP<double> dHBardu;//The first derivatives
    ClusterSetCPP<double> d2HBardu2;//The second derivatives
    size_t M, numH, i;
    mxArray *CSRetVal;
    mxArray *clusterElsMATLAB,*clusterSizesMATLAB, *offsetArrayMATLAB;
    
    if(nrhs!=3){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    u=getDoubleFromMatlab(prhs[0]);
    M=getSizeTFromMatlab(prhs[1]);
    scalFactor=getDoubleFromMatlab(prhs[2]);
    
    if(M<3) {
       mexErrMsgTxt("The maximum order should be at least 3.");
       return; 
    }
    
    numH=(M+1)*(M+2)/2;
    
    //Allocate space for the results.
    clusterElsMATLAB=mxCreateDoubleMatrix(numH,1,mxREAL);
    clusterSizesMATLAB=allocUnsignedSizeMatInMatlab(M+1,1);
    offsetArrayMATLAB=allocUnsignedSizeMatInMatlab(M+1,1);
    
    HBar.numClust=M+1;
    HBar.totalNumEl=numH;
    HBar.clusterEls=reinterpret_cast<double*>(mxGetData(clusterElsMATLAB));
    HBar.offsetArray=reinterpret_cast<size_t*>(mxGetData(offsetArrayMATLAB));
    HBar.clusterSizes=reinterpret_cast<size_t*>(mxGetData(clusterSizesMATLAB));
    
    //Initialize the offset array and cluster sizes.
    HBar.offsetArray[0]=0;
    HBar.clusterSizes[0]=1;
    for(i=1;i<=M;i++){
        HBar.clusterSizes[i]=i+1;
        HBar.offsetArray[i]=HBar.offsetArray[i-1]+HBar.clusterSizes[i-1];
    }
    
    normHelmHoltzCPP(HBar,u,scalFactor);
    
    //Set the first return value
    mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
    mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
    mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
    mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);
    
    plhs[0]=CSRetVal;
    
    if(nlhs>1) {//Compute the first derivatives, if they are desired.
        mxArray *clusterEls1stDerivMATLAB=mxCreateDoubleMatrix(numH,1,mxREAL);
        
        dHBardu.numClust=M+1;
        dHBardu.totalNumEl=numH;
        dHBardu.clusterEls=reinterpret_cast<double*>(mxGetData(clusterEls1stDerivMATLAB));
        dHBardu.offsetArray=reinterpret_cast<size_t*>(mxGetData(offsetArrayMATLAB));
        dHBardu.clusterSizes=reinterpret_cast<size_t*>(mxGetData(clusterSizesMATLAB));
        
        normHelmHoltzDerivCPP(dHBardu,HBar);
        //Set the second return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls1stDerivMATLAB);
        mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
        mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);

        plhs[1]=CSRetVal;
        mxDestroyArray(clusterEls1stDerivMATLAB);
    }
    
    if(nlhs>2) {//Compute the second derivatives if they are desired.
        mxArray *clusterEls2ndDerivMATLAB=mxCreateDoubleMatrix(numH,1,mxREAL);
        
        d2HBardu2.numClust=M+1;
        d2HBardu2.totalNumEl=numH;
        d2HBardu2.clusterEls=reinterpret_cast<double*>(mxGetData(clusterEls2ndDerivMATLAB));
        d2HBardu2.offsetArray=reinterpret_cast<size_t*>(mxGetData(offsetArrayMATLAB));
        d2HBardu2.clusterSizes=reinterpret_cast<size_t*>(mxGetData(clusterSizesMATLAB));
        
        normHelmHoltzDeriv2CPP(d2HBardu2,HBar);
        
        //Set the third return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls2ndDerivMATLAB);
        mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
        mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);

        plhs[2]=CSRetVal;
        mxDestroyArray(clusterEls2ndDerivMATLAB);
    }
    
    //Free the buffers. The mxSetProperty command copied the data.
    mxDestroyArray(clusterElsMATLAB);
    mxDestroyArray(clusterSizesMATLAB);
    mxDestroyArray(offsetArrayMATLAB);
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
