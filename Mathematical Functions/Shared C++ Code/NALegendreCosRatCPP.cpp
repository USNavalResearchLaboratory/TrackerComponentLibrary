/*NALEGENDRECOSRATCPP C++ implementations of functions to evaluate the
 *                    ratio of fully normalized associated Legendre
 *                    functions of the cosine of a parameter to u^m,
 *                    where u is the sine of a parameter and m is the
 *                    order of the function.
 *
 *The three functions contained in this file correspond to the results of
 *the function NALegendreCosRat in Matlab. The function NALegendreCosRat
 *can be consulted for more information regarding the implementation and
 *the meaning of the results. The functions below place their results in
 *ClusterSet classes that must be preallocated and provied to the
 *functions and that hold the values indexed by degree and order. See the
 *Matlab implementation of the function NALegendreCosRat for more
 *information on the algorithm used.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"
//For sin and cos.
#include <math.h>

void NALegendreCosRatCPP(ClusterSetCPP<double> &PBarUVals, const double theta, const double scalFactor) {
//It is assumed that space for the results is already preallocated in 
//PBarUVals along with the proper offset array.
    const double u=sin(theta);
    const double t=cos(theta);
    const double jTerm=1/sqrt(2.0);
    const size_t M=PBarUVals.numClust-1;
    double g, h;
    size_t numPBarU;
    size_t n, m;
    //Floating point versions of the integer indices are needed for
    //computing coefficients.
    double mf, nf;

    numPBarU=(M+1)*(M+2)/2;
    
    //The value of PBar_{0,0}(cos(theta)) is independent of theta and is
    //one. 
    PBarUVals[0][0]=1.0*scalFactor;

    //Set the seed value for PBar_{1,1}(cos(theta))/u from which the other
    //values will be computed.
    PBarUVals[1][1]=sqrt(3.0)*scalFactor;

    //Compute the values along the main diagonal, where m=n starting from
    //m=n=2. This implements equation 28 in the first Holmes and
    //Featherstone paper for the normalized associated Legendre function 
    //ratio.
    mf=2.0;
    for(m=2;m<=M;m++) {
        PBarUVals[m][m]=sqrt((2*mf+1)/(2*mf))*PBarUVals[m-1][m-1];
        
        mf++;
    }

    //Recursively compute the values using Equation 27 from the first Holmes
    //and Featherstone paper, taking into account the fact that the first
    //element of the recursion only has one term.
    
    //First, deal with the case where n=1, m=0;
    n=1;
    nf=1.0;
    m=0;
    mf=0.0;
    //g is given in Equation 18 of the first Holmes and Featherstone paper.
    g=2*(mf+1)/sqrt((nf-mf)*(nf+mf+1));
    PBarUVals[n][m]=jTerm*g*t*PBarUVals[n][m+1];
    
    //Next, evaluate the values for all other valid n and m.
    nf=2.0;
    for(n=2;n<=M;n++) {
        //Deal with the first element of the recursion,which is  where m=n-1.
        m=n-1;
        mf=nf-1.0;
        
        g=2*(mf+1)/sqrt((nf-mf)*(nf+mf+1));
        PBarUVals[n][m]=g*t*PBarUVals[n][m+1];
        
        //Recursively compute the values of the rest of the m terms.
        mf=nf-2.0;
        for(m=n-2;m>0;m--) {
            g=2*(mf+1)/sqrt((nf-mf)*(nf+mf+1));
            //h is given in Equation 19 of the first Holmes and Featherstone
            //paper.
            h=sqrt((nf+mf+2)*(nf-mf-1)/((nf-mf)*(nf+mf+1)));
            PBarUVals[n][m]=g*t*PBarUVals[n][m+1]-h*u*u*PBarUVals[n][m+2];
            
            mf--;
        }

        //Deal with the special m=0 case.
        m=0;
        mf=0.0;
        g=2*(mf+1)/sqrt((nf-mf)*(nf+mf+1));
        h=sqrt((nf+mf+2)*(nf-mf-1)/((nf-mf)*(nf+mf+1)));
        PBarUVals[n][m]=jTerm*(g*t*PBarUVals[n][m+1]-h*u*u*PBarUVals[n][m+2]);
        
        nf++;
    }
}


void NALegendreCosRatDerivCPP(ClusterSetCPP<double> &dPBarUValsdTheta, const ClusterSetCPP<double> &PBarUVals, const double theta) {
    const double u=sin(theta);
    const double t=cos(theta);
    double e;
    const size_t M=PBarUVals.numClust-1;
    size_t n, m;
    double nf, mf;
    
    //The first derivative of PBar_{0,0}(cos(theta)) is just zero.
    dPBarUValsdTheta[0][0]=0;

    n=1;
    m=1;
    nf=1.0;
    mf=1.0;
    
    //From Equation 30 in the first Holmes and Featherstone paper. This
    //is the seed value from which other values will be computed.
    dPBarUValsdTheta[1][1]=mf*(t/u)*PBarUVals[1][1];

    //Compute the values along the main diagonal, where m=n starting 
    //from m=n=2. This implements Equation 30 in the first Holmes and
    //Featherstone paper for the ratio of the first derivative.
    mf=2.0;
    for(m=2;m<=M;m++) {
        dPBarUValsdTheta[m][m]=mf*(t/u)*PBarUVals[m][m];
        
        mf++;
    }

    n=1;
    m=0;
    nf=1.0;
    mf=0.0;
    
    //e is given in Equation 22 of the first Holmes and Featherstone
    //paper.
    e=sqrt((nf+mf+1)*(nf-mf)/2);
    //This is Equation 30 of the first Holmes and Featherstone paper for
    //m=0.
    dPBarUValsdTheta[n][m]=-e*u*PBarUVals[n][m+1];

    //Next, evaluate the values for all other valid n and m.
    nf=2.0;
    for(n=2;n<=M;n++) {
        //Recursively compute the values of the m terms for m>0.
        mf=nf-1;
        for(m=(n-1);m>0;m--){
            e=sqrt((nf+mf+1)*(nf-mf));
            dPBarUValsdTheta[n][m]=mf*(t/u)*PBarUVals[n][m]-e*u*PBarUVals[n][m+1];
            
            mf--;
        }
        //Deal with the special m=0 case.
        m=0;
        mf=0.0;
        e=sqrt((nf+mf+1)*(nf-mf)/2);
        dPBarUValsdTheta[n][m]=-e*u*PBarUVals[n][m+1];
        
        nf++;
    }
}

void NALegendreCosRatDeriv2CPP(ClusterSetCPP<double> &d2PBarUValsdTheta2, const ClusterSetCPP<double> &dPBarUValsdTheta, const ClusterSetCPP<double> &PBarUVals, const double theta) {    
    const double u=sin(theta);
    const double t=cos(theta);
    size_t n, m;
    double nf, mf;
    const size_t M=PBarUVals.numClust-1;

    nf=0.0;
    for(n=0;n<=M;n++) {
        mf=0.0;
        for(m=0;m<=n;m++) {
            //From the first (un-numbered) equation in the second Holmes
            //and Featherstone paper AFTER correction.
            d2PBarUValsdTheta2[n][m]=(mf*mf/(u*u)-nf*(nf+1))*PBarUVals[n][m]-(t/u)*dPBarUValsdTheta[n][m];
            
            mf++;
        }
        
        nf++;
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
