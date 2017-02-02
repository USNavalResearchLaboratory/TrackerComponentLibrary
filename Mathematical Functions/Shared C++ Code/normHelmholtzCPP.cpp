/*NORMHELMHOLTZCPP C++ implementations of functions to computes fully
 *                 normalized derived Legendre functions (also known as
 *                 fully normalized Helmholtz  polynomials) and their first
 *                 and second derivatives.
 *
 *The three functions contained in this file correspond to the results of
 *the function normHelmholtz in Matlab. The function normHelmholtz
 *can be consulted for more information regarding the implementation and
 *the meaning of the results. The functions below place their results in
 *ClusterSet classes that must be preallocated and provied to the
 *functions and that hold the values indexed by degree and order. 
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"
//For sin and cos.
#include <math.h>

void normHelmHoltzCPP(ClusterSetCPP<double> &HBar,const double u, const double scalFactor) {
    const size_t M=HBar.numClust-1;
    double g,h;
    size_t n,m;
    double nf, mf;
    
    //Set the first few terms explicitely.
    HBar[0][0]=1.0*scalFactor;
    HBar[1][1]=sqrt(3.0)*scalFactor;
    HBar[1][0]=sqrt(3.0)*u*scalFactor;

    //Compute all terms of the form (n,n) and (n,n-1).
    nf=2.0;
    for(n=2;n<=M;n++) {
        //Get the (n,n) term using Equation 55.
        double f=sqrt((2*nf+1)/(2*nf));
        HBar[n][n]=f*HBar[n-1][n-1];
        //Get the (n,n-1) term using Equation 56.
        HBar[n][n-1]=u*sqrt(2*nf)*HBar[n][n];
        
        nf++;
    }
    
    //Now, compute all of the other terms using Equation 58.
    mf=0;
    for(m=0;m<=M;m++) {
        nf=mf+2;
        for(n=(m+2);n<=M;n++) {
            g=sqrt((2*nf+1)*(2*nf-1)/((nf+mf)*(nf-mf)));
            h=sqrt((2*nf+1)*(nf-mf-1)*(nf+mf-1)/((2*nf-3)*(nf+mf)*(nf-mf)));
            
            HBar[n][m]=u*g*HBar[n-1][m]-h*HBar[n-2][m];
            
            nf++;
        }
                    
        mf++;
    }
}

void normHelmHoltzDerivCPP(ClusterSetCPP<double> &dHBardu,const ClusterSetCPP<double> &HBar) {
    const size_t M=HBar.numClust-1;
    double k;
    size_t n, m;
    double nf, mf;
    //The first derivatives are easily computed using Equation 61.

    //For m=0
    dHBardu[0][0]=0;
    nf=1.0;
    for(n=1;n<=M;n++) {
        k=sqrt(nf*(nf+1)/2);
        dHBardu[n][0]=k*HBar[n][1];
        
        nf++;
    }
    
    //For m=1 and n=1;
    n=1;
    m=1;
    dHBardu[n][m]=0;
    
    //For other n and m pairs
    nf=2.0;
    for(n=2;n<=M;n++) {
        mf=1.0;
        for(m=1;m<=(n-1);m++) {
            k=sqrt((nf-mf)*(nf+mf+1));
            
            dHBardu[n][m]=k*HBar[n][m+1];
            
            mf++;
        }
        dHBardu[n][n]=0;
        
        nf++;
    }
}

void normHelmHoltzDeriv2CPP(ClusterSetCPP<double> &d2HBardu2,const ClusterSetCPP<double> &HBar) {
    const size_t M=HBar.numClust-1;
    double k, kp;
    size_t n, m;
    double nf, mf;
    //The first derivatives are easily computed using Equation 62.
    
    //For m=0 and m=1
    d2HBardu2[0][0]=0;
    d2HBardu2[1][0]=0;
    d2HBardu2[1][1]=0;

    //For n=2, m=0.
    m=0;
    n=2;
    mf=0.0;
    nf=2.0;
    k=sqrt((nf-mf)*(nf+mf+1));
    kp=sqrt((nf-(mf+1))*(nf+(mf+1)+1));
    d2HBardu2[n][m]=k*kp*HBar[n][m+2];

    //For n=2, m=1;
    m=1;
    mf=1.0;
    d2HBardu2[n][m]=0;

    //For m=0 and m=1 and all other n.
    nf=3.0;
    for(n=3;n<=M;n++) {          
        m=0;
        mf=0.0;
        k=sqrt(nf*(nf+1)/2);
        kp=sqrt((nf-(mf+1))*(nf+(mf+1)+1));
        d2HBardu2[n][m]=k*kp*HBar[n][m+2];

        m=1;
        mf=1.0;
        k=sqrt((nf-mf)*(nf+mf+1));
        kp=sqrt((nf-(mf+1))*(nf+(mf+1)+1));
        d2HBardu2[n][m]=k*kp*HBar[n][m+2];
        
        nf++;
    }

    //For other n and m pairs.
    nf=3.0;
    for(n=3;n<=M;n++) {
        mf=2.0;
        for(m=2;m<=n-2;m++) {
            k=sqrt((nf-mf)*(nf+mf+1));
            kp=sqrt((nf-(mf+1))*(nf+(mf+1)+1));

            d2HBardu2[n][m]=k*kp*HBar[n][m+2];
            
            mf++;
        }

        d2HBardu2[n-1][n-1]=0;
        d2HBardu2[n][n-1]=0;
        d2HBardu2[n][n]=0;
        
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
