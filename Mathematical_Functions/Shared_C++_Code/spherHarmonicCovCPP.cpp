/*SPHERHARMONICCOVCPP A C++ implementation of a function to find the
 *                    standard deviation of a potential and the covariance
 *                    matrix of the gradient of the potential when then
 *                    potential is given by spherical harmonic
 *                    coefficients.
 *
 *This function is a C++ implementation of the main routine of the function
 *spherHarmonicEval in Matlab. If the correct files are compiled, the
 *Matlab function spherHarmonicCov will call this function through the C++
 *mex function spherHarmonicCovCPPInt in place of running the
 *significantly slower Matlab routine. The function spherHarmonicCov
 *can be consulted for more information regarding the implementation and
 *the meaning of the results. 
 *
 *April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"
#include "CoordFuncs.hpp"

//For the sin, cos and infinity.
#include <cmath>
#include <limits>
//for memset
#include <string.h>

//Many versions of Windows and some other operating systems do not support
//the isfinite function in C++, so we just define our own isFinite function
//to take its place in case it is not defined in cmath.
template<typename T> bool isFinite(T arg)
{
    return arg == arg && 
           arg != std::numeric_limits<T>::infinity() &&
           arg != -std::numeric_limits<T>::infinity();
}

bool spherHarmonicCovCPP(double *sigma2, double *Sigma, const CountingClusterSetCPP<double> &CStdDev,const CountingClusterSetCPP<double> &SStdDev, const double *point, const size_t numPoints, const double a, const double c, const double scalFactor) {
    //If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. Otherwise, a pointer to a buffer for 3
    //doubles should be passed. The return value indicates whether any
    //terms were discarded due to overflow errors.
    double r, *nCoeff;
    const size_t M=CStdDev.numClust-1;
    size_t n,m,curPoint;
    double nf,mf;
    double rPrevVal,thetaPrev;
    CountingClusterSetCPP<double> FuncVals;
    CountingClusterSetCPP<double> FuncDerivs;
    double *rm,*im;
    //A big chunk of memory will be allocated into a single buffer and
    //split between the variables that need it. That is faster than
    //allocating a bunch of small buffers, and all of the variables are of
    //the same type.
    double *buffer;
    bool didOverflow=false;//The return value
        
    //Initialize the ClusterSet classes for the coefficients. The space
    //for the elements will be allocated shortly.
    FuncVals.numClust=CStdDev.numClust;
    FuncVals.totalNumEl=CStdDev.totalNumEl;
    FuncDerivs.numClust=CStdDev.numClust;
    FuncDerivs.totalNumEl=CStdDev.totalNumEl;

    //Allocate the buffer and partition it between variables.
    if(Sigma==NULL){
        buffer = new double[CStdDev.totalNumEl+3*CStdDev.numClust];
    }else{
        buffer = new double[2*CStdDev.totalNumEl+3*CStdDev.numClust];
    }
    {
        double *tempPtr=buffer;
        //This stores all of the powers of a/r needed for the sum, regardless
        //of which algorithm is used.
        nCoeff=tempPtr;
        tempPtr+=CStdDev.numClust;
        rm=tempPtr;
        tempPtr+=CStdDev.numClust;
        im=tempPtr;
        tempPtr+=CStdDev.numClust;
        FuncVals.clusterEls=tempPtr;
             
        if(Sigma!=NULL) {
            tempPtr+=CStdDev.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;
        }
    }
        
    nCoeff[0]=1;
    
    rPrevVal=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double thetaCur;
        bool rChanged;
        bool thetaChanged;
        double s,t,u, CartPoint[3];
        
        r=point[0+3*curPoint];
        thetaCur=point[2+3*curPoint];
        
        rChanged=rPrevVal!=r;
        thetaChanged=thetaCur!=thetaPrev;
        rPrevVal=r;
        thetaPrev=thetaCur;
        
        if(rChanged) {
            double temp=a/r;
            for(n=1;n<=M;n++) {
                nCoeff[n]=nCoeff[n-1]*temp;
            }
        }

    //The non-singular algorithm of Pines using the fully normalized
    //Helmholtz equations from Fantino and Casotto is used. The algorithm
    //has been slightly modified so that the c/r term is out front and the
    //fully normalized Helmholtz polynomials can be scaled. Also, lumped
    //coefficients are not used. The Pines algorithm can suffer a loss of
    //precision near the equator. However, it is simple to just square the
    //terms in the sum.
        spher2CartCPP(CartPoint,point+3*curPoint,0);

        //Get the direction cosines used by Pines' algorithm.
        s=CartPoint[0]/r;
        t=CartPoint[1]/r;
        u=CartPoint[2]/r;

        //Compute the fully normalized Helmholtz polynomials.
        if(thetaChanged) {
            normHelmHoltzCPP(FuncVals,u, scalFactor);
            thetaPrev=thetaCur;
        }

        //Recursively compute the rm and im terms for the sums.
        rm[0]=1;
        im[0]=0;
        for(m=1;m<=M;m++) {
            //These are equation 49 in the Fantino and Casotto paper.
            rm[m]=s*rm[m-1]-t*im[m-1];
            im[m]=s*im[m-1]+t*rm[m-1];
        }

        //Perform the sum for the potential from Equation 44 in the
        //Fantino and Casotto paper.
        sigma2[curPoint]=0;
        for(n=0;n<=M;n++) {
            double innerTerm=0;
            for(m=0;m<=n;m++) {
                double CVal, SVal, rVal, iVal, FVal;
                double temp1, temp2;
                CVal=CStdDev[n][m];
                SVal=SStdDev[n][m];
                rVal=rm[m];
                iVal=im[m];
                FVal=FuncVals[n][m];
                
                temp1=CVal*rVal*FVal;
                temp2=SVal*iVal*FVal;
                
                innerTerm+=temp1*temp1+temp2*temp2;
            }
            if(isFinite(innerTerm)){
                sigma2[curPoint]+=nCoeff[n]*nCoeff[n]*innerTerm;
            } else {
                didOverflow=true;
            }
        }
        
        //The variance of the potential.
        sigma2[curPoint]=(c/r)*(c/r)*sigma2[curPoint]/(scalFactor*scalFactor);

        //Compute the covariance matrix of the gradient, if needed.
        if(Sigma!=NULL) {
            double a11=0;
            double a22=0;
            double a33=0;
            double a44=0;
            double a12=0;
            double a13=0;
            double a14=0;
            double a23=0;
            double a24=0;
            double a34=0;

            normHelmHoltzDerivCPP(FuncDerivs,FuncVals);

            //The equations in these loops are from Table 10.
            nf=0.0;
            for(n=0;n<=M;n++) {
                double CProdMN;
                double HVal;
                double dHVal;
                double Lmn;
                double CCur2,SCur2,rCur,iCur;
                
                double a11Loop=0;
                double a22Loop=0;
                double a33Loop;
                double a12Loop=0;
                double a13Loop=0;
                double a14Loop=0;
                double a23Loop=0;
                double a24Loop=0;
                double a34Loop;
                double a44Loop;

                //The m=0 case only applies to a3 and a4, so that means only to
                //a33, a34, and a44.
                m=0;
                mf=0.0;
                HVal=FuncVals[n][m];
                dHVal=FuncDerivs[n][m];
                CCur2=CStdDev[n][m]*CStdDev[n][m];
                SCur2=SStdDev[n][m]*SStdDev[n][m];
                rCur=rm[m];
                iCur=im[m];
                
                CProdMN=CCur2*rCur*rCur+SCur2*iCur*iCur;
                Lmn=(nf+mf+1)*HVal+u*dHVal;//Defined in Table 14.
                
                a33Loop=CProdMN*dHVal*dHVal;
                a44Loop=CProdMN*Lmn*Lmn;
                a34Loop=-CProdMN*Lmn*dHVal;
                
                mf=1.0;
                for(m=1;m<=n;m++) {
                    const double rPrev=rm[m-1];
                    const double iPrev=im[m-1];
                    
                    HVal=FuncVals[n][m];
                    dHVal=FuncDerivs[n][m];
                    CCur2=CStdDev[n][m]*CStdDev[n][m];
                    SCur2=SStdDev[n][m]*SStdDev[n][m];
                    rCur=rm[m];
                    iCur=im[m];
                    
                    CProdMN=CCur2*rCur*rCur+SCur2*iCur*iCur;
                    Lmn=(nf+mf+1)*HVal+u*dHVal;
                    
                    //These if-statements are to deal with numerical
                    //precision problems near the poles. We want to avoid
                    //0*Inf terms due to limitations in the valid range of
                    //double precision numbers. Of course, the loss of the
                    //terms where overflow occurs means that the covariance
                    //matrix will be underestimated.
                    
                    if(isFinite(HVal)) {
                        a11Loop+=mf*mf*(CCur2*rPrev*rPrev+SCur2*iPrev*iPrev)*HVal*HVal;
                        a12Loop+=mf*mf*rPrev*iPrev*(SCur2-CCur2)*HVal*HVal;
                        a22Loop+=mf*mf*(SCur2*rPrev*rPrev+CCur2*iPrev*iPrev)*HVal*HVal;
                    } else {
                        didOverflow=true;
                    }

                    if(isFinite(Lmn)) {
                        a44Loop+=CProdMN*Lmn*Lmn;
                        if(isFinite(HVal)) {
                            a14Loop-=mf*(CCur2*rPrev*rCur+SCur2*iPrev*iCur)*HVal*Lmn;
                            a24Loop-=mf*(-CCur2*iPrev*rCur+SCur2*rPrev*iCur)*HVal*Lmn;
                        }
                        if(isFinite(dHVal)){
                            a34Loop-=CProdMN*Lmn*dHVal;
                        }
                    } else {
                        didOverflow=true;
                    }
                    
                    if(isFinite(dHVal)) {
                        a33Loop+=CProdMN*dHVal*dHVal;
                        if(isFinite(HVal)) {
                            a13Loop+=mf*(CCur2*rPrev*rCur+SCur2*iPrev*iCur)*HVal*dHVal;
                            a23Loop+=mf*(-CCur2*iPrev*rCur+SCur2*rPrev*iCur)*HVal*dHVal;
                        }
                    } else {
                        didOverflow=true;
                    }
                    
                    mf++;
                }

                {
                    const double nCoeff2=nCoeff[n]*nCoeff[n];
                    
                    a11+=nCoeff2*a11Loop;
                    a22+=nCoeff2*a22Loop;
                    a33+=nCoeff2*a33Loop;
                    a44+=nCoeff2*a44Loop;
                    a12+=nCoeff2*a12Loop;
                    a13+=nCoeff2*a13Loop;
                    a14+=nCoeff2*a14Loop;
                    a23+=nCoeff2*a23Loop;
                    a24+=nCoeff2*a24Loop;
                    a34+=nCoeff2*a34Loop;
                }

                nf++;
            }

//These are based on squaring the terms in equation 70, removing cross
//terms. However, an additional 1/r (squared) term has been added,
//which the original paper omitted when going from Equation 68 to 70.
            {
                double temp=c/(r*r*scalFactor);
                double s11=a11+2*s*a14+s*s*a44;
                double s12=a12+s*a24+t*a14+s*t*a44;
                double s13=a13+s*a34+u*a14+s*u*a44;
                double s22=a22+2*t*a24+t*t*a44;
                double s23=a23+t*a34+u*a24+t*u*a44;
                double s33=a33+2*u*a34+u*u*a44;
                
                *Sigma=temp*(temp*s11);
                *(Sigma+1)=temp*(temp*s12);
                *(Sigma+2)=temp*(temp*s13);
                *(Sigma+3)=temp*(temp*s12);
                *(Sigma+4)=temp*(temp*s22);
                *(Sigma+5)=temp*(temp*s23);
                *(Sigma+6)=temp*(temp*s13);
                *(Sigma+7)=temp*(temp*s23);
                *(Sigma+8)=temp*(temp*s33);
                
                //Go to the next point.
                Sigma+=9;
            }
        }
    }
    
    delete[] buffer;
    
    return didOverflow;
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
