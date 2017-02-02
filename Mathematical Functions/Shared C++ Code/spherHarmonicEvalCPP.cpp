/*SPHERHARMONICEVALCPP A C++ implementation of a function to determine a
 *                     potential and/ or the gradient of the potential from
 *                     spherical harmonic coefficients.
 *
 *This function is a C++ implementation of the main routine of the function
 *spherHarmonicEval in Matlab. If the correct files are compiled, the
 *Matlab function spherHarmonicEval will call this function through the C++
 *mex function spherHarmonicEvalCPPInt in place of running the
 *significantly slower Matlab routine. The function spherHarmonicEval
 *can be consulted for more information regarding the implementation and
 *the meaning of the results. 
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"
#include "CoordFuncs.hpp"

//For the sin and cos.
#include <cmath>
#include <limits>
//for memset
#include <cstring>

void spherHarmonicEvalCPP(double *V, double *gradV, const ClusterSetCPP<double> &C,const ClusterSetCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const double scalFactor) {
    //If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. Otherwise, a pointer to a buffer for 3
    //doubles should be passed.
    double temp, r, lambda, *nCoeff;
    const size_t M=C.numClust-1;
    const double pi = 2*acos(0.0);
    size_t n,m,curPoint;
    double nf,mf;
    double rPrev,thetaPrev;
    ClusterSetCPP<double> FuncVals;
    ClusterSetCPP<double> FuncDerivs;
    //These are to store sin(m*lambda) and cos(m*lambda) for m=0->M.
    double *SinVec,*CosVec;//This are each length C.numClust.
    //These are never used at the same time as SinVec andCosVec and are the
    //same size, so they will point to the same memory.
    double *rm,*im;
    //These hold values for the modified forward row algorithm that stay
    //constant for points with the same range and latitude but different
    //longitudes.
    double *XC,*XS;
    //These values are only used if gradV!=NULL. There are initialized to
    //NULL here to suppress a warning if compiled using
    //-Wconditional-uninitialized
    double *XCdr=NULL;
    double *XSdr=NULL;
    double *XCdTheta=NULL;
    double *XSdTheta=NULL;
    //A big chunk of memory will be allocated into a single buffer and
    //split between the variables that need it. That is faster than
    //allocating a bunch of small buffers, and all of the variables are of
    //the same type.
    double *buffer;
        
    //Initialize the ClusterSet classes for the coefficients. The space
    //for the elements will be allocated shortly.
    FuncVals.numClust=C.numClust;
    FuncVals.totalNumEl=C.totalNumEl;
    FuncVals.offsetArray=C.offsetArray;
    FuncVals.clusterSizes=C.clusterSizes;
    FuncDerivs.numClust=C.numClust;
    FuncDerivs.totalNumEl=C.totalNumEl;
    FuncDerivs.offsetArray=C.offsetArray;
    FuncDerivs.clusterSizes=C.clusterSizes;

    //Allocate the buffer and partition it between variables.
    if(gradV==NULL){
        buffer = new double[C.totalNumEl+5*C.numClust];
    }else{
        buffer = new double[2*C.totalNumEl+9*C.numClust];
    }
    {
        double *tempPtr=buffer;
        //This stores all of the powers of a/r needed for the sum, regardless
        //of which algorithm is used.
        nCoeff=tempPtr;
        tempPtr+=C.numClust;
        //The sin and cosine values use the same memory as the im and rm
        //values, because only one will be used depending on which
        //algorithm is executed.
        SinVec=tempPtr;
        rm=SinVec;
        tempPtr+=C.numClust;
        CosVec=tempPtr;
        im=CosVec;
        tempPtr+=C.numClust;
        FuncVals.clusterEls=tempPtr;
        tempPtr+=C.totalNumEl;
        XC=tempPtr;
        tempPtr+=C.numClust;
        XS=tempPtr;
                
        if(gradV!=NULL) {
            tempPtr+=C.numClust;
            FuncDerivs.clusterEls=tempPtr;
            tempPtr+=C.totalNumEl;
            XCdr=tempPtr;
            tempPtr+=C.numClust;
            XSdr=tempPtr;
            tempPtr+=C.numClust;
            XCdTheta=tempPtr;
            tempPtr+=C.numClust;
            XSdTheta=tempPtr;
        }
    }
        
    nCoeff[0]=1;
    
    rPrev=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double thetaCur;
        bool rChanged;
        bool thetaChanged;
        
        r=point[0+3*curPoint];
        lambda=point[1+3*curPoint];
        thetaCur=point[2+3*curPoint];
        
        rChanged=rPrev!=r;
        thetaChanged=thetaCur!=thetaPrev;
        rPrev=r;
        thetaPrev=thetaCur;
        
        if(rChanged) {
            temp=a/r;
            for(n=1;n<=M;n++) {
                nCoeff[n]=nCoeff[n-1]*temp;
            }
        }

        if(fabs(thetaCur)<88*pi/180||gradV==NULL) {
        //At latitudes that are not near the poles, the algorithm of Holmes and
        //Featherstone is used. It can not be used for the gradient near the
        //poles, because of the singularity of the spherical coordinate system.
            double u, theta;
            
            //Compute the sine and cosine terms 
            //Explicitely set the first two terms.
            SinVec[0]=0;
            CosVec[0]=1;
            SinVec[1]=sin(lambda);
            CosVec[1]=cos(lambda);
            //Use a double angle identity to get the second order term.
            SinVec[2]=2*SinVec[1]*CosVec[1];
            CosVec[2]=1-2*SinVec[1]*SinVec[1];
            //Use a two-part recursion for the rest of the terms.
            for(m=3;m<=M;m++){
                SinVec[m]=2*CosVec[1]*SinVec[m-1]-SinVec[m-2];
                CosVec[m]=2*CosVec[1]*CosVec[m-1]-CosVec[m-2];
            }
                        
        //The spherical coordinate system used in ellips2Sphere uses azimuth
        //and elevation (latitude). However the formulae for spherical harmonic
        //synthesis in the Holmes and Featherstone paper use, pi/2-elevation
        //(colatitude). Thus, the point must be transformed.
            theta=pi/2-thetaCur;
            u=sin(theta);
            if(thetaChanged) {
                //Get the associated Legendre function ratios.
                NALegendreCosRatCPP(FuncVals,  theta, scalFactor);

                //Get the derivatives of the ratios if the gradient is desired.
                if(gradV!=NULL) {
                    NALegendreCosRatDerivCPP(FuncDerivs,FuncVals,theta);
                    
                }
            }
            
            //Evaluate Equation 7 from the Holmes and Featherstone paper.
            if(rChanged||thetaChanged) {
                //Zero the arrays.
                memset(XC,0,sizeof(double)*C.numClust);
                memset(XS,0,sizeof(double)*C.numClust);

                //Compute the X coefficients for the sum
                for(m=0;m<=M;m++) {
                    for(n=m;n<=M;n++) {
                        XC[m]+=nCoeff[n]*C[n][m]*FuncVals[n][m];
                        XS[m]+=nCoeff[n]*S[n][m]*FuncVals[n][m];
                    }
                }
            }
            
            //Use Horner's method to compute V.
            V[curPoint]=0;
            m=M+1;
            do {
                m--;
                
                V[curPoint]=V[curPoint]*u+XC[m]*CosVec[m]+XS[m]*SinVec[m];
            } while(m>0);

            //Multiple by the constant in front of the sum and get rid of the
            //scale factor.
            V[curPoint]=(c/r)*V[curPoint]/scalFactor;

            //Compute the gradient, if it is desired.
            if(gradV!=NULL) {
                double J[9];
                double dVdr=0;
                double dVdLambda=0;
                double dVdTheta=0;
                
                if(rChanged||thetaChanged) {
                    memset(XCdr,0,sizeof(double)*C.numClust);
                    memset(XSdr,0,sizeof(double)*C.numClust);
                    memset(XCdTheta,0,sizeof(double)*C.numClust);
                    memset(XSdTheta,0,sizeof(double)*C.numClust);

                    //Evaluate Equation 7 from the Holmes and Featherstone paper.
                    mf=0;
                    for(m=0;m<=M;m++) {
                        nf=mf;
                        for(n=m;n<=M;n++) {
                            double CScal=nCoeff[n]*C[n][m];
                            double SScal=nCoeff[n]*S[n][m];

                            XCdr[m]+=(nf+1)*CScal*FuncVals[n][m];
                            XSdr[m]+=(nf+1)*SScal*FuncVals[n][m];

                            XCdTheta[m]+=CScal*FuncDerivs[n][m];
                            XSdTheta[m]+=SScal*FuncDerivs[n][m];
                                                        
                            nf++;
                        }
                                    
                        mf++;
                    }
                 }
                
                //Use Horner's method to compute the partials.
                m=M+1;
                mf=static_cast<double>(m);
                do {
                    m--;
                    mf--;
                    
                    dVdr=dVdr*u+XCdr[m]*CosVec[m]+XSdr[m]*SinVec[m];
                    dVdLambda=dVdLambda*u+mf*(-XC[m]*SinVec[m]+XS[m]*CosVec[m]);
                    dVdTheta=dVdTheta*u+XCdTheta[m]*CosVec[m]+XSdTheta[m]*SinVec[m];
                } while(m>0);            

                dVdr=-(c/(r*r))*dVdr/scalFactor;
                dVdLambda=(c/r)*dVdLambda/scalFactor;
            //The minus sign is because the input coordinate was with respect
            //to latitude, not the co-latitude that the NALegendreCosRat
            //function uses.
                dVdTheta=-(c/r)*dVdTheta/scalFactor;

                calcSpherJacobCPP(J, point+3*curPoint,0);

                //Now, multiply the transpose of the Jacobian Matrix by the
                //vector of [dVdr;dVdLambda;dVdTheta]
                gradV[0+3*curPoint]=dVdr*J[0]+dVdLambda*J[1]+dVdTheta*J[2];
                gradV[1+3*curPoint]=dVdr*J[3]+dVdLambda*J[4]+dVdTheta*J[5];
                gradV[2+3*curPoint]=dVdr*J[6]+dVdLambda*J[7]+dVdTheta*J[8];
            }
        } else {  
        //At latitudes that are near the poles, the non-singular algorithm of
        //Pines using the fully normalized Helmholtz equations from Fantino and
        //Casotto is used. The algorithm has been slightly modified so that the
        //c/r term is out front and the fully normalized Helmholtz polynomials
        //can be scaled. Also, lumped coefficients are not used. The Pines
        //algorithm is generally slower than the algorithm of Holmes and
        //Featherstone and it is suffers a loss of precision near the equator.
        //Thus, the Pines algorithm is only used near the poles where the other
        //algorithm has issues with a singularity.
            double s,t,u, CartPoint[3];

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
            V[curPoint]=0;
            for(n=0;n<=M;n++) {
                double innerTerm=0;
                for(m=0;m<=n;m++) {
                    innerTerm+=(C[n][m]*rm[m]+S[n][m]*im[m])*FuncVals[n][m];
                }
                V[curPoint]+=nCoeff[n]*innerTerm;
            }

            V[curPoint]=(c/r)*V[curPoint]/scalFactor;

            //Compute the gradient.
            if(gradV!=NULL) {
                double a1=0;
                double a2=0;
                double a3=0;
                double a4=0;

                normHelmHoltzDerivCPP(FuncDerivs,FuncVals);

                //The equations in these loops are from Table 10.
                nf=0.0;
                for(n=0;n<=M;n++) {
                    double a1Loop=0;
                    double a2Loop=0;
                    double a3Loop;
                    double a4Loop;
                    double CProdMN;
                    double HVal;
                    double dHVal;
                    double Lmn;

                    //The m=0 case only applies to a3 and a4.
                    m=0;
                    mf=0.0;
                    HVal=FuncVals[n][m];
                    dHVal=FuncDerivs[n][m];
                    CProdMN=C[n][m]*rm[m]+S[n][m]*im[m];

                    a3Loop=CProdMN*dHVal;
                    Lmn=(nf+mf+1)*HVal+u*dHVal;//Defined in Table 14.
                    a4Loop=-CProdMN*Lmn;

                    mf=1.0;
                    for(m=1;m<=n;m++) {
                        HVal=FuncVals[n][m];
                        dHVal=FuncDerivs[n][m];

                        a1Loop+=mf*(C[n][m]*rm[m-1]+S[n][m]*im[m-1])*HVal;
                        a2Loop+=mf*(S[n][m]*rm[m-1]-C[n][m]*im[m-1])*HVal;

                        CProdMN=C[n][m]*rm[m]+S[n][m]*im[m];

                        a3Loop+=CProdMN*dHVal;
                        Lmn=(nf+mf+1)*HVal+u*dHVal;
                        a4Loop-=CProdMN*Lmn;

                        mf++;
                    }

                    a1+=nCoeff[n]*a1Loop;
                    a2+=nCoeff[n]*a2Loop;
                    a3+=nCoeff[n]*a3Loop;
                    a4+=nCoeff[n]*a4Loop;

                    nf++;
                }

       //These are equation 70. However, an additional 1/r term has been added,
       //which the original paper omitted when going from Equation 68 to 70.
                {
                    temp=c/(r*r);
                    gradV[0+3*curPoint]=temp*(a1+s*a4)/scalFactor;
                    gradV[1+3*curPoint]=temp*(a2+t*a4)/scalFactor;
                    gradV[2+3*curPoint]=temp*(a3+u*a4)/scalFactor;
                }
            }
        }
    }
    
    delete[] buffer;
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
