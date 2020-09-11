/*SPHERHARMONICEVAL A C++ implementation of functions to determine a
 *                  potential and/ or the gradient of the potential from
 *                  ral or complex spherical harmonic coefficients.
 *
 *This function is a C++ implementation of the main routine of the function
 *spherHarmonicEval in Matlab. If the correct files are compiled, the
 *Matlab function spherHarmonicEval will call this function through the C++
 *mex function spherHarmonicEvalCPP in place of running the
 *significantly slower Matlab routine. The function spherHarmonicEval
 *can be consulted for more information regarding the implementation and
 *the meaning of the results. 
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathFuncs.hpp"
#include "CoordFuncs.hpp"

//For the sin and cos.
#include <cmath>
#include <limits>
//for memset
#include <cstring>
//For the implementation in the complex domain.
#include <complex>

using namespace std;

void spherHarmonicEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetCPP<double> &C,const CountingClusterSetCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm) {
    //If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. If a NULL pointer is passed for HessianV,
    //then it is assumed that the Hessian is not desired.
    const size_t M=C.numClust-1;
    const size_t M1=C.numClust;
    const double pi=2.0*acos(0.0);
    size_t curPoint;
    double rPrev, thetaPrev, crScal;
    CountingClusterSetCPP<double> FuncVals;
    CountingClusterSetCPP<double> FuncDerivs;
    CountingClusterSetCPP<double> FuncDerivs2;
    //A big chunk of memory will be allocated into a single buffer and
    //split between the variables that need it. That is faster than
    //allocating a bunch of small buffers, and all of the variables are of
    //the same type.
    double *buffer;
    //This stores all of the powers of a/r needed for the sum, regardless
    //of which algorithm is used.
    double *nCoeff;//Length M+1
    //These are to store sin(m*lambda) and cos(m*lambda) for m=0->M for the
    //Legendre algorithm.
    double *SinVec,*CosVec;//These are each length M+1.
    //rm and im are never used at the same time as SinVec andCosVec,
    //because rm and im are used by Pines' algorithm. They will be the same
    //size as SinVec and CosVec, and thus will point to the same memory.
    double *rm,*im;
    //The following are needed if the Legendre algorithm is used.
    double *A=NULL;//Length M1
    double *B=NULL;//Length M1
    //The following are only used if the gradient or Hessian is desired and
    //the Legendre algorithm is used.
    double *Ar=NULL;//Length M1
    double *Br=NULL;//Length M1
    double *ATheta=NULL;//Length M1
    double *BTheta=NULL;//Length M1
    //The following are only used if the Hessian is desired and the
    //Legendre algorithm is used.
    double *AThetaTheta=NULL;//Length M1
    double *BThetaTheta=NULL;//Length M1
    double *Arr=NULL;//Length M1
    double *Brr=NULL;//Length M1
    double *AThetar=NULL;//Length M1
    double *BThetar=NULL;//Length M1

    //Initialize the ClusterSet classes for the coefficients. The space
    //for the elements will be allocated shortly.
    //FuncVals will be used for HBar and PBarUVals.
    //FuncDerivs will be used for dPBarUValsdTheta and dHBardu.
    //FuncDerivs2 will be used for d2PBarUValsdTheta2 and d2HBardu2
    FuncVals.numClust=C.numClust;
    FuncVals.totalNumEl=C.totalNumEl;
    FuncDerivs.numClust=C.numClust;
    FuncDerivs.totalNumEl=C.totalNumEl;
    FuncDerivs2.numClust=C.numClust;
    FuncDerivs2.totalNumEl=C.totalNumEl;
    
    if(algorithm==2) {//If only Pines' algorithm is used.
        //This pointer is used to help parition the memory in buffer.
        double *tempPtr;

        if(gradV==NULL&&HessianV==NULL) {
            buffer=new double[3*M1+C.totalNumEl];
        }else{
            buffer=new double[3*M1+3*C.totalNumEl];
        }
        
        tempPtr=buffer;
        nCoeff=tempPtr;//Length M1

        tempPtr+=M1;
        rm=tempPtr;//Length M1

        tempPtr+=M1;
        im=tempPtr;//Length M1

        tempPtr+=M1;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        if(!(gradV==NULL&&HessianV==NULL)) {
            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl
            
            tempPtr+=C.totalNumEl;
            FuncDerivs2.clusterEls=tempPtr;// Length C.totalNumEl
        }
    } else {
        double *tempPtr;

        //If the Legendre algorithm might be used, then the amount of
        //memory needed varies depending on the highest order derivative
        //that is needed.
        if(HessianV!=NULL) {
            buffer=new double[15*M1+3*C.totalNumEl];
        } else if(gradV!=NULL) {
            if(algorithm==0) {//If Pines' algorithm might be used
                buffer=new double[9*M1+3*C.totalNumEl];
            } else {
                buffer=new double[9*M1+2*C.totalNumEl];
            }
        } else {
            buffer=new double[5*M1+C.totalNumEl];
        }
        
        tempPtr=buffer;
        nCoeff=tempPtr;//Length M1

        tempPtr+=M1;
        rm=tempPtr;
        SinVec=tempPtr;//Length M1

        tempPtr+=M1;
        im=tempPtr;
        CosVec=tempPtr;//Length M1

        tempPtr+=M1;
        A=tempPtr;//Length M1

        tempPtr+=M1;
        B=tempPtr;//Length M1

        tempPtr+=M1;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
            
        if(gradV!=NULL||HessianV!=NULL) {
            tempPtr+=C.totalNumEl;
            Ar=tempPtr;//Length M1

            tempPtr+=M1;
            Br=tempPtr;//Length M1

            tempPtr+=M1;
            ATheta=tempPtr;//Length M1

            tempPtr+=M1;
            BTheta=tempPtr;//Length M1

            tempPtr+=M1;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl

            if(HessianV!=NULL||algorithm!=0) {
                tempPtr+=C.totalNumEl;
                FuncDerivs2.clusterEls=tempPtr;//Length C.totalNumEl
            }

            if(HessianV!=NULL) {
                tempPtr+=C.totalNumEl;
                AThetaTheta=tempPtr;//Length M1

                tempPtr+=M1;
                BThetaTheta=tempPtr;//Length M1

                tempPtr+=M1;
                Arr=tempPtr;//Length M1

                tempPtr+=M1;
                Brr=tempPtr;//Length M1

                tempPtr+=M1;
                AThetar=tempPtr;//Length M1

                tempPtr+=M1;
                BThetar=tempPtr;//Length M1
            }
        }
    }

    nCoeff[0]=1;
    
    rPrev=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double pointCur[3];
        double r, lambda, thetaCur;
        size_t n,m;
        double nf,mf;
        bool useLegendre;
        bool rChanged;
        bool thetaChanged;
        
        pointCur[0]=point[0+3*curPoint];
        pointCur[1]=point[1+3*curPoint];
        if(systemType==2) {
            pointCur[2]=pi/2.0-point[2+3*curPoint];
        } else {
            pointCur[2]=point[2+3*curPoint];
        }
        
        r=pointCur[0];
        lambda=pointCur[1];
        thetaCur=pointCur[2];
        
        rChanged=(rPrev!=r);
        thetaChanged=(thetaCur!=thetaPrev);
        rPrev=r;
        thetaPrev=thetaCur;
        
        if(rChanged) {
            double temp=a/r;
            crScal=(c/r)/scalFactor;

            for(n=1;n<=M;n++) {
                nCoeff[n]=nCoeff[n-1]*temp;
            }
        }
        
        if(algorithm==0) {
        //At latitudes that are not near the poles, the Legendre method is
        //used. It cannot be used for the gradient or Hessian near the
        //poles, because of the singularity of the spherical coordinate
        //system.
            useLegendre=(fabs(thetaCur)<88.0*pi/180.0||(gradV==NULL&&HessianV==NULL));
        } else if (algorithm==1) {
            useLegendre=true;
        } else {
            useLegendre=false;
        }

        if(useLegendre) {
            double u, theta;
            
            //Compute the sine and cosine terms 
            //Explicitly set the first two terms.
            SinVec[0]=0;
            CosVec[0]=1;
            SinVec[1]=sin(lambda);
            CosVec[1]=cos(lambda);
            //Use a double angle identity to get the second order term.
            SinVec[2]=2.0*SinVec[1]*CosVec[1];
            CosVec[2]=1.0-2.0*SinVec[1]*SinVec[1];
            //Use a two-part recursion for the rest of the terms.
            for(m=3;m<=M;m++){
                SinVec[m]=2.0*CosVec[1]*SinVec[m-1]-SinVec[m-2];
                CosVec[m]=2.0*CosVec[1]*CosVec[m-1]-CosVec[m-2];
            }
                        
            //The formulae for spherical harmonic synthesis with legendre's
            //method uses clolatitude, pi/2-elevation
            theta=pi/2-thetaCur;
            u=sin(theta);
            if(thetaChanged) {
                //Get the associated Legendre function ratios.
                NALegendreCosRatCPP(FuncVals,  theta, scalFactor);

                //Get Legendre ratio values if the gradient or Hessian
                //terms are desired.
                if(gradV!=NULL||HessianV!=NULL) {
                    NALegendreCosRatDerivCPP(FuncDerivs,FuncVals,theta);
                    
                    if(HessianV!=NULL) {
                        NALegendreCosRatDeriv2CPP(FuncDerivs2,FuncDerivs,FuncVals,theta);
                    }
                }
            }
            
            //Evaluate Equation 7 from the Holmes and Featherstone paper.
            if(rChanged||thetaChanged) {
                //Zero the arrays.
                memset(A,0,sizeof(double)*M1);
                memset(B,0,sizeof(double)*M1);

                //Compute the X coefficients for the sum
                for(m=0;m<=M;m++) {
                    for(n=m;n<=M;n++) {
                        A[m]+=nCoeff[n]*C[n][m]*FuncVals[n][m];
                        B[m]+=nCoeff[n]*S[n][m]*FuncVals[n][m];
                    }
                }
                
                //If additional terms should be computed so a gradient or
                //Hessian can be computed.
                if(gradV!=NULL||HessianV!=NULL) {
                    memset(Ar,0,sizeof(double)*M1);
                    memset(Br,0,sizeof(double)*M1);
                    memset(ATheta,0,sizeof(double)*M1);
                    memset(BTheta,0,sizeof(double)*M1);

                    mf=0;
                    for(m=0;m<=M;m++) {
                        nf=mf;
                        for(n=m;n<=M;n++) {
                            const double CScal=nCoeff[n]*C[n][m];
                            const double SScal=nCoeff[n]*S[n][m];
                            double curVal;
                            
                            curVal=FuncVals[n][m];
                            Ar[m]+=(nf+1)*CScal*curVal;
                            Br[m]+=(nf+1)*SScal*curVal;
                            
                            curVal=FuncDerivs[n][m];
                            ATheta[m]+=CScal*curVal;
                            BTheta[m]+=SScal*curVal;
                                                        
                            nf++;
                        }
                                    
                        mf++;
                    }
                    
                    if(HessianV!=NULL) {
                        memset(Arr,0,sizeof(double)*M1);
                        memset(Brr,0,sizeof(double)*M1);
                        memset(AThetar,0,sizeof(double)*M1);
                        memset(BThetar,0,sizeof(double)*M1);
                        memset(AThetaTheta,0,sizeof(double)*M1);
                        memset(BThetaTheta,0,sizeof(double)*M1);
                        
                        mf=0;
                        for(m=0;m<=M;m++) {
                            nf=mf;
                            for(n=m;n<=M;n++) {
                                const double CScal=nCoeff[n]*C[n][m];
                                const double SScal=nCoeff[n]*S[n][m];
                                double curVal;
                                
                                curVal=FuncVals[n][m];
                                //From Table 5, with the correction from the
                                //erratum.
                                Arr[m]+=(nf+1)*(nf+2)*CScal*curVal;
                                Brr[m]+=(nf+1)*(nf+2)*SScal*curVal;;
                                
                                curVal=FuncDerivs[n][m];
                                //From Table 5
                                AThetar[m]+=(nf+1)*CScal*curVal;
                                BThetar[m]+=(nf+1)*SScal*curVal;
                                
                                curVal=FuncDerivs2[n][m];
                                //From Table 5
                                AThetaTheta[m]+=CScal*curVal;
                                BThetaTheta[m]+=SScal*curVal;

                                nf++;
                            }
                            
                            mf++;
                        }
                    }
                }
            }
 
            //Use Horner's method to compute V.
            V[curPoint]=0;
            m=M+1;
            do {
                m--;
                
                V[curPoint]=V[curPoint]*u+A[m]*CosVec[m]+B[m]*SinVec[m];
            } while(m>0);

            //Multiply by the constant in front of the sum and get rid of the
            //scale factor.
            V[curPoint]=crScal*V[curPoint];

            //Compute the gradient, if it is desired.
            if(gradV!=NULL||HessianV!=NULL) {
                double J[9];
                double dVdr=0;
                double dVdLambda=0;
                double dVdTheta=0;

                //Use Horner's method to compute all dV values.
                m=M+1;
                mf=static_cast<double>(m);
                do {
                    m--;
                    mf--;
                    
                    dVdr=dVdr*u-(Ar[m]*CosVec[m]+Br[m]*SinVec[m]);
                    dVdLambda=dVdLambda*u-mf*(A[m]*SinVec[m]-B[m]*CosVec[m]);
                    dVdTheta=dVdTheta*u+(ATheta[m]*CosVec[m]+BTheta[m]*SinVec[m]);
                } while(m>0);            

                dVdr=crScal*dVdr/r;
                dVdLambda=crScal*dVdLambda;
                //The minus sign adjusts for the coordinate system change.
                dVdTheta=-crScal*dVdTheta;
                
                if(spherDerivs&&gradV!=NULL) {
                    const size_t idx=3*curPoint;
                    gradV[0+idx]=dVdr;
                    gradV[1+idx]=dVdLambda;
                    gradV[2+idx]=dVdTheta;
                } else {
                    const size_t idx=3*curPoint;
                    calcSpherConvJacobCPP(J, pointCur,0);

                    if(gradV!=NULL) {
                        //Multiply the transpose of the Jacobian Matrix by
                        //the vector of [dVdr;dVdLambda;dVdTheta]
                        gradV[0+idx]=dVdr*J[0]+dVdLambda*J[1]+dVdTheta*J[2];
                        gradV[1+idx]=dVdr*J[3]+dVdLambda*J[4]+dVdTheta*J[5];
                        gradV[2+idx]=dVdr*J[6]+dVdLambda*J[7]+dVdTheta*J[8];
                    }
                }

                if(HessianV!=NULL) {
                     
                    //Use Horner's method to compute d2V all values.
                    double d2VdLambdadLambda=0;
                    double d2VdLambdadTheta=0;
                    double d2VdrdLambda=0;
                    double d2VdThetadTheta=0;
                    double d2VdrdTheta=0;
                    double d2Vdrdr=0;
                    
                    //The following second-order derivative formulae are
                    //from Table 2 (expressed using Horner's method).
                    m=M1;
                    mf=static_cast<double>(m);
                    do {
                        m--;
                        mf--;
                        d2VdLambdadLambda=d2VdLambdadLambda*u-mf*mf*(A[m]*CosVec[m]+B[m]*SinVec[m]);
                        d2VdLambdadTheta=d2VdLambdadTheta*u+mf*(ATheta[m]*SinVec[m]-BTheta[m]*CosVec[m]);
                        d2VdrdLambda=d2VdrdLambda*u+mf*(Ar[m]*SinVec[m]-Br[m]*CosVec[m]);
                        d2VdThetadTheta=d2VdThetadTheta*u+(AThetaTheta[m]*CosVec[m]+BThetaTheta[m]*SinVec[m]);
                        d2VdrdTheta=d2VdrdTheta*u-(AThetar[m]*CosVec[m]+BThetar[m]*SinVec[m]);
                        d2Vdrdr=d2Vdrdr*u+(Arr[m]*CosVec[m]+Brr[m]*SinVec[m]);
                    } while(m>0);

                    //The minus signs in the following equations adjust for the
                    //spherical coordinate system difference.
                    d2Vdrdr=crScal*d2Vdrdr/(r*r);
                    d2VdLambdadLambda=crScal*d2VdLambdadLambda;
                    d2VdThetadTheta=crScal*d2VdThetadTheta;
                    d2VdrdLambda=crScal*d2VdrdLambda/r;
                    d2VdrdTheta=-crScal*d2VdrdTheta/r;
                    d2VdLambdadTheta=crScal*d2VdLambdadTheta;

                    if(spherDerivs) {
                        const size_t idx=9*curPoint;
                        
                        HessianV[0+idx]=d2Vdrdr;
                        HessianV[1+idx]=d2VdrdLambda;
                        HessianV[2+idx]=d2VdrdTheta;
                        
                        HessianV[3+idx]=d2VdrdLambda;
                        HessianV[4+idx]=d2VdLambdadLambda;
                        HessianV[5+idx]=d2VdLambdadTheta;
                        
                        HessianV[6+idx]=d2VdrdTheta;
                        HessianV[7+idx]=d2VdLambdadTheta;
                        HessianV[8+idx]=d2VdThetadTheta;

                    } else {
                        double H[27];
                        const double drdx=J[0];
                        const double drdy=J[3];
                        const double drdz=J[6];
                        const double dLambdadx=J[1];
                        const double dLambdady=J[4];
                        const double dLambdadz=J[7];
                        const double dPhidx=J[2];
                        const double dPhidy=J[5];
                        const double dPhidz=J[8];
                        const size_t idx=9*curPoint;

                        calcSpherConvHessianCPP(H,pointCur,0,true);

                        //Commented indices on H are what would be used in
                        //Matlab.
                        const double drdxdx=H[0];//H(1,1,1)
                        const double drdydy=H[4];//H(2,2,1)
                        const double drdzdz=H[8];//H(3,3,1)
                        const double drdxdy=H[3];//H(1,2,1)
                        const double drdxdz=H[6];//H(1,3,1)
                        const double drdydz=H[7];//H(2,3,1)

                        const double dLambdadxdx=H[9];//H(1,1,2)
                        const double dLambdadydy=H[13];//H(2,2,2)
                        const double dLambdadzdz=H[17];//H(3,3,2)
                        const double dLambdadxdy=H[12];//H(1,2,2)
                        const double dLambdadxdz=H[15];//H(1,3,2)
                        const double dLambdadydz=H[16];//H(2,3,2)

                        const double dPhidxdx=H[18];//H(1,1,3)
                        const double dPhidydy=H[22];//H(2,2,3)
                        const double dPhidzdz=H[26];//H(3,3,3)
                        const double dPhidxdy=H[21];//H(1,2,3)
                        const double dPhidxdz=H[24];//H(1,3,3)
                        const double dPhidydz=H[25];//H(2,3,3)

                        //HessianV(1,1,curPoint)
                        HessianV[0+idx]=d2VdLambdadLambda*(dLambdadx*dLambdadx)+2.0*d2VdLambdadTheta*dLambdadx*dPhidx+d2VdThetadTheta*(dPhidx*dPhidx)+2.0*dPhidx*drdx*d2VdrdTheta+2.0*dLambdadx*drdx*d2VdrdLambda+dVdLambda*dLambdadxdx+dVdTheta*dPhidxdx+dVdr*drdxdx+(drdx*drdx)*d2Vdrdr;
                        //HessianV(2,2,curPoint)
                        HessianV[4+idx]=d2VdThetadTheta*(dPhidy*dPhidy)+2.0*dLambdady*dPhidy*d2VdLambdadTheta+dVdLambda*dLambdadydy+dVdTheta*dPhidydy+(dLambdady*dLambdady)*d2VdLambdadLambda+drdydy*dVdr+2.0*dPhidy*drdy*d2VdrdTheta+2.0*dLambdady*drdy*d2VdrdLambda+(drdy*drdy)*d2Vdrdr;
                        //HessianV(3,3,curPoint)
                        HessianV[8+idx]=dVdTheta*dPhidzdz+(dPhidz*dPhidz)*d2VdThetadTheta+dLambdadzdz*dVdLambda+2.0*dLambdadz*dPhidz*d2VdLambdadTheta+(dLambdadz*dLambdadz)*d2VdLambdadLambda+drdzdz*dVdr+2.0*dPhidz*drdz*d2VdrdTheta+2.0*dLambdadz*drdz*d2VdrdLambda+(drdz*drdz)*d2Vdrdr;
                        //HessianV(1,2,curPoint)
                        HessianV[3+idx]=dPhidy*d2VdLambdadTheta*dLambdadx+dLambdady*d2VdLambdadLambda*dLambdadx+d2VdThetadTheta*dPhidy*dPhidx+dLambdady*d2VdLambdadTheta*dPhidx+drdy*dPhidx*d2VdrdTheta+dPhidy*drdx*d2VdrdTheta+dVdLambda*dLambdadxdy+dVdTheta*dPhidxdy+dVdr*drdxdy+drdy*dLambdadx*d2VdrdLambda+dLambdady*drdx*d2VdrdLambda+drdy*drdx*d2Vdrdr;
                        //HessianV(2,1,curPoint)
                        HessianV[1+idx]=HessianV[3+idx];
                        //HessianV(1,3,curPoint)
                        HessianV[6+idx]=dPhidz*d2VdLambdadTheta*dLambdadx+dLambdadz*d2VdLambdadLambda*dLambdadx+dPhidz*d2VdThetadTheta*dPhidx+dLambdadz*d2VdLambdadTheta*dPhidx+dVdLambda*dLambdadxdz+dVdTheta*dPhidxdz+dVdr*drdxdz+drdz*dPhidx*d2VdrdTheta+dPhidz*drdx*d2VdrdTheta+drdz*dLambdadx*d2VdrdLambda+dLambdadz*drdx*d2VdrdLambda+drdz*drdx*d2Vdrdr;
                        //HessianV(3,1,curPoint)
                        HessianV[2+idx]=HessianV[6+idx];
                        //HessianV(2,3,curPoint)
                        HessianV[7+idx]=dPhidz*d2VdThetadTheta*dPhidy+dVdLambda*dLambdadydz+dVdTheta*dPhidydz+dPhidz*dLambdady*d2VdLambdadTheta+dLambdadz*dPhidy*d2VdLambdadTheta+dLambdadz*dLambdady*d2VdLambdadLambda+drdydz*dVdr+drdz*dPhidy*d2VdrdTheta+dPhidz*drdy*d2VdrdTheta+drdz*dLambdady*d2VdrdLambda+dLambdadz*drdy*d2VdrdLambda+drdz*drdy*d2Vdrdr;
                        //HessianV(3,2,curPoint)
                        HessianV[5+idx]=HessianV[7+idx];
                    }
                }
            }
        } else {  
        //At latitudes that are near the poles, the non-singular algorithm
        //of Pines using the fully normalized Helmholtz equations from
        //Fantino and Casotto is used. The algorithm has been slightly
        //modified so that the c/r term is out front and the fully
        //normalized Helmholtz polynomials can be scaled. Also, lumped
        //coefficients are not used. The Pines algorithm is generally
        //slower than the algorithm of Holmes and Featherstone and it
        //suffers a loss of precision near the equator. Thus, the Pines
        //algorithm is only used near the poles where the other algorithm
        //has issues with a singularity.
            double CartPoint[3];

            spher2CartCPP(CartPoint,pointCur,0);

            //Get the direction cosines used by Pines' algorithm.
            const double s=CartPoint[0]/r;
            const double t=CartPoint[1]/r;
            const double u=CartPoint[2]/r;

            //Compute the fully normalized Helmholtz polynomials.
            if(thetaChanged) {
                normHelmHoltzCPP(FuncVals,u, scalFactor);
                
                if(gradV!=NULL||HessianV!=NULL) {
                    normHelmHoltzDerivCPP(FuncDerivs,FuncVals);
                }
                
                if(HessianV!=NULL) {
                    normHelmHoltzDeriv2CPP(FuncDerivs2,FuncVals);   
                }
                
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

            V[curPoint]=crScal*V[curPoint];
            
            //If only the gradient is desired as the next output       
            if(gradV!=NULL&&HessianV==NULL) {
                double a1=0;
                double a2=0;
                double a3=0;
                double a4=0;

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    double A1=0;
                    double A2=0;
                    double A3=0;
                    double B1=0;
                    double B2=0;
                    double B3=0;
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        //The expressions for Lmn, is from Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;       
                        const double rhoC=nCoeff[n]*C[n][m];
                        const double rhoS=nCoeff[n]*S[n][m];

                        A1+=rhoC*HVal;
                        A2+=rhoC*dHVal;
                        A3+=rhoC*Lmn;

                        B1+=rhoS*HVal;
                        B2+=rhoS*dHVal;
                        B3+=rhoS*Lmn;
                        
                        nf++;
                    }
                    if(M>=1) {
                        a1+=mf*(A1*rm[m-1]+B1*im[m-1]);
                        a2+=mf*(B1*rm[m-1]-A1*im[m-1]);
                    }
                    a3+=(A2*rm[m]+B2*im[m]);
                    a4-=(A3*rm[m]+B3*im[m]);

                    mf++;
                }
                a1/=r;
                a2/=r;
                a3/=r;
                a4/=r;
                
                {
                    const double dVdx=crScal*(a1+s*a4);
                    const double dVdy=crScal*(a2+t*a4);
                    const double dVdz=crScal*(a3+u*a4);

                    if(spherDerivs) {
                        double J[9];
                        size_t idx=3*curPoint;
                        //Convert the derivatives to spherical coordinates.
                        calcSpherInvJacobCPP(J,pointCur,0);

                        //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                        gradV[0+idx]=dVdx*J[0]+dVdy*J[1]+dVdz*J[2];
                        gradV[1+idx]=dVdx*J[3]+dVdy*J[4]+dVdz*J[5];
                        gradV[2+idx]=dVdx*J[6]+dVdy*J[7]+dVdz*J[8];
                    } else {//If a gradient in Cartesian coordinates is desired.
                        size_t idx=3*curPoint;

                        gradV[0+idx]=dVdx;
                        gradV[1+idx]=dVdy;
                        gradV[2+idx]=dVdz;
                    }
                }
            } else if(HessianV!=NULL) {
            //Compute the gradient and Hessian
                double a1=0;
                double a2=0;
                double a3=0;
                double a4=0;
                double a11=0;
                double a12=0;
                double a13=0;
                double a14=0;
                double a23=0;
                double a24=0;
                double a33=0;
                double a34=0;
                double a44=0;
                const double r2=r*r;
                const double s2=s*s;
                const double u2=u*u;
                const double t2=t*t;
                double a22;

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    double A1=0;
                    double A2=0;
                    double A3=0;
                    double A4=0;
                    double A5=0;
                    double A6=0;
                    double B1=0;
                    double B2=0;
                    double B3=0;
                    double B4=0;
                    double B5=0;
                    double B6=0;
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        const double d2HVal=FuncDerivs2[n][m];

                        //The expressions for Lmn, dLmn, and Omn are from
                        //Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;
                        const double dLmn=(nf+mf+2)*dHVal+u*d2HVal;
                        const double Omn=(nf+mf+1)*(n+m+2)*HVal+2.0*u*(nf+mf+2)*dHVal+u*u*d2HVal;

                        const double rhoC=nCoeff[n]*C[n][m];
                        const double rhoS=nCoeff[n]*S[n][m];

                        A1+=rhoC*HVal;
                        A2+=rhoC*dHVal;
                        A3+=rhoC*Lmn;
                        A4+=rhoC*d2HVal;
                        A5+=rhoC*dLmn;
                        A6+=rhoC*Omn;

                        B1+=rhoS*HVal;
                        B2+=rhoS*dHVal;
                        B3+=rhoS*Lmn;
                        B4+=rhoS*d2HVal;
                        B5+=rhoS*dLmn;
                        B6+=rhoS*Omn;
                        
                        nf++;
                    }

                    if(m>=1) {
                        a1+=mf*(A1*rm[m-1]+B1*im[m-1]);
                        a2+=mf*(B1*rm[m-1]-A1*im[m-1]);
                        
                        a13+=mf*(A2*rm[m-1]+B2*im[m-1]);
                        a14-=mf*(A3*rm[m-1]+B3*im[m-1]);
                        a23+=mf*(B2*rm[m-1]-A2*im[m-1]);
                        a24-=mf*(B3*rm[m-1]-A3*im[m-1]);
                        
                        if(m>=2) {
                            a11+=mf*(mf-1)*(A1*rm[m-2]+B1*im[m-2]);
                            a12+=mf*(mf-1)*(B1*rm[m-2]-A1*im[m-2]);
                        }
                    }
                    a3+=(A2*rm[m]+B2*im[m]);
                    a4-=(A3*rm[m]+B3*im[m]);

                    a33+=(A4*rm[m]+B4*im[m]);
                    a34-=(A5*rm[m]+B5*im[m]);
                    a44+=(A6*rm[m]+B6*im[m]);

                    mf++;
                }
                a1/=r;
                a2/=r;
                a3/=r;
                a4/=r;
                a11/=r2;
                a12/=r2;
                a13/=r2;
                a14/=r2;
                a23/=r2;
                a24/=r2;
                a33/=r2;
                a34/=r2;
                a44/=r2;
                a22=-a11;

                {
                    double dxdr;
                    double dxdAz;
                    double dxdEl;
                    double dydr;
                    double dydAz;
                    double dydEl;
                    double dzdr;
                    double dzdAz;
                    double dzdEl;
                    const double dVdx=crScal*(a1+s*a4);
                    const double dVdy=crScal*(a2+t*a4);
                    const double dVdz=crScal*(a3+u*a4);
                    const double d2Vdxdx=crScal*(a11+2.0*s*a14+a4/r+s2*a44-s2*a4/r);
                    const double d2Vdydy=crScal*(a22+2.0*t*a24+a4/r+t2*a44-t2*a4/r);
                    const double d2Vdzdz=crScal*(a33+2.0*u*a34+a4/r+u2*a44-u2*a4/r);
                    const double d2Vdxdy=crScal*(a12+s*t*a44+s*a24+t*a14-s*t*a4/r);
                    const double d2Vdxdz=crScal*(a13+s*u*a44+s*a34+u*a14-s*u*a4/r);
                    const double d2Vdydz=crScal*(a23+t*u*a44+t*a34+u*a24-t*u*a4/r);

                    if(spherDerivs) {
                        double J[9];
                        calcSpherInvJacobCPP(J,pointCur,0);

                        dxdr=J[0];//J(1,1);
                        dxdAz=J[3];//J(1,2);
                        dxdEl=J[6];//J(1,3);

                        dydr=J[1];//J(2,1);
                        dydAz=J[4];//J(2,2);
                        dydEl=J[7];//J(2,3);

                        dzdr=J[2];//J(3,1);
                        dzdAz=J[5];//J(3,2);
                        dzdEl=J[8];//J(3,3);
                    }
                    
                    //Compute the gradient, if desired.
                    if(gradV!=NULL) {
                        if(spherDerivs) {
                            const size_t idx=3*curPoint;
                            
                            //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                            gradV[0+idx]=dVdx*dxdr+dVdy*dydr+dVdz*dzdr;
                            gradV[1+idx]=dVdx*dxdAz+dVdy*dydAz+dVdz*dzdAz;
                            gradV[2+idx]=dVdx*dxdEl+dVdy*dydEl+dVdz*dzdEl;
                        } else {
                            const size_t idx=3*curPoint;

                            gradV[0+idx]=dVdx;
                            gradV[1+idx]=dVdy;
                            gradV[2+idx]=dVdz;
                        }
                    }
                    
                    //Compute the Hessian.
                    if(spherDerivs) {
                        const size_t idx=9*curPoint;
                        double H[27];
   
                        calcSpherInvHessianCPP(H,pointCur,0);

                        const double d2xdrdr=H[0];//H(1,1,1);
                        const double d2xdAzdAz=H[4];//H(2,2,1);
                        const double d2xdEldEl=H[8];//H(3,3,1);
                        const double d2xdrdAz=H[3];//H(1,2,1);
                        const double d2xdrdEl=H[6];//H(1,3,1);
                        const double d2xdAzdEl=H[7];//H(2,3,1);

                        const double d2ydrdr=H[9];//H(1,1,2);
                        const double d2ydAzdAz=H[13];//H(2,2,2);
                        const double d2ydEldEl=H[17];//H(3,3,2);
                        const double d2ydrdAz=H[12];//H(1,2,2);
                        const double d2ydrdEl=H[15];//H(1,3,2);
                        const double d2ydAzdEl=H[16];//H(2,3,2);

                        const double d2zdrdr=H[18];//H(1,1,3);
                        const double d2zdAzdAz=H[22];//H(2,2,3);
                        const double d2zdEldEl=H[26];//H(3,3,3);
                        const double d2zdrdAz=H[21];//H(1,2,3);
                        const double d2zdrdEl=H[24];//H(1,3,3);
                        const double d2zdAzdEl=H[25];//H(2,3,3);

                        //d2Vdrdr, HessianV(1,1,curPoint)
                        HessianV[0+idx]=d2Vdydy*dydr*dydr+2.0*d2Vdydz*dydr*dzdr+d2Vdzdz*dzdr*dzdr+2.0*dxdr*dzdr*d2Vdxdz+2.0*dxdr*dydr*d2Vdxdy+dxdr*dxdr*d2Vdxdx+dVdx*d2xdrdr+dVdy*d2ydrdr+dVdz*d2zdrdr;
                        //d2VdAzdAz, HessianV(2,2,curPoint)
                        HessianV[4+idx]=d2Vdzdz*dzdAz*dzdAz+2.0*dydAz*dzdAz*d2Vdydz+dydAz*dydAz*d2Vdydy+dVdy*d2ydAzdAz+dVdz*d2zdAzdAz+d2xdAzdAz*dVdx+2.0*dxdAz*dzdAz*d2Vdxdz+2.0*dxdAz*dydAz*d2Vdxdy+dxdAz*dxdAz*d2Vdxdx;
                        //d2VdEldEl, HessianV(3,3,curPoint)
                        HessianV[8+idx]=dzdEl*dzdEl*d2Vdzdz+dVdz*d2zdEldEl+d2ydEldEl*dVdy+2.0*dydEl*dzdEl*d2Vdydz+dydEl*dydEl*d2Vdydy+d2xdEldEl*dVdx+2.0*dxdEl*dzdEl*d2Vdxdz+2.0*dxdEl*dydEl*d2Vdxdy+dxdEl*dxdEl*d2Vdxdx;
                        //d2VdrdAz, HessianV(1,2,curPoint)
                        HessianV[3+idx]=dzdAz*d2Vdydz*dydr+dydAz*d2Vdydy*dydr+d2Vdzdz*dzdAz*dzdr+dydAz*d2Vdydz*dzdr+dzdAz*dxdr*d2Vdxdz+dxdAz*dzdr*d2Vdxdz+dydAz*dxdr*d2Vdxdy+dxdAz*dydr*d2Vdxdy+dVdx*d2xdrdAz+dVdy*d2ydrdAz+dVdz*d2zdrdAz+dxdAz*dxdr*d2Vdxdx;
                        //d2VdAzdr, HessianV(2,1,curPoint)
                        HessianV[1+idx]=HessianV[3+idx];
                        //d2VdrdEl, HessianV(1,3,curPoint)
                        HessianV[6+idx]=dzdEl*d2Vdydz*dydr+dydEl*d2Vdydy*dydr+dzdEl*d2Vdzdz*dzdr+dydEl*d2Vdydz*dzdr+dzdEl*dxdr*d2Vdxdz+dxdEl*dzdr*d2Vdxdz+dVdx*d2xdrdEl+dVdy*d2ydrdEl+dVdz*d2zdrdEl+dydEl*dxdr*d2Vdxdy+dxdEl*dydr*d2Vdxdy+dxdEl*dxdr*d2Vdxdx;
                        //d2VdEldr, HessianV(3,1,curPoint)
                        HessianV[2+idx]=HessianV[6+idx];
                        //d2VdAzdEl, HessianV(2,3,curPoint)
                        HessianV[7+idx]=dzdEl*d2Vdzdz*dzdAz+dzdEl*dydAz*d2Vdydz+dydEl*dzdAz*d2Vdydz+dVdy*d2ydAzdEl+dVdz*d2zdAzdEl+dydEl*dydAz*d2Vdydy+d2xdAzdEl*dVdx+dzdEl*dxdAz*d2Vdxdz+dxdEl*dzdAz*d2Vdxdz+dydEl*dxdAz*d2Vdxdy+dxdEl*dydAz*d2Vdxdy+dxdEl*dxdAz*d2Vdxdx;
                        //d2VdEldAz, HessianV(3,2,curPoint)
                        HessianV[5+idx]=HessianV[7+idx];
                    } else {
                        const size_t idx=9*curPoint;
                        
                        //HessianV(1,1,curPoint)
                        HessianV[0+idx]=d2Vdxdx;
                        //HessianV(2,2,curPoint)
                        HessianV[4+idx]=d2Vdydy;
                        //HessianV(3,3,curPoint)
                        HessianV[8+idx]=d2Vdzdz;
                        //HessianV(1,2,curPoint)
                        HessianV[3+idx]=d2Vdxdy;
                        //HessianV(2,1,curPoint)
                        HessianV[1+idx]=HessianV[3+idx];
                        //HessianV(1,3,curPoint)
                        HessianV[6+idx]=d2Vdxdz;
                        //HessianV(3,1,curPoint)
                        HessianV[2+idx]=HessianV[6+idx];
                        //HessianV(2,3,curPoint)
                        HessianV[7+idx]=d2Vdydz;
                        //HessianV(3,2,curPoint)
                        HessianV[5+idx]=HessianV[7+idx];
                    }
                }
            }
        }
        
        if(spherDerivs&&systemType==2) {
            size_t idx;
            //Flip signs of the elevation terms reflecting the
            //difference in definition between systems 0 and 2.
            
            if(gradV!=NULL) {
                idx=3*curPoint;
                gradV[2+idx]=-gradV[2+idx];
            }

            if(HessianV!=NULL) {
                idx=9*curPoint;
                //HessianV(1,3,curPoint)
                HessianV[6+idx]=-HessianV[6+idx];
                //HessianV(3,1,curPoint)
                HessianV[2+idx]=-HessianV[2+idx];
                //HessianV(2,3,curPoint)
                HessianV[7+idx]=-HessianV[7+idx];
                //HessianV(3,2,curPoint)
                HessianV[5+idx]=-HessianV[5+idx];
            }
        }
    }
    
    delete[] buffer;
}

void spherHarmonicEvalCPPComplex(complex<double> *VRet, complex<double> *gradVRet, complex<double> *HessianVRet,const CountingClusterSetCPP<complex<double>> &C,const CountingClusterSetCPP<complex<double>> &S, const double *point, const size_t numPoints, const complex <double> a, const complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm) {
    //If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. If a NULL pointer is passed for HessianV,
    //then it is assumed that the Hessian is not desired.
    const size_t M=C.numClust-1;
    const size_t M1=C.numClust;
    const double pi=2.0*acos(0.0);
    size_t curPoint;
    double rPrev, thetaPrev;
    complex <double> crScal;
    CountingClusterSetCPP<double> FuncVals;
    CountingClusterSetCPP<double> FuncDerivs;
    CountingClusterSetCPP<double> FuncDerivs2;
    //A big chunk of memory will be allocated into two buffers, one for
    //real and one for complex variables, and split between the variables
    //that need it. That is faster than allocating a bunch of small
    //buffers.
    double *buffer;
    complex <double> *bufferComplex;
    //This stores all of the powers of a/r needed for the sum, regardless
    //of which algorithm is used.
    complex <double> *nCoeff;//Length M+1
    //These are to store sin(m*lambda) and cos(m*lambda) for m=0->M for the
    //Legendre algorithm.
    double *SinVec,*CosVec;//These are each length M+1.
    //rm and im are never used at the same time as SinVec andCosVec,
    //because rm and im are used by Pines' algorithm. They will be the same
    //size as SinVec and CosVec, and thus will point to the same memory.
    double *rm,*im;
    //The following are needed if the Legendre algorithm is used.
    complex<double> *A=NULL;//Length M1
    complex<double> *B=NULL;//Length M1
    //The following are only used if the gradient or Hessian is desired and
    //the Legendre algorithm is used.
    complex<double> *Ar=NULL;//Length M1
    complex<double> *Br=NULL;//Length M1
    complex<double> *ATheta=NULL;//Length M1
    complex<double> *BTheta=NULL;//Length M1

    //The following are only used if the Hessian is desired and the
    //Legendre algorithm is used.
    complex<double> *AThetaTheta=NULL;//Length M1
    complex<double> *BThetaTheta=NULL;//Length M1
    complex<double> *Arr=NULL;//Length M1
    complex<double> *Brr=NULL;//Length M1
    complex<double> *AThetar=NULL;//Length M1
    complex<double> *BThetar;//Length M1

    //Initialize the ClusterSet classes for the coefficients. The space
    //for the elements will be allocated shortly.
    //FuncVals will be used for HBar and PBarUVals.
    //FuncDerivs will be used for dPBarUValsdTheta and dHBardu.
    //FuncDerivs2 will be used for d2PBarUValsdTheta2 and d2HBardu2
    FuncVals.numClust=C.numClust;
    FuncVals.totalNumEl=C.totalNumEl;
    FuncDerivs.numClust=C.numClust;
    FuncDerivs.totalNumEl=C.totalNumEl;
    FuncDerivs2.numClust=C.numClust;
    FuncDerivs2.totalNumEl=C.totalNumEl;
    
    if(algorithm==2) {//If only Pines' algorithm is used.
        //This pointer is used to help parition the memory in the real
        //buffer.
        double *tempPtr;

        if(gradVRet==NULL&&HessianVRet==NULL) {
            buffer=new double[2*M1+C.totalNumEl];
            bufferComplex=new complex <double>[M1];
        } else{
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex <double>[M1];
        }
        
        //The complex buffer
        nCoeff=bufferComplex;//Length M1

        //The real buffers
        tempPtr=buffer;
        rm=tempPtr;//Length M1

        tempPtr+=M1;
        im=tempPtr;//Length M1

        tempPtr+=M1;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        if(!(gradVRet==NULL&&HessianVRet==NULL)) {
            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl
            
            tempPtr+=C.totalNumEl;
            FuncDerivs2.clusterEls=tempPtr;// Length C.totalNumEl
        }
    } else {
        double *tempPtr;
        complex<double> *tempPtrComplex;
        
        //If the Legendre algorithm might be used, then the amount of
        //memory needed varies depending on the highest order derivative
        //that is needed.
        if(HessianVRet!=NULL) {
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex <double>[13*M1];
        } else if(gradVRet!=NULL) {
            if(algorithm==0) {//If Pines' algorithm might be used
                buffer=new double[2*M1+3*C.totalNumEl];
                bufferComplex=new complex <double>[7*M1];
            } else {
                buffer=new double[2*M1+2*C.totalNumEl];
                bufferComplex=new complex <double>[7*M1];
            }
        } else {
            buffer=new double[2*M1+C.totalNumEl];
            bufferComplex=new complex <double>[3*M1];
        }
        
        //The complex buffer
        tempPtrComplex=bufferComplex;
        nCoeff=bufferComplex;//Length M1
        
        tempPtrComplex+=M1;
        A=tempPtrComplex;//Length M1

        tempPtrComplex+=M1;
        B=tempPtrComplex;//Length M1
        
        tempPtr=buffer;
        rm=tempPtr;
        SinVec=tempPtr;//Length M1

        tempPtr+=M1;
        im=tempPtr;
        CosVec=tempPtr;//Length M1

        tempPtr+=M1;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
            
        if(gradVRet!=NULL||HessianVRet!=NULL) {
            tempPtrComplex+=M1;
            Ar=tempPtrComplex;//Length M1

            tempPtrComplex+=M1;
            Br=tempPtrComplex;//Length M1

            tempPtrComplex+=M1;
            ATheta=tempPtrComplex;//Length M1

            tempPtrComplex+=M1;
            BTheta=tempPtrComplex;//Length M1

            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl

            if(HessianVRet!=NULL||algorithm!=0) {
                tempPtr+=C.totalNumEl;
                FuncDerivs2.clusterEls=tempPtr;//Length C.totalNumEl
            }

            if(HessianVRet!=NULL) {
                tempPtrComplex+=M1;
                AThetaTheta=tempPtrComplex;//Length M1

                tempPtrComplex+=M1;
                BThetaTheta=tempPtrComplex;//Length M1

                tempPtrComplex+=M1;
                Arr=tempPtrComplex;//Length M1

                tempPtrComplex+=M1;
                Brr=tempPtrComplex;//Length M1

                tempPtrComplex+=M1;
                AThetar=tempPtrComplex;//Length M1

                tempPtrComplex+=M1;
                BThetar=tempPtrComplex;//Length M1
            }
        }
    }

    nCoeff[0]=1;
    
    rPrev=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double pointCur[3];
        double r, lambda, thetaCur;
        size_t n,m;
        double nf,mf;
        bool useLegendre;
        bool rChanged;
        bool thetaChanged;
        
        pointCur[0]=point[0+3*curPoint];
        pointCur[1]=point[1+3*curPoint];
        if(systemType==2) {
            pointCur[2]=pi/2.0-point[2+3*curPoint];
        } else {
            pointCur[2]=point[2+3*curPoint];
        }

        r=pointCur[0];
        lambda=pointCur[1];
        thetaCur=pointCur[2];
        
        rChanged=(rPrev!=r);
        thetaChanged=(thetaCur!=thetaPrev);
        rPrev=r;
        thetaPrev=thetaCur;
        
        if(rChanged) {
            complex <double> temp;
            crScal=(c/r)/scalFactor;

            temp=a/r;
            for(n=1;n<=M;n++) {
                nCoeff[n]=nCoeff[n-1]*temp;
            }
        }

        if(algorithm==0) {
        //At latitudes that are not near the poles, the Legendre method is
        //used. It cannot be used for the gradient or Hessian near the
        //poles, because of the singularity of the spherical coordinate
        //system.
            useLegendre=(fabs(thetaCur)<88.0*pi/180.0||(gradVRet==NULL&&HessianVRet==NULL));
        } else if (algorithm==1) {
            useLegendre=true;
        } else {
            useLegendre=false;
        }

        if(useLegendre) {
            double u, theta;
            
            //Compute the sine and cosine terms 
            //Explicitly set the first two terms.
            SinVec[0]=0;
            CosVec[0]=1;
            SinVec[1]=sin(lambda);
            CosVec[1]=cos(lambda);
            //Use a double angle identity to get the second order term.
            SinVec[2]=2.0*SinVec[1]*CosVec[1];
            CosVec[2]=1.0-2.0*SinVec[1]*SinVec[1];
            //Use a two-part recursion for the rest of the terms.
            for(m=3;m<=M;m++){
                SinVec[m]=2.0*CosVec[1]*SinVec[m-1]-SinVec[m-2];
                CosVec[m]=2.0*CosVec[1]*CosVec[m-1]-CosVec[m-2];
            }
                        
            //The formulae for spherical harmonic synthesis with legendre's
            //method uses clolatitude, pi/2-elevation
            theta=pi/2-thetaCur;
            u=sin(theta);
            if(thetaChanged) {
                //Get the associated Legendre function ratios.
                NALegendreCosRatCPP(FuncVals,  theta, scalFactor);

                //Get Legendre ratio values if the gradient or Hessian
                //terms are desired.
                if(gradVRet!=NULL||HessianVRet!=NULL) {
                    NALegendreCosRatDerivCPP(FuncDerivs,FuncVals,theta);
                    
                    if(HessianVRet!=NULL) {
                        NALegendreCosRatDeriv2CPP(FuncDerivs2,FuncDerivs,FuncVals,theta);
                    }
                }
            }
            
            //Evaluate Equation 7 from the Holmes and Featherstone paper.
            if(rChanged||thetaChanged) {
                //Zero the arrays.
                fill(A,A+M1,complex<double>(0.0,0.0));
                fill(B,B+M1,complex<double>(0.0,0.0));

                //Compute the X coefficients for the sum
                for(m=0;m<=M;m++) {
                    for(n=m;n<=M;n++) {
                        A[m]+=nCoeff[n]*C[n][m]*FuncVals[n][m];
                        B[m]+=nCoeff[n]*S[n][m]*FuncVals[n][m];
                    }
                }
                
                //If additional terms should be computed so a gradient or
                //Hessian can be computed.
                if(gradVRet!=NULL||HessianVRet!=NULL) {
                    fill(Ar,Ar+M1,complex<double>(0.0,0.0));
                    fill(Br,Br+M1,complex<double>(0.0,0.0));
                    fill(ATheta,ATheta+M1,complex<double>(0.0,0.0));
                    fill(BTheta,BTheta+M1,complex<double>(0.0,0.0));

                    mf=0;
                    for(m=0;m<=M;m++) {
                        nf=mf;
                        for(n=m;n<=M;n++) {
                            const complex<double> CScal=nCoeff[n]*C[n][m];
                            const complex<double> SScal=nCoeff[n]*S[n][m];
                            double curVal;
                            
                            curVal=FuncVals[n][m];
                            Ar[m]+=(nf+1)*CScal*curVal;
                            Br[m]+=(nf+1)*SScal*curVal;
                            
                            curVal=FuncDerivs[n][m];
                            ATheta[m]+=CScal*curVal;
                            BTheta[m]+=SScal*curVal;
                                                        
                            nf++;
                        }
                                    
                        mf++;
                    }
                    
                    if(HessianVRet!=NULL) {
                        fill(Arr,Arr+M1,complex<double>(0.0,0.0));
                        fill(Brr,Brr+M1,complex<double>(0.0,0.0));
                        fill(AThetar,AThetar+M1,complex<double>(0.0,0.0));
                        fill(BThetar,BThetar+M1,complex<double>(0.0,0.0));
                        fill(AThetaTheta,AThetaTheta+M1,complex<double>(0.0,0.0));
                        fill(BThetaTheta,BThetaTheta+M1,complex<double>(0.0,0.0));
                        
                        mf=0;
                        for(m=0;m<=M;m++) {
                            nf=mf;
                            for(n=m;n<=M;n++) {
                                const complex<double> CScal=nCoeff[n]*C[n][m];
                                const complex<double> SScal=nCoeff[n]*S[n][m];
                                double curVal;
                                
                                curVal=FuncVals[n][m];
                                //From Table 5, with the correction from the
                                //erratum.
                                Arr[m]+=(nf+1)*(nf+2)*CScal*curVal;
                                Brr[m]+=(nf+1)*(nf+2)*SScal*curVal;;
                                
                                curVal=FuncDerivs[n][m];
                                //From Table 5
                                AThetar[m]+=(nf+1)*CScal*curVal;
                                BThetar[m]+=(nf+1)*SScal*curVal;
                                
                                curVal=FuncDerivs2[n][m];
                                //From Table 5
                                AThetaTheta[m]+=CScal*curVal;
                                BThetaTheta[m]+=SScal*curVal;

                                nf++;
                            }
                            
                            mf++;
                        }
                    }
                }
            }
            
            //Use Horner's method to compute V.
            {
                complex <double> V(0,0);
                m=M+1;
                do {
                    m--;

                    V=V*u+A[m]*CosVec[m]+B[m]*SinVec[m];
                } while(m>0);
                
                //Multiply by the constant in front of the sum and get rid
                //of the scale factor.
                V*=crScal;
                VRet[curPoint]=V;
            }

            //Compute the gradient, if it is desired.
            if(gradVRet!=NULL||HessianVRet!=NULL) {
                double J[9];
                complex<double> dVdr(0,0);
                complex<double> dVdLambda(0,0);
                complex<double> dVdTheta(0,0);

                //Use Horner's method to compute all dV values.
                m=M+1;
                mf=static_cast<double>(m);
                do {
                    m--;
                    mf--;
                    
                    dVdr=dVdr*u-(Ar[m]*CosVec[m]+Br[m]*SinVec[m]);
                    dVdLambda=dVdLambda*u-mf*(A[m]*SinVec[m]-B[m]*CosVec[m]);
                    dVdTheta=dVdTheta*u+(ATheta[m]*CosVec[m]+BTheta[m]*SinVec[m]);
                } while(m>0);            

                dVdr=crScal*dVdr/r;
                dVdLambda=crScal*dVdLambda;
                //The minus sign adjusts for the coordinate system change.
                dVdTheta=-crScal*dVdTheta;
                
                if(spherDerivs&&gradVRet!=NULL) {
                    const size_t idx=3*curPoint;
                    
                    gradVRet[0+idx]=dVdr;
                    gradVRet[1+idx]=dVdLambda;
                    gradVRet[2+idx]=dVdTheta;
                } else {
                    calcSpherConvJacobCPP(J, pointCur,0);

                    if(gradVRet!=NULL) {
                        const size_t idx=3*curPoint;
                        
                        //Multiply the transpose of the Jacobian Matrix by
                        //the vector of [dVdr;dVdLambda;dVdTheta]
                        gradVRet[0+idx]=dVdr*J[0]+dVdLambda*J[1]+dVdTheta*J[2];
                        gradVRet[1+idx]=dVdr*J[3]+dVdLambda*J[4]+dVdTheta*J[5];
                        gradVRet[2+idx]=dVdr*J[6]+dVdLambda*J[7]+dVdTheta*J[8];
                    }
                }

                if(HessianVRet!=NULL) {
                    //Use Horner's method to compute d2V all values.
                    complex<double> d2VdLambdadLambda(0,0);
                    complex<double> d2VdLambdadTheta(0,0);
                    complex<double> d2VdrdLambda(0,0);
                    complex<double> d2VdThetadTheta(0,0);
                    complex<double> d2VdrdTheta(0,0);
                    complex<double> d2Vdrdr(0,0);
                    
                    //The following second-order derivative formulae are
                    //from Table 2 (expressed using Horner's method).
                    m=M1;
                    mf=static_cast<double>(m);
                    do {
                        m--;
                        mf--;
                        d2VdLambdadLambda=d2VdLambdadLambda*u-mf*mf*(A[m]*CosVec[m]+B[m]*SinVec[m]);
                        d2VdLambdadTheta=d2VdLambdadTheta*u+mf*(ATheta[m]*SinVec[m]-BTheta[m]*CosVec[m]);
                        d2VdrdLambda=d2VdrdLambda*u+mf*(Ar[m]*SinVec[m]-Br[m]*CosVec[m]);
                        d2VdThetadTheta=d2VdThetadTheta*u+(AThetaTheta[m]*CosVec[m]+BThetaTheta[m]*SinVec[m]);
                        d2VdrdTheta=d2VdrdTheta*u-(AThetar[m]*CosVec[m]+BThetar[m]*SinVec[m]);
                        d2Vdrdr=d2Vdrdr*u+(Arr[m]*CosVec[m]+Brr[m]*SinVec[m]);
                    } while(m>0);

                    //The minus signs in the following equations adjust for the
                    //spherical coordinate system difference.
                    d2Vdrdr=crScal*d2Vdrdr/(r*r);
                    d2VdLambdadLambda=crScal*d2VdLambdadLambda;
                    d2VdThetadTheta=crScal*d2VdThetadTheta;
                    d2VdrdLambda=crScal*d2VdrdLambda/r;
                    d2VdrdTheta=-crScal*d2VdrdTheta/r;
                    d2VdLambdadTheta=crScal*d2VdLambdadTheta;

                    if(spherDerivs) {
                        const size_t idx=9*curPoint;
                        
                        HessianVRet[0+idx]=d2Vdrdr;
                        HessianVRet[1+idx]=d2VdrdLambda;
                        HessianVRet[2+idx]=d2VdrdTheta;
                        
                        HessianVRet[3+idx]=d2VdrdLambda;
                        HessianVRet[4+idx]=d2VdLambdadLambda;
                        HessianVRet[5+idx]=d2VdLambdadTheta;
                        
                        HessianVRet[6+idx]=d2VdrdTheta;
                        HessianVRet[7+idx]=d2VdLambdadTheta;
                        HessianVRet[8+idx]=d2VdThetadTheta;

                    } else {
                        double H[27];
                        const double drdx=J[0];
                        const double drdy=J[3];
                        const double drdz=J[6];
                        const double dLambdadx=J[1];
                        const double dLambdady=J[4];
                        const double dLambdadz=J[7];
                        const double dPhidx=J[2];
                        const double dPhidy=J[5];
                        const double dPhidz=J[8];
                        const size_t idx=9*curPoint;

                        calcSpherConvHessianCPP(H,pointCur,0,true);

                        //Commented indices on H are what would be used in
                        //Matlab.
                        const double drdxdx=H[0];//H(1,1,1)
                        const double drdydy=H[4];//H(2,2,1)
                        const double drdzdz=H[8];//H(3,3,1)
                        const double drdxdy=H[3];//H(1,2,1)
                        const double drdxdz=H[6];//H(1,3,1)
                        const double drdydz=H[7];//H(2,3,1)

                        const double dLambdadxdx=H[9];//H(1,1,2)
                        const double dLambdadydy=H[13];//H(2,2,2)
                        const double dLambdadzdz=H[17];//H(3,3,2)
                        const double dLambdadxdy=H[12];//H(1,2,2)
                        const double dLambdadxdz=H[15];//H(1,3,2)
                        const double dLambdadydz=H[16];//H(2,3,2)

                        const double dPhidxdx=H[18];//H(1,1,3)
                        const double dPhidydy=H[22];//H(2,2,3)
                        const double dPhidzdz=H[26];//H(3,3,3)
                        const double dPhidxdy=H[21];//H(1,2,3)
                        const double dPhidxdz=H[24];//H(1,3,3)
                        const double dPhidydz=H[25];//H(2,3,3)

                        //HessianV(1,1,curPoint)
                        HessianVRet[0+idx]=d2VdLambdadLambda*(dLambdadx*dLambdadx)+2.0*d2VdLambdadTheta*dLambdadx*dPhidx+d2VdThetadTheta*(dPhidx*dPhidx)+2.0*dPhidx*drdx*d2VdrdTheta+2.0*dLambdadx*drdx*d2VdrdLambda+dVdLambda*dLambdadxdx+dVdTheta*dPhidxdx+dVdr*drdxdx+(drdx*drdx)*d2Vdrdr;
                        //HessianV(2,2,curPoint)
                        HessianVRet[4+idx]=d2VdThetadTheta*(dPhidy*dPhidy)+2.0*dLambdady*dPhidy*d2VdLambdadTheta+dVdLambda*dLambdadydy+dVdTheta*dPhidydy+(dLambdady*dLambdady)*d2VdLambdadLambda+drdydy*dVdr+2.0*dPhidy*drdy*d2VdrdTheta+2.0*dLambdady*drdy*d2VdrdLambda+(drdy*drdy)*d2Vdrdr;
                        //HessianV(3,3,curPoint)
                        HessianVRet[8+idx]=dVdTheta*dPhidzdz+(dPhidz*dPhidz)*d2VdThetadTheta+dLambdadzdz*dVdLambda+2.0*dLambdadz*dPhidz*d2VdLambdadTheta+(dLambdadz*dLambdadz)*d2VdLambdadLambda+drdzdz*dVdr+2.0*dPhidz*drdz*d2VdrdTheta+2.0*dLambdadz*drdz*d2VdrdLambda+(drdz*drdz)*d2Vdrdr;
                        //HessianV(1,2,curPoint)
                        HessianVRet[3+idx]=dPhidy*d2VdLambdadTheta*dLambdadx+dLambdady*d2VdLambdadLambda*dLambdadx+d2VdThetadTheta*dPhidy*dPhidx+dLambdady*d2VdLambdadTheta*dPhidx+drdy*dPhidx*d2VdrdTheta+dPhidy*drdx*d2VdrdTheta+dVdLambda*dLambdadxdy+dVdTheta*dPhidxdy+dVdr*drdxdy+drdy*dLambdadx*d2VdrdLambda+dLambdady*drdx*d2VdrdLambda+drdy*drdx*d2Vdrdr;
                        //HessianV(1,3,curPoint)
                        HessianVRet[6+idx]=dPhidz*d2VdLambdadTheta*dLambdadx+dLambdadz*d2VdLambdadLambda*dLambdadx+dPhidz*d2VdThetadTheta*dPhidx+dLambdadz*d2VdLambdadTheta*dPhidx+dVdLambda*dLambdadxdz+dVdTheta*dPhidxdz+dVdr*drdxdz+drdz*dPhidx*d2VdrdTheta+dPhidz*drdx*d2VdrdTheta+drdz*dLambdadx*d2VdrdLambda+dLambdadz*drdx*d2VdrdLambda+drdz*drdx*d2Vdrdr;
                        //HessianV(2,3,curPoint)
                        HessianVRet[7+idx]=dPhidz*d2VdThetadTheta*dPhidy+dVdLambda*dLambdadydz+dVdTheta*dPhidydz+dPhidz*dLambdady*d2VdLambdadTheta+dLambdadz*dPhidy*d2VdLambdadTheta+dLambdadz*dLambdady*d2VdLambdadLambda+drdydz*dVdr+drdz*dPhidy*d2VdrdTheta+dPhidz*drdy*d2VdrdTheta+drdz*dLambdady*d2VdrdLambda+dLambdadz*drdy*d2VdrdLambda+drdz*drdy*d2Vdrdr;

                        //HessianV(2,1,curPoint)
                        HessianVRet[1+idx]=HessianVRet[3+idx];
                        //HessianV(3,1,curPoint)
                        HessianVRet[2+idx]=HessianVRet[6+idx];
                        //HessianV(3,2,curPoint)
                        HessianVRet[5+idx]=HessianVRet[7+idx];
                    }
                }
            }
        } else {  
        //At latitudes that are near the poles, the non-singular algorithm
        //of Pines using the fully normalized Helmholtz equations from
        //Fantino and Casotto is used. The algorithm has been slightly
        //modified so that the c/r term is out front and the fully
        //normalized Helmholtz polynomials can be scaled. Also, lumped
        //coefficients are not used. The Pines algorithm is generally
        //slower than the algorithm of Holmes and Featherstone and it
        //suffers a loss of precision near the equator. Thus, the Pines
        //algorithm is only used near the poles where the other algorithm
        //has issues with a singularity.
            double CartPoint[3];

            spher2CartCPP(CartPoint,pointCur,0);

            //Get the direction cosines used by Pines' algorithm.
            const double s=CartPoint[0]/r;
            const double t=CartPoint[1]/r;
            const double u=CartPoint[2]/r;

            //Compute the fully normalized Helmholtz polynomials.
            if(thetaChanged) {
                normHelmHoltzCPP(FuncVals,u, scalFactor);
                
                if(gradVRet!=NULL||HessianVRet!=NULL) {
                    normHelmHoltzDerivCPP(FuncDerivs,FuncVals);
                }
                
                if(HessianVRet!=NULL) {
                    normHelmHoltzDeriv2CPP(FuncDerivs2,FuncVals);   
                }
                
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
            {
                complex <double> V(0,0);
                for(n=0;n<=M;n++) {
                    complex<double> innerTerm(0,0);
                    for(m=0;m<=n;m++) {
                        innerTerm+=(C[n][m]*rm[m]+S[n][m]*im[m])*FuncVals[n][m];
                    }
                    
                    
                    V+=nCoeff[n]*innerTerm;
                }

                V*=crScal;
                VRet[curPoint]=V;
            }

            //If only the gradient is desired as the next output       
            if(gradVRet!=NULL&&HessianVRet==NULL) {
                complex<double> a1(0,0);
                complex<double> a2(0,0);
                complex<double> a3(0,0);
                complex<double> a4(0,0);

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    complex <double> A1(0,0);
                    complex <double> A2(0,0);
                    complex <double> A3(0,0);
                    complex <double> B1(0,0);
                    complex <double> B2(0,0);
                    complex <double> B3(0,0);
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        //The expressions for Lmn, is from Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;
                        const complex<double> rhoC=nCoeff[n]*C[n][m];
                        const complex<double> rhoS=nCoeff[n]*S[n][m];

                        A1+=rhoC*HVal;
                        A2+=rhoC*dHVal;
                        A3+=rhoC*Lmn;

                        B1+=rhoS*HVal;
                        B2+=rhoS*dHVal;
                        B3+=rhoS*Lmn;
                        
                        nf++;
                    }
                    if(M>=1) {
                        a1+=mf*(A1*rm[m-1]+B1*im[m-1]);
                        a2+=mf*(B1*rm[m-1]-A1*im[m-1]);
                    }
                    a3+=(A2*rm[m]+B2*im[m]);
                    a4-=(A3*rm[m]+B3*im[m]);

                    mf++;
                }
                a1/=r;
                a2/=r;
                a3/=r;
                a4/=r;
                
                {
                    const complex<double> dVdx=crScal*(a1+s*a4);
                    const complex<double> dVdy=crScal*(a2+t*a4);
                    const complex<double> dVdz=crScal*(a3+u*a4);

                    if(spherDerivs) {
                        double J[9];
                        complex <double> grad[3];
                        size_t idx=3*curPoint;
                        //Convert the derivatives to spherical coordinates.
                        calcSpherInvJacobCPP(J,pointCur,0);

                        //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                        grad[0]=dVdx*J[0]+dVdy*J[1]+dVdz*J[2];
                        grad[1]=dVdx*J[3]+dVdy*J[4]+dVdz*J[5];
                        grad[2]=dVdx*J[6]+dVdy*J[7]+dVdz*J[8];
                        
                        gradVRet[0+idx]=grad[0];
                        gradVRet[1+idx]=grad[1];
                        gradVRet[2+idx]=grad[2];  
                    } else {//If a gradient in Cartesian coordinates is
                            //desired.
                        size_t idx=3*curPoint;

                        gradVRet[0+idx]=dVdx;
                        gradVRet[1+idx]=dVdy;
                        gradVRet[2+idx]=dVdz;
                    }
                }
            } else if(HessianVRet!=NULL) {
            //Compute the gradient and Hessian
                complex <double> a1(0,0);
                complex <double> a2(0,0);
                complex <double> a3(0,0);
                complex <double> a4(0,0);
                complex <double> a11(0,0);
                complex <double> a12(0,0);
                complex <double> a13(0,0);
                complex <double> a14(0,0);
                complex <double> a23(0,0);
                complex <double> a24(0,0);
                complex <double> a33(0,0);
                complex <double> a34(0,0);
                complex <double> a44(0,0);
                const double r2=r*r;
                const double s2=s*s;
                const double u2=u*u;
                const double t2=t*t;
                complex <double> a22;

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    complex <double> A1(0,0);
                    complex <double> A2(0,0);
                    complex <double> A3(0,0);
                    complex <double> A4(0,0);
                    complex <double> A5(0,0);
                    complex <double> A6(0,0);
                    complex <double> B1(0,0);
                    complex <double> B2(0,0);
                    complex <double> B3(0,0);
                    complex <double> B4(0,0);
                    complex <double> B5(0,0);
                    complex <double> B6(0,0);
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        const double d2HVal=FuncDerivs2[n][m];

                        //The expressions for Lmn, dLmn, and Omn are from
                        //Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;
                        const double dLmn=(nf+mf+2)*dHVal+u*d2HVal;
                        const double Omn=(nf+mf+1)*(n+m+2)*HVal+2.0*u*(nf+mf+2)*dHVal+u*u*d2HVal;

                        const complex <double> rhoC=nCoeff[n]*C[n][m];
                        const complex <double> rhoS=nCoeff[n]*S[n][m];

                        A1+=rhoC*HVal;
                        A2+=rhoC*dHVal;
                        A3+=rhoC*Lmn;
                        A4+=rhoC*d2HVal;
                        A5+=rhoC*dLmn;
                        A6+=rhoC*Omn;

                        B1+=rhoS*HVal;
                        B2+=rhoS*dHVal;
                        B3+=rhoS*Lmn;
                        B4+=rhoS*d2HVal;
                        B5+=rhoS*dLmn;
                        B6+=rhoS*Omn;
                        
                        nf++;
                    }

                    if(m>=1) {
                        a1+=mf*(A1*rm[m-1]+B1*im[m-1]);
                        a2+=mf*(B1*rm[m-1]-A1*im[m-1]);
                        
                        a13+=mf*(A2*rm[m-1]+B2*im[m-1]);
                        a14-=mf*(A3*rm[m-1]+B3*im[m-1]);
                        a23+=mf*(B2*rm[m-1]-A2*im[m-1]);
                        a24-=mf*(B3*rm[m-1]-A3*im[m-1]);
                        
                        if(m>=2) {
                            a11+=mf*(mf-1)*(A1*rm[m-2]+B1*im[m-2]);
                            a12+=mf*(mf-1)*(B1*rm[m-2]-A1*im[m-2]);
                        }
                    }
                    a3+=(A2*rm[m]+B2*im[m]);
                    a4-=(A3*rm[m]+B3*im[m]);

                    a33+=(A4*rm[m]+B4*im[m]);
                    a34-=(A5*rm[m]+B5*im[m]);
                    a44+=(A6*rm[m]+B6*im[m]);

                    mf++;
                }
                a1/=r;
                a2/=r;
                a3/=r;
                a4/=r;
                a11/=r2;
                a12/=r2;
                a13/=r2;
                a14/=r2;
                a23/=r2;
                a24/=r2;
                a33/=r2;
                a34/=r2;
                a44/=r2;
                a22=-a11;

                {
                    double dxdr;
                    double dxdAz;
                    double dxdEl;
                    double dydr;
                    double dydAz;
                    double dydEl;
                    double dzdr;
                    double dzdAz;
                    double dzdEl;
                    const complex<double> dVdx=crScal*(a1+s*a4);
                    const complex<double> dVdy=crScal*(a2+t*a4);
                    const complex<double> dVdz=crScal*(a3+u*a4);
                    const complex<double> d2Vdxdx=crScal*(a11+2.0*s*a14+a4/r+s2*a44-s2*a4/r);
                    const complex<double> d2Vdydy=crScal*(a22+2.0*t*a24+a4/r+t2*a44-t2*a4/r);
                    const complex<double> d2Vdzdz=crScal*(a33+2.0*u*a34+a4/r+u2*a44-u2*a4/r);
                    const complex<double> d2Vdxdy=crScal*(a12+s*t*a44+s*a24+t*a14-s*t*a4/r);
                    const complex<double> d2Vdxdz=crScal*(a13+s*u*a44+s*a34+u*a14-s*u*a4/r);
                    const complex<double> d2Vdydz=crScal*(a23+t*u*a44+t*a34+u*a24-t*u*a4/r);

                    if(spherDerivs) {
                        double J[9];
                        calcSpherInvJacobCPP(J,pointCur,0);

                        dxdr=J[0];//J(1,1);
                        dxdAz=J[3];//J(1,2);
                        dxdEl=J[6];//J(1,3);

                        dydr=J[1];//J(2,1);
                        dydAz=J[4];//J(2,2);
                        dydEl=J[7];//J(2,3);

                        dzdr=J[2];//J(3,1);
                        dzdAz=J[5];//J(3,2);
                        dzdEl=J[8];//J(3,3);
                    }
                    
                    //Compute the gradient, if desired.
                    if(gradVRet!=NULL) {
                        if(spherDerivs) {
                            const size_t idx=3*curPoint;
                            
                            //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                            gradVRet[0+idx]=dVdx*dxdr+dVdy*dydr+dVdz*dzdr;
                            gradVRet[1+idx]=dVdx*dxdAz+dVdy*dydAz+dVdz*dzdAz;
                            gradVRet[2+idx]=dVdx*dxdEl+dVdy*dydEl+dVdz*dzdEl;
                        } else {
                            const size_t idx=3*curPoint;

                            gradVRet[0+idx]=dVdx;
                            gradVRet[1+idx]=dVdy;
                            gradVRet[2+idx]=dVdz;
                        }
                    }
                    
                    //Compute the Hessian.
                    if(spherDerivs) {
                        const size_t idx=9*curPoint;
                        double H[27];
   
                        calcSpherInvHessianCPP(H,pointCur,0);

                        const double d2xdrdr=H[0];//H(1,1,1);
                        const double d2xdAzdAz=H[4];//H(2,2,1);
                        const double d2xdEldEl=H[8];//H(3,3,1);
                        const double d2xdrdAz=H[3];//H(1,2,1);
                        const double d2xdrdEl=H[6];//H(1,3,1);
                        const double d2xdAzdEl=H[7];//H(2,3,1);

                        const double d2ydrdr=H[9];//H(1,1,2);
                        const double d2ydAzdAz=H[13];//H(2,2,2);
                        const double d2ydEldEl=H[17];//H(3,3,2);
                        const double d2ydrdAz=H[12];//H(1,2,2);
                        const double d2ydrdEl=H[15];//H(1,3,2);
                        const double d2ydAzdEl=H[16];//H(2,3,2);

                        const double d2zdrdr=H[18];//H(1,1,3);
                        const double d2zdAzdAz=H[22];//H(2,2,3);
                        const double d2zdEldEl=H[26];//H(3,3,3);
                        const double d2zdrdAz=H[21];//H(1,2,3);
                        const double d2zdrdEl=H[24];//H(1,3,3);
                        const double d2zdAzdEl=H[25];//H(2,3,3);
                        
                        //d2Vdrdr, HessianV(1,1,curPoint)
                        HessianVRet[0+idx]=d2Vdydy*dydr*dydr+2.0*d2Vdydz*dydr*dzdr+d2Vdzdz*dzdr*dzdr+2.0*dxdr*dzdr*d2Vdxdz+2.0*dxdr*dydr*d2Vdxdy+dxdr*dxdr*d2Vdxdx+dVdx*d2xdrdr+dVdy*d2ydrdr+dVdz*d2zdrdr;
                        //d2VdAzdAz, HessianV(2,2,curPoint)
                        HessianVRet[4+idx]=d2Vdzdz*dzdAz*dzdAz+2.0*dydAz*dzdAz*d2Vdydz+dydAz*dydAz*d2Vdydy+dVdy*d2ydAzdAz+dVdz*d2zdAzdAz+d2xdAzdAz*dVdx+2.0*dxdAz*dzdAz*d2Vdxdz+2.0*dxdAz*dydAz*d2Vdxdy+dxdAz*dxdAz*d2Vdxdx;
                        //d2VdEldEl, HessianV(3,3,curPoint)
                        HessianVRet[8+idx]=dzdEl*dzdEl*d2Vdzdz+dVdz*d2zdEldEl+d2ydEldEl*dVdy+2.0*dydEl*dzdEl*d2Vdydz+dydEl*dydEl*d2Vdydy+d2xdEldEl*dVdx+2.0*dxdEl*dzdEl*d2Vdxdz+2.0*dxdEl*dydEl*d2Vdxdy+dxdEl*dxdEl*d2Vdxdx;
                        //d2VdrdAz, HessianV(1,2,curPoint)
                        HessianVRet[3+idx]=dzdAz*d2Vdydz*dydr+dydAz*d2Vdydy*dydr+d2Vdzdz*dzdAz*dzdr+dydAz*d2Vdydz*dzdr+dzdAz*dxdr*d2Vdxdz+dxdAz*dzdr*d2Vdxdz+dydAz*dxdr*d2Vdxdy+dxdAz*dydr*d2Vdxdy+dVdx*d2xdrdAz+dVdy*d2ydrdAz+dVdz*d2zdrdAz+dxdAz*dxdr*d2Vdxdx;
                        //d2VdrdEl, HessianV(1,3,curPoint)
                        HessianVRet[6+idx]=dzdEl*d2Vdydz*dydr+dydEl*d2Vdydy*dydr+dzdEl*d2Vdzdz*dzdr+dydEl*d2Vdydz*dzdr+dzdEl*dxdr*d2Vdxdz+dxdEl*dzdr*d2Vdxdz+dVdx*d2xdrdEl+dVdy*d2ydrdEl+dVdz*d2zdrdEl+dydEl*dxdr*d2Vdxdy+dxdEl*dydr*d2Vdxdy+dxdEl*dxdr*d2Vdxdx;
                        //d2VdAzdEl, HessianV(2,3,curPoint)
                        HessianVRet[7+idx]=dzdEl*d2Vdzdz*dzdAz+dzdEl*dydAz*d2Vdydz+dydEl*dzdAz*d2Vdydz+dVdy*d2ydAzdEl+dVdz*d2zdAzdEl+dydEl*dydAz*d2Vdydy+d2xdAzdEl*dVdx+dzdEl*dxdAz*d2Vdxdz+dxdEl*dzdAz*d2Vdxdz+dydEl*dxdAz*d2Vdxdy+dxdEl*dydAz*d2Vdxdy+dxdEl*dxdAz*d2Vdxdx;

 
                        //d2VdAzdr, HessianV(2,1,curPoint)
                        HessianVRet[1+idx]=HessianVRet[3+idx];
                        //d2VdEldr, HessianV(3,1,curPoint)
                        HessianVRet[2+idx]=HessianVRet[6+idx];
                        //d2VdEldAz, HessianV(3,2,curPoint)
                        HessianVRet[5+idx]=HessianVRet[7+idx];
                    } else {
                        const size_t idx=9*curPoint;
                        
                        //HessianV(1,1,curPoint)
                        HessianVRet[0+idx]=d2Vdxdx;
                        //HessianV(2,2,curPoint)
                        HessianVRet[4+idx]=d2Vdydy;
                        //HessianV(3,3,curPoint)
                        HessianVRet[8+idx]=d2Vdzdz;
                        //HessianV(1,2,curPoint)
                        HessianVRet[3+idx]=d2Vdxdy;
                        //HessianV(2,1,curPoint)
                        HessianVRet[1+idx]=HessianVRet[3+idx];
                        //HessianV(1,3,curPoint)
                        HessianVRet[6+idx]=d2Vdxdz;
                        //HessianV(3,1,curPoint)
                        HessianVRet[2+idx]=HessianVRet[6+idx];
                        //HessianV(2,3,curPoint)
                        HessianVRet[7+idx]=d2Vdydz;
                        //HessianV(3,2,curPoint)
                        HessianVRet[5+idx]=HessianVRet[7+idx];
                    }
                }
            }
        }
        
        if(spherDerivs&&systemType==2) {
            size_t idx;
            //Flip signs of the elevation terms reflecting the
            //difference in definition between systems 0 and 2.
            
            if(gradVRet!=NULL) {
                idx=3*curPoint;
                gradVRet[2+idx]=-gradVRet[2+idx];
            }

            if(HessianVRet!=NULL) {
                idx=9*curPoint;
                //HessianV(1,3,curPoint)
                HessianVRet[6+idx]=-HessianVRet[6+idx];
                //HessianV(3,1,curPoint)
                HessianVRet[2+idx]=-HessianVRet[2+idx];
                //HessianV(2,3,curPoint)
                HessianVRet[7+idx]=-HessianVRet[7+idx];
                //HessianV(3,2,curPoint)
                HessianVRet[5+idx]=-HessianVRet[5+idx];
            }
        }
    }
    
    delete[] buffer;
    delete[] bufferComplex;
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
