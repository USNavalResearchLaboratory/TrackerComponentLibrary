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
//For max
#include <algorithm>

using namespace std;

void spherHarmonicSetEvalCPPReal(double *V, double *gradV, double *HessianV,const CountingClusterSetVecCPP<double> &C,const CountingClusterSetVecCPP<double> &S, const double *point, const size_t numPoints, const double a, const double c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm) {
    //If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. If a NULL pointer is passed for HessianV,
    //then it is assumed that the Hessian is not desired.
    const size_t numSets=C.numSets;
    const size_t M=C.numClust-1;
    const size_t M1=C.numClust;
    const size_t M1NumSets=M1*numSets;
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
    double *nCoeff;//Length M1
    //These are to store sin(m*lambda) and cos(m*lambda) for m=0->M for the
    //Legendre algorithm.
    double *SinVec,*CosVec;//These are each length M+1.
    //rm and im are never used at the same time as SinVec andCosVec,
    //because rm and im are used by Pines' algorithm. They will be the same
    //size as SinVec and CosVec, and thus will point to the same memory.
    double *rm,*im;
    //The following are needed if the Legendre algorithm is used.
    double *A=NULL;//Length M1NumSets
    double *B=NULL;//Length M1NumSets
    //The following are only used if the gradient or Hessian is desired and
    //the Legendre algorithm is used.
    double *Ar=NULL;//Length M1NumSets
    double *Br=NULL;//Length M1NumSets
    double *ATheta=NULL;//Length M1NumSets
    double *BTheta=NULL;//Length M1NumSets
    //The following are only used if the Hessian is desired and the
    //Legendre algorithm is used.
    double *AThetaTheta=NULL;//Length M1NumSets
    double *BThetaTheta=NULL;//Length M1NumSets
    double *Arr=NULL;//Length M1NumSets
    double *Brr=NULL;//Length M1NumSets
    double *AThetar=NULL;//Length M1NumSets
    double *BThetar=NULL;//Length M1NumSets
    //The following are only needed if Pines' algorithm is chosen and the
    //gradient and/or Hessian is desired
    double *a1;//Length numSets
    double *a2;//Length numSets
    double *a3;//Length numSets
    double *a4;//Length numSets
    double *A1=NULL;//Length numSets
    double *A2=NULL;//Length numSets
    double *A3=NULL;//Length numSets
    double *B1=NULL;//Length numSets
    double *B2=NULL;//Length numSets
    double *B3=NULL;//Length numSets
    //The following are only needed if Pines' algorithm is chosen and the
    //Hessian is desired. 
    double *a11=NULL;//Length numSets
    double *a12=NULL;//Length numSets
    double *a13=NULL;//Length numSets
    double *a14=NULL;//Length numSets
    double *a23=NULL;//Length numSets
    double *a24=NULL;//Length numSets
    double *a33=NULL;//Length numSets
    double *a34=NULL;//Length numSets
    double *a44=NULL;//Length numSets
    double *A4=NULL;//Length numSets
    double *A5=NULL;//Length numSets
    double *A6=NULL;//Length numSets
    double *B4=NULL;//Length numSets
    double *B5=NULL;//Length numSets
    double *B6=NULL;//Length numSets

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
        }else if(gradV!=NULL&&HessianV==NULL) {
            buffer=new double[3*M1+3*C.totalNumEl+10*numSets];
        } else {
            buffer=new double[3*M1+3*C.totalNumEl+25*numSets];
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
            
            tempPtr+=C.totalNumEl;
            a1=tempPtr;//Length numSets

            tempPtr+=numSets;
            a2=tempPtr;//Length numSets

            tempPtr+=numSets;
            a3=tempPtr;//Length numSets

            tempPtr+=numSets;
            a4=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            A1=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            A2=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            A3=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            B1=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            B2=tempPtr;//Length numSets
            
            tempPtr+=numSets;
            B3=tempPtr;//Length numSets

            if(HessianV!=NULL) {
                tempPtr+=numSets;
                a11=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a12=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a13=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a14=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a23=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a24=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a33=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a34=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                a44=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                A4=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                A5=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                A6=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                B4=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                B5=tempPtr;//Length numSets
                
                tempPtr+=numSets;
                B6=tempPtr;//Length numSets
            }
        }
    } else if(algorithm==1) {//If only Legendre's algorithm is used.
        double *tempPtr;
        
        if(gradV==NULL&&HessianV==NULL) {
            buffer=new double[3*M1+2*M1NumSets+C.totalNumEl];
        }else if(gradV!=NULL&&HessianV==NULL) {
            buffer=new double[3*M1+6*M1NumSets+2*C.totalNumEl];
        } else {
            buffer=new double[3*M1+12*M1NumSets+3*C.totalNumEl];
        }

        tempPtr=buffer;
        nCoeff=tempPtr;//Length M1

        tempPtr+=M1;
        SinVec=tempPtr;//Length M1

        tempPtr+=M1;
        CosVec=tempPtr;//Length M1

        tempPtr+=M1;
        A=tempPtr;//Length M1NumSets

        tempPtr+=M1NumSets;
        B=tempPtr;//Length M1NumSets

        tempPtr+=M1NumSets;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        if(gradV!=NULL||HessianV!=NULL) {
            tempPtr+=C.totalNumEl;
            Ar=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            Br=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            ATheta=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            BTheta=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl
            
             if(HessianV!=NULL) {
                tempPtr+=C.totalNumEl;
                FuncDerivs2.clusterEls=tempPtr;//Length C.totalNumEl

                tempPtr+=C.totalNumEl;
                AThetaTheta=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                BThetaTheta=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                Arr=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                Brr=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                AThetar=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                BThetar=tempPtr;//Length M1NumSets
            }
        }
    } else {
        double *tempPtr, *tempPtrPines;
        //If the Legendre and Pines algorithm can be used.
        
        if(HessianV!=NULL) {
            buffer=new double[3*M1+2*M1NumSets+3*C.totalNumEl+max(4*M1NumSets,10*numSets)+max(6*M1NumSets,15*numSets)];
        } else if(gradV!=NULL) {
            //The 4*M1NumSets is memory for the Legendre algorithm.
            //The 10NumSets is memory for Pines' algorithm.
            buffer=new double[3*M1+2*M1NumSets+3*C.totalNumEl+max(4*M1NumSets,10*numSets)];
        } else {
            buffer=new double[3*M1+2*M1NumSets+C.totalNumEl];
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
        A=tempPtr;//Length M1NumSets

        tempPtr+=M1NumSets;
        B=tempPtr;//Length M1NumSets

        tempPtr+=M1NumSets;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        if(gradV!=NULL||HessianV!=NULL) {
            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl

            tempPtr+=C.totalNumEl;
            FuncDerivs2.clusterEls=tempPtr;// Length C.totalNumEl 

            tempPtr+=C.totalNumEl;
            tempPtrPines=tempPtr;
            Ar=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            Br=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            ATheta=tempPtr;//Length M1NumSets

            tempPtr+=M1NumSets;
            BTheta=tempPtr;//Length M1NumSets

            a1=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            a2=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            a3=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            a4=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            A1=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            A2=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            A3=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            B1=tempPtrPines;//Length numSets

            tempPtrPines+=numSets;
            B2=tempPtrPines;//Length numSets

            tempPtr+=numSets;
            B3=tempPtrPines;//Length numSets
            
            if(HessianV!=NULL) {
                tempPtr+=M1NumSets;
                AThetaTheta=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                BThetaTheta=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                Arr=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                Brr=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                AThetar=tempPtr;//Length M1NumSets

                tempPtr+=M1NumSets;
                BThetar=tempPtr;//Length M1NumSets

                //Pines Part

                tempPtrPines+=numSets;
                a11=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a12=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a13=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a14=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a23=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a24=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a33=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a34=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                a44=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                A4=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                A5=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                A6=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                B4=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                B5=tempPtrPines;//Length numSets

                tempPtrPines+=numSets;
                B6=tempPtrPines;//Length numSets
            }
        }
    }

    nCoeff[0]=1;
    
    rPrev=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double pointCur[3];
        double r, lambda, thetaCur;
        size_t n,m,curSet;
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
                memset(A,0,sizeof(double)*M1NumSets);
                memset(B,0,sizeof(double)*M1NumSets);

                //Compute the X coefficients for the sum
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t curSetOffset=curSet*M1;
                    
                    for(m=0;m<=M;m++) {
                        for(n=m;n<=M;n++) {
                            const size_t k=m+curSetOffset;
                            A[k]+=nCoeff[n]*C(n,m,curSet)*FuncVals[n][m];
                            B[k]+=nCoeff[n]*S(n,m,curSet)*FuncVals[n][m];
                        }
                    }
                }

                //If additional terms should be computed so a gradient or
                //Hessian can be computed.
                if(gradV!=NULL||HessianV!=NULL) {
                    memset(Ar,0,sizeof(double)*M1NumSets);
                    memset(Br,0,sizeof(double)*M1NumSets);
                    memset(ATheta,0,sizeof(double)*M1NumSets);
                    memset(BTheta,0,sizeof(double)*M1NumSets);
 
                    for(curSet=0;curSet<numSets;curSet++) {
                        const size_t curSetOffset=curSet*M1;
                        
                        mf=0;
                        for(m=0;m<=M;m++) {
                            const size_t k=m+curSetOffset;
                            nf=mf;
                            for(n=m;n<=M;n++) {
                                const double CScal=nCoeff[n]*C(n,m,curSet);
                                const double SScal=nCoeff[n]*S(n,m,curSet);
                                double curVal;

                                curVal=FuncVals[n][m];
                                Ar[k]+=(nf+1)*CScal*curVal;
                                Br[k]+=(nf+1)*SScal*curVal;

                                curVal=FuncDerivs[n][m];
                                ATheta[k]+=CScal*curVal;
                                BTheta[k]+=SScal*curVal;

                                nf++;
                            }

                            mf++;
                        }
                    }
                    
                    if(HessianV!=NULL) {
                        memset(Arr,0,sizeof(double)*M1NumSets);
                        memset(Brr,0,sizeof(double)*M1NumSets);
                        memset(AThetar,0,sizeof(double)*M1NumSets);
                        memset(BThetar,0,sizeof(double)*M1NumSets);
                        memset(AThetaTheta,0,sizeof(double)*M1NumSets);
                        memset(BThetaTheta,0,sizeof(double)*M1NumSets);
                        
                        for(curSet=0;curSet<numSets;curSet++) {
                            const size_t curSetOffset=curSet*M1;
                            
                            mf=0;
                            for(m=0;m<=M;m++) {
                                const size_t k=m+curSetOffset;
                                nf=mf;
                                for(n=m;n<=M;n++) {
                                    const double CScal=nCoeff[n]*C(n,m,curSet);
                                    const double SScal=nCoeff[n]*S(n,m,curSet);
                                    double curVal;

                                    curVal=FuncVals[n][m];
                                    //From Table 5, with the correction
                                    //from the erratum.
                                    Arr[k]+=(nf+1)*(nf+2)*CScal*curVal;
                                    Brr[k]+=(nf+1)*(nf+2)*SScal*curVal;;

                                    curVal=FuncDerivs[n][m];
                                    //From Table 5
                                    AThetar[k]+=(nf+1)*CScal*curVal;
                                    BThetar[k]+=(nf+1)*SScal*curVal;

                                    curVal=FuncDerivs2[n][m];
                                    //From Table 5
                                    AThetaTheta[k]+=CScal*curVal;
                                    BThetaTheta[k]+=SScal*curVal;

                                    nf++;
                                }

                                mf++;
                            }
                        }
                    }
                }
            }
            
            //Use Horner's method to compute V.
            for(curSet=0;curSet<numSets;curSet++) {
                const size_t curSetOffset=curSet*M1;
                const size_t idx=curSet+curPoint*numSets;

                V[idx]=0;
                m=M+1;
                do {
                    m--;
                    const size_t k=m+curSetOffset;

                    V[idx]=V[idx]*u+A[k]*CosVec[m]+B[k]*SinVec[m];
                } while(m>0);
                
                //Multiply by the constant in front of the sum and get rid
                //of the scale factor.
                V[idx]=crScal*V[idx];
            }

            //Compute the gradient, if it is desired.
            if(gradV!=NULL||HessianV!=NULL) {
                double J[9];
                double H[27];
                double dVdr;
                double dVdLambda;
                double dVdTheta;
                double d2VdLambdadLambda;
                double d2VdLambdadTheta;
                double d2VdrdLambda;
                double d2VdThetadTheta;
                double d2VdrdTheta;
                double d2Vdrdr;                                                
                double drdx;
                double drdy;
                double drdz;
                double dLambdadx;
                double dLambdady;
                double dLambdadz;
                double dPhidx;
                double dPhidy;
                double dPhidz;
                double drdxdx;
                double drdydy;
                double drdzdz;
                double drdxdy;
                double drdxdz;
                double drdydz;
                double dLambdadxdx;
                double dLambdadydy;
                double dLambdadzdz;
                double dLambdadxdy;
                double dLambdadxdz;
                double dLambdadydz;
                double dPhidxdx;
                double dPhidydy;
                double dPhidzdz;
                double dPhidxdy;
                double dPhidxdz;
                double dPhidydz;
                
                if(spherDerivs==false) {
                    calcSpherConvJacobCPP(J, pointCur,0);
                    
                    if(HessianV!=NULL) {
                        drdx=J[0];
                        drdy=J[3];
                        drdz=J[6];
                        dLambdadx=J[1];
                        dLambdady=J[4];
                        dLambdadz=J[7];
                        dPhidx=J[2];
                        dPhidy=J[5];
                        dPhidz=J[8];

                        calcSpherConvHessianCPP(H,pointCur,0,true);

                        //Commented indices on H are what would be used in
                        //Matlab.
                        drdxdx=H[0];//H(1,1,1)
                        drdydy=H[4];//H(2,2,1)
                        drdzdz=H[8];//H(3,3,1)
                        drdxdy=H[3];//H(1,2,1)
                        drdxdz=H[6];//H(1,3,1)
                        drdydz=H[7];//H(2,3,1)

                        dLambdadxdx=H[9];//H(1,1,2)
                        dLambdadydy=H[13];//H(2,2,2)
                        dLambdadzdz=H[17];//H(3,3,2)
                        dLambdadxdy=H[12];//H(1,2,2)
                        dLambdadxdz=H[15];//H(1,3,2)
                        dLambdadydz=H[16];//H(2,3,2)

                        dPhidxdx=H[18];//H(1,1,3)
                        dPhidydy=H[22];//H(2,2,3)
                        dPhidzdz=H[26];//H(3,3,3)
                        dPhidxdy=H[21];//H(1,2,3)
                        dPhidxdz=H[24];//H(1,3,3)
                        dPhidydz=H[25];//H(2,3,3)
                    }
                }

                for(curSet=0;curSet<numSets;curSet++) {
                    dVdr=0;
                    dVdLambda=0;
                    dVdTheta=0;
                    //Use Horner's method to compute all dV values.
                    m=M+1;
                    mf=static_cast<double>(m);
                    do {
                        m--;
                        mf--;
                        const size_t k=m+curSet*M1;

                        dVdr=dVdr*u-(Ar[k]*CosVec[m]+Br[k]*SinVec[m]);
                        dVdLambda=dVdLambda*u-mf*(A[k]*SinVec[m]-B[k]*CosVec[m]);
                        dVdTheta=dVdTheta*u+(ATheta[k]*CosVec[m]+BTheta[k]*SinVec[m]);
                    } while(m>0);            

                    dVdr=crScal*dVdr/r;
                    dVdLambda=crScal*dVdLambda;
                    //The minus sign adjusts for the coordinate system change.
                    dVdTheta=-crScal*dVdTheta;

                    if(spherDerivs&&gradV!=NULL) {
                        const size_t idx=3*curSet+3*numSets*curPoint;
                        gradV[0+idx]=dVdr;
                        gradV[1+idx]=dVdLambda;
                        gradV[2+idx]=dVdTheta;
                    } else {
                        const size_t idx=3*curSet+3*numSets*curPoint;
                        calcSpherConvJacobCPP(J,pointCur,0);

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
                        d2VdLambdadLambda=0;
                        d2VdLambdadTheta=0;
                        d2VdrdLambda=0;
                        d2VdThetadTheta=0;
                        d2VdrdTheta=0;
                        d2Vdrdr=0;

                        //The following second-order derivative formulae are
                        //from Table 2 (expressed using Horner's method).
                        m=M1;
                        mf=static_cast<double>(m);
                        do {
                            m--;
                            mf--;
                            const size_t k=m+curSet*M1;
                            
                            d2VdLambdadLambda=d2VdLambdadLambda*u-mf*mf*(A[k]*CosVec[m]+B[k]*SinVec[m]);
                            d2VdLambdadTheta=d2VdLambdadTheta*u+mf*(ATheta[k]*SinVec[m]-BTheta[k]*CosVec[m]);
                            d2VdrdLambda=d2VdrdLambda*u+mf*(Ar[k]*SinVec[m]-Br[k]*CosVec[m]);
                            d2VdThetadTheta=d2VdThetadTheta*u+(AThetaTheta[k]*CosVec[m]+BThetaTheta[k]*SinVec[m]);
                            d2VdrdTheta=d2VdrdTheta*u-(AThetar[k]*CosVec[m]+BThetar[k]*SinVec[m]);
                            d2Vdrdr=d2Vdrdr*u+(Arr[k]*CosVec[m]+Brr[k]*SinVec[m]);
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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
            for(curSet=0;curSet<numSets;curSet++) {
                const size_t idx=curSet+curPoint*numSets;
                
                V[idx]=0;
                for(n=0;n<=M;n++) {
                    double innerTerm=0;
                    for(m=0;m<=n;m++) {
                        innerTerm+=(C(n,m,curSet)*rm[m]+S(n,m,curSet)*im[m])*FuncVals[n][m];
                    }
                    V[idx]+=nCoeff[n]*innerTerm;
                }

                V[idx]=crScal*V[idx];
            }
            
            //If only the gradient is desired as the next output       
            if(gradV!=NULL&&HessianV==NULL) {
                memset(a1,0,sizeof(double)*numSets);
                memset(a2,0,sizeof(double)*numSets);
                memset(a3,0,sizeof(double)*numSets);
                memset(a4,0,sizeof(double)*numSets);

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    memset(A1,0,sizeof(double)*numSets);
                    memset(A2,0,sizeof(double)*numSets);
                    memset(A3,0,sizeof(double)*numSets);
                    memset(B1,0,sizeof(double)*numSets);
                    memset(B2,0,sizeof(double)*numSets);
                    memset(B3,0,sizeof(double)*numSets);
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        //The expressions for Lmn, is from Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;  
                        for(curSet=0;curSet<numSets;curSet++) {
                            const double rhoC=nCoeff[n]*C(n,m,curSet);
                            const double rhoS=nCoeff[n]*S(n,m,curSet);

                            A1[curSet]+=rhoC*HVal;
                            A2[curSet]+=rhoC*dHVal;
                            A3[curSet]+=rhoC*Lmn;

                            B1[curSet]+=rhoS*HVal;
                            B2[curSet]+=rhoS*dHVal;
                            B3[curSet]+=rhoS*Lmn;
                        }
                        
                        nf++;
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        if(M>=1) {
                            a1[curSet]+=mf*(A1[curSet]*rm[m-1]+B1[curSet]*im[m-1]);
                            a2[curSet]+=mf*(B1[curSet]*rm[m-1]-A1[curSet]*im[m-1]);
                        }
                        a3[curSet]+=(A2[curSet]*rm[m]+B2[curSet]*im[m]);
                        a4[curSet]-=(A3[curSet]*rm[m]+B3[curSet]*im[m]);
                    }

                    mf++;
                }
                
                {
                    double J[9];
                    if(spherDerivs) {
                        calcSpherInvJacobCPP(J,pointCur,0);
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        a1[curSet]=a1[curSet]/r;
                        a2[curSet]=a2[curSet]/r;
                        a3[curSet]=a3[curSet]/r;
                        a4[curSet]=a4[curSet]/r;

                        const double dVdx=crScal*(a1[curSet]+s*a4[curSet]);
                        const double dVdy=crScal*(a2[curSet]+t*a4[curSet]);
                        const double dVdz=crScal*(a3[curSet]+u*a4[curSet]);
                    
                        if(spherDerivs) {
                            //Convert the derivatives to spherical coordinates.
                            const size_t idx=3*curSet+3*numSets*curPoint;
                            //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];

                            gradV[0+idx]=dVdx*J[0]+dVdy*J[1]+dVdz*J[2];
                            gradV[1+idx]=dVdx*J[3]+dVdy*J[4]+dVdz*J[5];
                            gradV[2+idx]=dVdx*J[6]+dVdy*J[7]+dVdz*J[8];

                        } else {//If a gradient in Cartesian coordinates is desired.
                            const size_t idx=3*curSet+3*numSets*curPoint;

                            gradV[0+idx]=dVdx;
                            gradV[1+idx]=dVdy;
                            gradV[2+idx]=dVdz;
                        }
                    }
                }
            } else if(HessianV!=NULL) {
            //Compute the gradient and Hessian
                memset(a1,0,sizeof(double)*numSets);
                memset(a2,0,sizeof(double)*numSets);
                memset(a3,0,sizeof(double)*numSets);
                memset(a4,0,sizeof(double)*numSets);
                memset(a11,0,sizeof(double)*numSets);
                memset(a12,0,sizeof(double)*numSets);
                memset(a13,0,sizeof(double)*numSets);
                memset(a14,0,sizeof(double)*numSets);
                memset(a23,0,sizeof(double)*numSets);
                memset(a24,0,sizeof(double)*numSets);
                memset(a33,0,sizeof(double)*numSets);
                memset(a34,0,sizeof(double)*numSets);
                memset(a44,0,sizeof(double)*numSets);
                
                const double r2=r*r;
                const double s2=s*s;
                const double u2=u*u;
                const double t2=t*t;

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    memset(A1,0,sizeof(double)*numSets);
                    memset(A2,0,sizeof(double)*numSets);
                    memset(A3,0,sizeof(double)*numSets);
                    memset(A4,0,sizeof(double)*numSets);
                    memset(A5,0,sizeof(double)*numSets);
                    memset(A6,0,sizeof(double)*numSets);
                    memset(B1,0,sizeof(double)*numSets);
                    memset(B2,0,sizeof(double)*numSets);
                    memset(B3,0,sizeof(double)*numSets);
                    memset(B4,0,sizeof(double)*numSets);
                    memset(B5,0,sizeof(double)*numSets);
                    memset(B6,0,sizeof(double)*numSets);
                    
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
                        
                        for(curSet=0;curSet<numSets;curSet++) {
                            const double rhoC=nCoeff[n]*C(n,m,curSet);
                            const double rhoS=nCoeff[n]*S(n,m,curSet);

                            A1[curSet]+=rhoC*HVal;
                            A2[curSet]+=rhoC*dHVal;
                            A3[curSet]+=rhoC*Lmn;
                            A4[curSet]+=rhoC*d2HVal;
                            A5[curSet]+=rhoC*dLmn;
                            A6[curSet]+=rhoC*Omn;

                            B1[curSet]+=rhoS*HVal;
                            B2[curSet]+=rhoS*dHVal;
                            B3[curSet]+=rhoS*Lmn;
                            B4[curSet]+=rhoS*d2HVal;
                            B5[curSet]+=rhoS*dLmn;
                            B6[curSet]+=rhoS*Omn;
                        }
                        
                        nf++;
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        if(m>=1) {
                            a1[curSet]+=mf*(A1[curSet]*rm[m-1]+B1[curSet]*im[m-1]);
                            a2[curSet]+=mf*(B1[curSet]*rm[m-1]-A1[curSet]*im[m-1]);

                            a13[curSet]+=mf*(A2[curSet]*rm[m-1]+B2[curSet]*im[m-1]);
                            a14[curSet]-=mf*(A3[curSet]*rm[m-1]+B3[curSet]*im[m-1]);
                            a23[curSet]+=mf*(B2[curSet]*rm[m-1]-A2[curSet]*im[m-1]);
                            a24[curSet]-=mf*(B3[curSet]*rm[m-1]-A3[curSet]*im[m-1]);

                            if(m>=2) {
                                a11[curSet]+=mf*(mf-1)*(A1[curSet]*rm[m-2]+B1[curSet]*im[m-2]);
                                a12[curSet]+=mf*(mf-1)*(B1[curSet]*rm[m-2]-A1[curSet]*im[m-2]);
                            }
                        }
                        a3[curSet]+=(A2[curSet]*rm[m]+B2[curSet]*im[m]);
                        a4[curSet]-=(A3[curSet]*rm[m]+B3[curSet]*im[m]);

                        a33[curSet]+=(A4[curSet]*rm[m]+B4[curSet]*im[m]);
                        a34[curSet]-=(A5[curSet]*rm[m]+B5[curSet]*im[m]);
                        a44[curSet]+=(A6[curSet]*rm[m]+B6[curSet]*im[m]);
                    }

                    mf++;
                }

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
                    double d2xdrdr;
                    double d2xdAzdAz;
                    double d2xdEldEl;
                    double d2xdrdAz;
                    double d2xdrdEl;
                    double d2xdAzdEl;
                    double d2ydrdr;
                    double d2ydAzdAz;
                    double d2ydEldEl;
                    double d2ydrdAz;
                    double d2ydrdEl;
                    double d2ydAzdEl;
                    double d2zdrdr;
                    double d2zdAzdAz;
                    double d2zdEldEl;
                    double d2zdrdAz;
                    double d2zdrdEl;
                    double d2zdAzdEl;
                    
                    if(spherDerivs) {
                        double J[9];
                        double H[27];

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
                        
                        calcSpherInvHessianCPP(H,pointCur,0);
                        d2xdrdr=H[0];//H(1,1,1);
                        d2xdAzdAz=H[4];//H(2,2,1);
                        d2xdEldEl=H[8];//H(3,3,1);
                        d2xdrdAz=H[3];//H(1,2,1);
                        d2xdrdEl=H[6];//H(1,3,1);
                        d2xdAzdEl=H[7];//H(2,3,1);

                        d2ydrdr=H[9];//H(1,1,2);
                        d2ydAzdAz=H[13];//H(2,2,2);
                        d2ydEldEl=H[17];//H(3,3,2);
                        d2ydrdAz=H[12];//H(1,2,2);
                        d2ydrdEl=H[15];//H(1,3,2);
                        d2ydAzdEl=H[16];//H(2,3,2);

                        d2zdrdr=H[18];//H(1,1,3);
                        d2zdAzdAz=H[22];//H(2,2,3);
                        d2zdEldEl=H[26];//H(3,3,3);
                        d2zdrdAz=H[21];//H(1,2,3);
                        d2zdrdEl=H[24];//H(1,3,3);
                        d2zdAzdEl=H[25];//H(2,3,3);
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        a1[curSet]/=r;
                        a2[curSet]/=r;
                        a3[curSet]/=r;
                        a4[curSet]/=r;
                        a11[curSet]/=r2;
                        a12[curSet]/=r2;
                        a13[curSet]/=r2;
                        a14[curSet]/=r2;
                        a23[curSet]/=r2;
                        a24[curSet]/=r2;
                        a33[curSet]/=r2;
                        a34[curSet]/=r2;
                        a44[curSet]/=r2;

                        const double dVdx=crScal*(a1[curSet]+s*a4[curSet]);
                        const double dVdy=crScal*(a2[curSet]+t*a4[curSet]);
                        const double dVdz=crScal*(a3[curSet]+u*a4[curSet]);
                        const double d2Vdxdx=crScal*(a11[curSet]+2.0*s*a14[curSet]+a4[curSet]/r+s2*a44[curSet]-s2*a4[curSet]/r);
                        const double d2Vdydy=crScal*(-a11[curSet]+2.0*t*a24[curSet]+a4[curSet]/r+t2*a44[curSet]-t2*a4[curSet]/r);
                        const double d2Vdzdz=crScal*(a33[curSet]+2.0*u*a34[curSet]+a4[curSet]/r+u2*a44[curSet]-u2*a4[curSet]/r);
                        const double d2Vdxdy=crScal*(a12[curSet]+s*t*a44[curSet]+s*a24[curSet]+t*a14[curSet]-s*t*a4[curSet]/r);
                        const double d2Vdxdz=crScal*(a13[curSet]+s*u*a44[curSet]+s*a34[curSet]+u*a14[curSet]-s*u*a4[curSet]/r);
                        const double d2Vdydz=crScal*(a23[curSet]+t*u*a44[curSet]+t*a34[curSet]+u*a24[curSet]-t*u*a4[curSet]/r);
                    
                        //Compute the gradient, if desired.
                        if(gradV!=NULL) {
                            if(spherDerivs) {
                                const size_t idx=3*curSet+3*numSets*curPoint;

                                //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                                gradV[0+idx]=dVdx*dxdr+dVdy*dydr+dVdz*dzdr;
                                gradV[1+idx]=dVdx*dxdAz+dVdy*dydAz+dVdz*dzdAz;
                                gradV[2+idx]=dVdx*dxdEl+dVdy*dydEl+dVdz*dzdEl;
                            } else {
                                const size_t idx=3*curSet+3*numSets*curPoint;

                                gradV[0+idx]=dVdx;
                                gradV[1+idx]=dVdy;
                                gradV[2+idx]=dVdz;
                            }
                        }
                        
                        //Compute the Hessian.
                        if(spherDerivs) {
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
        }
        
        if(spherDerivs&&systemType==2) {
            //Flip signs of the elevation terms reflecting the
            //difference in definition between systems 0 and 2.
            
            if(gradV!=NULL) {
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t idx=3*curSet+3*numSets*curPoint;
                    gradV[2+idx]=-gradV[2+idx];
                }
            }

            if(HessianV!=NULL) {
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t idx=9*curSet+9*numSets*curPoint;
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
    }
    
    delete[] buffer;
}

void spherHarmonicSetEvalCPPComplex(complex<double> *VRet, complex<double> *gradVRet, complex<double> *HessianVRet,const CountingClusterSetVecCPP<complex<double>> &C,const CountingClusterSetVecCPP<complex<double>> &S, const double *point, const size_t numPoints, const complex <double> a, const complex <double> c, const size_t systemType, const bool spherDerivs,const double scalFactor,const size_t algorithm) {
//If a NULL pointer is passed for gradV, then it is assumed that the
    //gradient is not desired. If a NULL pointer is passed for HessianV,
    //then it is assumed that the Hessian is not desired.
    const size_t numSets=C.numSets;
    const size_t M=C.numClust-1;
    const size_t M1=C.numClust;
    const size_t M1NumSets=M1*numSets;
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
    complex<double> *A=NULL;//Length M1NumSets
    complex<double> *B=NULL;//Length M1NumSets
    //The following are only used if the gradient or Hessian is desired and
    //the Legendre algorithm is used.
    complex<double> *Ar=NULL;//Length M1NumSets
    complex<double> *Br=NULL;//Length M1NumSets
    complex<double> *ATheta=NULL;//Length M1NumSets
    complex<double> *BTheta=NULL;//Length M1NumSets
    //The following are only used if the Hessian is desired and the
    //Legendre algorithm is used.
    complex<double> *AThetaTheta=NULL;//Length M1NumSets
    complex<double> *BThetaTheta=NULL;//Length M1NumSets
    complex<double> *Arr=NULL;//Length M1NumSets
    complex<double> *Brr=NULL;//Length M1NumSets
    complex<double> *AThetar=NULL;//Length M1NumSets
    complex<double> *BThetar;//Length M1NumSets
    //The following are only needed if Pines' algorithm is chosen and the
    //gradient and/or Hessian is desired
    complex<double> *a1;//Length numSets
    complex<double> *a2;//Length numSets
    complex<double> *a3;//Length numSets
    complex<double> *a4;//Length numSets
    complex<double> *A1=NULL;//Length numSets
    complex<double> *A2=NULL;//Length numSets
    complex<double> *A3=NULL;//Length numSets
    complex<double> *B1=NULL;//Length numSets
    complex<double> *B2=NULL;//Length numSets
    complex<double> *B3=NULL;//Length numSets
    //The following are only needed if Pines' algorithm is chosen and the
    //Hessian is desired. 
    complex<double> *a11=NULL;//Length numSets
    complex<double> *a12=NULL;//Length numSets
    complex<double> *a13=NULL;//Length numSets
    complex<double> *a14=NULL;//Length numSets
    complex<double> *a23=NULL;//Length numSets
    complex<double> *a24=NULL;//Length numSets
    complex<double> *a33=NULL;//Length numSets
    complex<double> *a34=NULL;//Length numSets
    complex<double> *a44=NULL;//Length numSets
    complex<double> *A4=NULL;//Length numSets
    complex<double> *A5=NULL;//Length numSets
    complex<double> *A6=NULL;//Length numSets
    complex<double> *B4=NULL;//Length numSets
    complex<double> *B5=NULL;//Length numSets
    complex<double> *B6=NULL;//Length numSets
    
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
        complex<double> *tempPtrComplex;

        if(gradVRet==NULL&&HessianVRet==NULL) {
            buffer=new double[3*M1+C.totalNumEl];
            bufferComplex=new complex<double>[M1];
        }else if(gradVRet!=NULL&&HessianVRet==NULL) {
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex<double>[M1+10*numSets];
        } else {
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex<double>[M1+25*numSets];
        }
        
        tempPtrComplex=bufferComplex;
        nCoeff=tempPtrComplex;//Length M1

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

            tempPtrComplex+=M1;
            a1=tempPtrComplex;//Length numSets (complex)

            tempPtrComplex+=numSets;
            a2=tempPtrComplex;//Length numSets (complex)

            tempPtrComplex+=numSets;
            a3=tempPtrComplex;//Length numSets (complex)

            tempPtrComplex+=numSets;
            a4=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            A1=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            A2=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            A3=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            B1=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            B2=tempPtrComplex;//Length numSets (complex)
            
            tempPtrComplex+=numSets;
            B3=tempPtrComplex;//Length numSets (complex)

            if(HessianVRet!=NULL) {
                tempPtrComplex+=numSets;
                a11=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a12=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a13=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a14=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a23=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a24=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a33=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a34=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                a44=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                A4=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                A5=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                A6=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                B4=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                B5=tempPtrComplex;//Length numSets (complex)
                
                tempPtrComplex+=numSets;
                B6=tempPtrComplex;//Length numSets (complex)
            }
        }
    } else if(algorithm==1) {//If only Legendre's algorithm is used.
        double *tempPtr;
        complex<double> *tempPtrComplex;
        
        if(gradVRet==NULL&&HessianVRet==NULL) {
            buffer=new double[2*M1+C.totalNumEl];
            bufferComplex= new complex<double>[M1+2*M1NumSets];
        }else if(gradVRet!=NULL&&HessianVRet==NULL) {
            buffer=new double[2*M1+2*C.totalNumEl];
            bufferComplex=new complex<double>[M1+6*M1NumSets];
        } else {
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex<double>[M1+12*M1NumSets];
        }

        tempPtrComplex=bufferComplex;
        nCoeff=tempPtrComplex;//Length M1

        tempPtr=buffer;
        SinVec=tempPtr;//Length M1

        tempPtr+=M1;
        CosVec=tempPtr;//Length M1

        tempPtr+=M1NumSets;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        tempPtrComplex+=M1;
        A=tempPtrComplex;//Length M1NumSets (complex)

        tempPtrComplex+=M1NumSets;
        B=tempPtrComplex;//Length M1NumSets (complex)

        if(gradVRet!=NULL||HessianVRet!=NULL) {
            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl
 
            tempPtrComplex+=M1NumSets;
            Ar=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            Br=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            ATheta=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            BTheta=tempPtrComplex;//Length M1NumSets (complex)

            if(HessianVRet!=NULL) {
                tempPtr+=C.totalNumEl;
                FuncDerivs2.clusterEls=tempPtr;//Length C.totalNumEl

                tempPtrComplex+=M1NumSets;
                AThetaTheta=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                BThetaTheta=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                Arr=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                Brr=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                AThetar=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                BThetar=tempPtrComplex;//Length M1NumSets (complex)
            }
        }
    } else {
        double *tempPtr;
        complex <double> *tempPtrComplex, *tempPtrPines;
        //If the Legendre and Pines algorithm can be used.
        
        if(HessianVRet!=NULL) {
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex<double>[M1+2*M1NumSets+max(4*M1NumSets,10*numSets)+max(6*M1NumSets,15*numSets)];
        } else if(gradVRet!=NULL) {
            //The 4*M1NumSets is memory for the Legendre algorithm.
            //The 10NumSets is memory for Pines' algorithm.
            buffer=new double[2*M1+3*C.totalNumEl];
            bufferComplex=new complex<double>[M1+2*M1NumSets+max(4*M1NumSets,10*numSets)];
        } else {
            buffer=new double[2*M1+C.totalNumEl];
            bufferComplex=new complex<double>[M1+2*M1NumSets];
        }
        
        tempPtrComplex=bufferComplex;
        nCoeff=tempPtrComplex;//Length M1

        tempPtr=buffer;
        rm=tempPtr;
        SinVec=tempPtr;//Length M1

        tempPtr+=M1;
        im=tempPtr;
        CosVec=tempPtr;//Length M1

        tempPtr+=M1;
        FuncVals.clusterEls=tempPtr;//Length C.totalNumEl
        
        tempPtrComplex+=M1;
        A=tempPtrComplex;//Length M1NumSets (complex)

        tempPtrComplex+=M1NumSets;
        B=tempPtrComplex;//Length M1NumSets (complex)

        if(gradVRet!=NULL||HessianVRet!=NULL) {
            tempPtr+=C.totalNumEl;
            FuncDerivs.clusterEls=tempPtr;//Length C.totalNumEl

            tempPtr+=C.totalNumEl;
            FuncDerivs2.clusterEls=tempPtr;// Length C.totalNumEl 

            tempPtrComplex+=M1NumSets;
            tempPtrPines=tempPtrComplex;
            Ar=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            Br=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            ATheta=tempPtrComplex;//Length M1NumSets (complex)

            tempPtrComplex+=M1NumSets;
            BTheta=tempPtrComplex;//Length M1NumSets (complex)

            a1=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            a2=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            a3=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            a4=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            A1=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            A2=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            A3=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            B1=tempPtrPines;//Length numSets (complex)

            tempPtrPines+=numSets;
            B2=tempPtrPines;//Length numSets (complex)

            tempPtr+=numSets;
            B3=tempPtrPines;//Length numSets (complex)
            
            if(HessianVRet!=NULL) {
                tempPtrComplex+=M1NumSets;
                AThetaTheta=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                BThetaTheta=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                Arr=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                Brr=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                AThetar=tempPtrComplex;//Length M1NumSets (complex)

                tempPtrComplex+=M1NumSets;
                BThetar=tempPtrComplex;//Length M1NumSets (complex)

                //Pines' Part
                tempPtrPines+=numSets;
                a11=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a12=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a13=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a14=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a23=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a24=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a33=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a34=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                a44=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                A4=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                A5=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                A6=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                B4=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                B5=tempPtrPines;//Length numSets (complex)

                tempPtrPines+=numSets;
                B6=tempPtrPines;//Length numSets (complex)
            }
        }
    }

    nCoeff[0]=1;
    
    rPrev=std::numeric_limits<double>::infinity();
    thetaPrev=std::numeric_limits<double>::infinity();
    for(curPoint=0;curPoint<numPoints;curPoint++) {
        double pointCur[3];
        double r, lambda, thetaCur;
        size_t n,m, curSet;
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
                fill(A,A+M1NumSets,complex<double>(0.0,0.0));
                fill(B,B+M1NumSets,complex<double>(0.0,0.0));

                //Compute the X coefficients for the sum
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t curSetOffset=curSet*M1;
                    for(m=0;m<=M;m++) {
                        const size_t k=m+curSetOffset;
                        for(n=m;n<=M;n++) {
                            A[k]+=nCoeff[n]*C(n,m,curSet)*FuncVals[n][m];
                            B[k]+=nCoeff[n]*S(n,m,curSet)*FuncVals[n][m];
                        }
                    }
                }
                
                //If additional terms should be computed so a gradient or
                //Hessian can be computed.
                if(gradVRet!=NULL||HessianVRet!=NULL) {
                    fill(Ar,Ar+M1NumSets,complex<double>(0.0,0.0));
                    fill(Br,Br+M1NumSets,complex<double>(0.0,0.0));
                    fill(ATheta,ATheta+M1NumSets,complex<double>(0.0,0.0));
                    fill(BTheta,BTheta+M1NumSets,complex<double>(0.0,0.0));
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        const size_t curSetOffset=curSet*M1;
                        
                        mf=0;
                        for(m=0;m<=M;m++) {
                            const size_t k=m+curSetOffset;
                            nf=mf;
                            for(n=m;n<=M;n++) {
                                const complex<double> CScal=nCoeff[n]*C(n,m,curSet);
                                const complex<double> SScal=nCoeff[n]*S(n,m,curSet);
                                double curVal;

                                curVal=FuncVals[n][m];
                                Ar[k]+=(nf+1)*CScal*curVal;
                                Br[k]+=(nf+1)*SScal*curVal;

                                curVal=FuncDerivs[n][m];
                                ATheta[k]+=CScal*curVal;
                                BTheta[k]+=SScal*curVal;

                                nf++;
                            }

                            mf++;
                        }
                    }
                    
                    if(HessianVRet!=NULL) {
                        fill(Arr,Arr+M1NumSets,complex<double>(0.0,0.0));
                        fill(Brr,Brr+M1NumSets,complex<double>(0.0,0.0));
                        fill(AThetar,AThetar+M1NumSets,complex<double>(0.0,0.0));
                        fill(BThetar,BThetar+M1NumSets,complex<double>(0.0,0.0));
                        fill(AThetaTheta,AThetaTheta+M1NumSets,complex<double>(0.0,0.0));
                        fill(BThetaTheta,BThetaTheta+M1NumSets,complex<double>(0.0,0.0));
                        
                        for(curSet=0;curSet<numSets;curSet++) {
                            const size_t curSetOffset=curSet*M1;
                            mf=0;
                            for(m=0;m<=M;m++) {
                                const size_t k=m+curSetOffset;
                                
                                nf=mf;
                                for(n=m;n<=M;n++) {
                                    const complex<double> CScal=nCoeff[n]*C(n,m,curSet);
                                    const complex<double> SScal=nCoeff[n]*S(n,m,curSet);
                                    double curVal;

                                    curVal=FuncVals[n][m];
                                    //From Table 5, with the correction from the
                                    //erratum.
                                    Arr[k]+=(nf+1)*(nf+2)*CScal*curVal;
                                    Brr[k]+=(nf+1)*(nf+2)*SScal*curVal;;

                                    curVal=FuncDerivs[n][m];
                                    //From Table 5
                                    AThetar[k]+=(nf+1)*CScal*curVal;
                                    BThetar[k]+=(nf+1)*SScal*curVal;

                                    curVal=FuncDerivs2[n][m];
                                    //From Table 5
                                    AThetaTheta[k]+=CScal*curVal;
                                    BThetaTheta[k]+=SScal*curVal;

                                    nf++;
                                }

                                mf++;
                            }
                        }
                    }
                }
            }
            
            //Use Horner's method to compute V.
            {
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t curSetOffset=M1*curSet;
                    const size_t idx=curSet+numSets*curPoint;
                    
                    complex <double> V(0,0);
                    m=M+1;
                    do {
                        m--;
                        const size_t k=m+curSetOffset;

                        V=V*u+A[k]*CosVec[m]+B[k]*SinVec[m];
                    } while(m>0);

                    //Multiply by the constant in front of the sum and get rid
                    //of the scale factor.
                    V*=crScal;
                    VRet[idx]=V;
                }
            }

            //Compute the gradient, if it is desired.
            if(gradVRet!=NULL||HessianVRet!=NULL) {
                double J[9];
                double H[27];
                double drdx;
                double drdy;
                double drdz;
                double dLambdadx;
                double dLambdady;
                double dLambdadz;
                double dPhidx;
                double dPhidy;
                double dPhidz;
                double drdxdx;
                double drdydy;
                double drdzdz;
                double drdxdy;
                double drdxdz;
                double drdydz;
                double dLambdadxdx;
                double dLambdadydy;
                double dLambdadzdz;
                double dLambdadxdy;
                double dLambdadxdz;
                double dLambdadydz;
                double dPhidxdx;
                double dPhidydy;
                double dPhidzdz;
                double dPhidxdy;
                double dPhidxdz;
                double dPhidydz;
                
                if(spherDerivs==false) {
                    calcSpherConvJacobCPP(J, pointCur,0);
                    
                    if(HessianVRet!=NULL) {    
                        calcSpherConvHessianCPP(H,pointCur,0,true);
                        
                        drdx=J[0];
                        drdy=J[3];
                        drdz=J[6];
                        dLambdadx=J[1];
                        dLambdady=J[4];
                        dLambdadz=J[7];
                        dPhidx=J[2];
                        dPhidy=J[5];
                        dPhidz=J[8];

                        //Commented indices on H are what would be used in
                        //Matlab.
                        drdxdx=H[0];//H(1,1,1)
                        drdydy=H[4];//H(2,2,1)
                        drdzdz=H[8];//H(3,3,1)
                        drdxdy=H[3];//H(1,2,1)
                        drdxdz=H[6];//H(1,3,1)
                        drdydz=H[7];//H(2,3,1)

                        dLambdadxdx=H[9];//H(1,1,2)
                        dLambdadydy=H[13];//H(2,2,2)
                        dLambdadzdz=H[17];//H(3,3,2)
                        dLambdadxdy=H[12];//H(1,2,2)
                        dLambdadxdz=H[15];//H(1,3,2)
                        dLambdadydz=H[16];//H(2,3,2)

                        dPhidxdx=H[18];//H(1,1,3)
                        dPhidydy=H[22];//H(2,2,3)
                        dPhidzdz=H[26];//H(3,3,3)
                        dPhidxdy=H[21];//H(1,2,3)
                        dPhidxdz=H[24];//H(1,3,3)
                        dPhidydz=H[25];//H(2,3,3)
                    }
                }
                
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t curSetOffset=M1*curSet;
                    
                    complex<double> dVdr(0,0);
                    complex<double> dVdLambda(0,0);
                    complex<double> dVdTheta(0,0);

                    //Use Horner's method to compute all dV values.
                    m=M+1;
                    mf=static_cast<double>(m);
                    do {
                        m--;
                        mf--;
                        const size_t k=m+curSetOffset;

                        dVdr=dVdr*u-(Ar[k]*CosVec[m]+Br[k]*SinVec[m]);
                        dVdLambda=dVdLambda*u-mf*(A[k]*SinVec[m]-B[k]*CosVec[m]);
                        dVdTheta=dVdTheta*u+(ATheta[k]*CosVec[m]+BTheta[k]*SinVec[m]);
                    } while(m>0);            

                    dVdr=crScal*dVdr/r;
                    dVdLambda=crScal*dVdLambda;
                    //The minus sign adjusts for the coordinate system change.
                    dVdTheta=-crScal*dVdTheta;

                    if(spherDerivs&&gradVRet!=NULL) {
                        const size_t idx=3*curSet+3*numSets*curPoint;
                        //Now, multiply the transpose of the Jacobian Matrix by the
                        //vector of [dVdr;dVdLambda;dVdTheta]
                        gradVRet[0+idx]=dVdr;
                        gradVRet[1+idx]=dVdLambda;
                        gradVRet[2+idx]=dVdTheta;
                    } else {
                        if(gradVRet!=NULL) {
                            const size_t idx=3*curSet+3*numSets*curPoint;
                            
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
                            const size_t k=m+curSetOffset;
                            
                            d2VdLambdadLambda=d2VdLambdadLambda*u-mf*mf*(A[k]*CosVec[m]+B[k]*SinVec[m]);
                            d2VdLambdadTheta=d2VdLambdadTheta*u+mf*(ATheta[k]*SinVec[m]-BTheta[k]*CosVec[m]);
                            d2VdrdLambda=d2VdrdLambda*u+mf*(Ar[k]*SinVec[m]-Br[k]*CosVec[m]);
                            d2VdThetadTheta=d2VdThetadTheta*u+(AThetaTheta[k]*CosVec[m]+BThetaTheta[k]*SinVec[m]);
                            d2VdrdTheta=d2VdrdTheta*u-(AThetar[k]*CosVec[m]+BThetar[k]*SinVec[m]);
                            d2Vdrdr=d2Vdrdr*u+(Arr[k]*CosVec[m]+Brr[k]*SinVec[m]);
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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
            for(curSet=0;curSet<numSets;curSet++) {
                const size_t idx=curSet+numSets*curPoint;
                
                complex <double> V(0,0);
                for(n=0;n<=M;n++) {
                    complex<double> innerTerm(0,0);
                    for(m=0;m<=n;m++) {
                        innerTerm+=(C(n,m,curSet)*rm[m]+S(n,m,curSet)*im[m])*FuncVals[n][m];
                    }

                    V+=nCoeff[n]*innerTerm;
                }

                V*=crScal;
                VRet[idx]=V;
            }

            //If only the gradient is desired as the next output       
            if(gradVRet!=NULL&&HessianVRet==NULL) {
                fill(a1,a1+numSets,complex<double>(0.0,0.0));
                fill(a2,a2+numSets,complex<double>(0.0,0.0));
                fill(a3,a3+numSets,complex<double>(0.0,0.0));
                fill(a4,a4+numSets,complex<double>(0.0,0.0));

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    fill(A1,A1+numSets,complex<double>(0.0,0.0));
                    fill(A2,A2+numSets,complex<double>(0.0,0.0));
                    fill(A3,A3+numSets,complex<double>(0.0,0.0));
                    fill(B1,B1+numSets,complex<double>(0.0,0.0));
                    fill(B2,B2+numSets,complex<double>(0.0,0.0));
                    fill(B3,B3+numSets,complex<double>(0.0,0.0));
                    
                    //Compute the lumped coefficients for Pine's method from
                    //Table 13 for the current m.
                    nf=mf;
                    for(n=m;n<=M;n++) {
                        const double HVal=FuncVals[n][m];
                        const double dHVal=FuncDerivs[n][m];
                        //The expressions for Lmn, is from Table 14
                        const double Lmn=(nf+mf+1)*HVal+u*dHVal;
                        
                        for(curSet=0;curSet<numSets;curSet++) {
                            const complex<double> rhoC=nCoeff[n]*C(n,m,curSet);
                            const complex<double> rhoS=nCoeff[n]*S(n,m,curSet);

                            A1[curSet]+=rhoC*HVal;
                            A2[curSet]+=rhoC*dHVal;
                            A3[curSet]+=rhoC*Lmn;

                            B1[curSet]+=rhoS*HVal;
                            B2[curSet]+=rhoS*dHVal;
                            B3[curSet]+=rhoS*Lmn;
                        }
                        nf++;
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        if(M>=1) {
                            a1[curSet]+=mf*(A1[curSet]*rm[m-1]+B1[curSet]*im[m-1]);
                            a2[curSet]+=mf*(B1[curSet]*rm[m-1]-A1[curSet]*im[m-1]);
                        }
                        a3[curSet]+=(A2[curSet]*rm[m]+B2[curSet]*im[m]);
                        a4[curSet]-=(A3[curSet]*rm[m]+B3[curSet]*im[m]);
                    }
                    
                    mf++;
                }
                
                {
                    double J[9];
                    if(spherDerivs) {
                        calcSpherInvJacobCPP(J,pointCur,0);
                    }
                
                    for(curSet=0;curSet<numSets;curSet++) {                
                        a1[curSet]/=r;
                        a2[curSet]/=r;
                        a3[curSet]/=r;
                        a4[curSet]/=r;

                        const complex<double> dVdx=crScal*(a1[curSet]+s*a4[curSet]);
                        const complex<double> dVdy=crScal*(a2[curSet]+t*a4[curSet]);
                        const complex<double> dVdz=crScal*(a3[curSet]+u*a4[curSet]);

                        if(spherDerivs) {
                            const size_t idx=3*curSet+3*numSets*curPoint;

                            //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                            gradVRet[0+idx]=dVdx*J[0]+dVdy*J[1]+dVdz*J[2];
                            gradVRet[1+idx]=dVdx*J[3]+dVdy*J[4]+dVdz*J[5];
                            gradVRet[2+idx]=dVdx*J[6]+dVdy*J[7]+dVdz*J[8];
                        } else {//If a gradient in Cartesian coordinates is
                                //desired.
                            const size_t idx=3*curSet+3*numSets*curPoint;

                            gradVRet[0+idx]=dVdx;
                            gradVRet[1+idx]=dVdy;
                            gradVRet[2+idx]=dVdz;
                        }
                    }
                
                }
            } else if(HessianVRet!=NULL) {
            //Compute the gradient and Hessian
                fill(a1,a1+numSets,complex<double>(0.0,0.0));
                fill(a2,a2+numSets,complex<double>(0.0,0.0));
                fill(a3,a3+numSets,complex<double>(0.0,0.0));
                fill(a4,a4+numSets,complex<double>(0.0,0.0));
                fill(a11,a11+numSets,complex<double>(0.0,0.0));
                fill(a12,a12+numSets,complex<double>(0.0,0.0));
                fill(a13,a13+numSets,complex<double>(0.0,0.0));
                fill(a14,a14+numSets,complex<double>(0.0,0.0));
                fill(a23,a23+numSets,complex<double>(0.0,0.0));
                fill(a24,a24+numSets,complex<double>(0.0,0.0));
                fill(a33,a33+numSets,complex<double>(0.0,0.0));
                fill(a34,a34+numSets,complex<double>(0.0,0.0));
                fill(a44,a44+numSets,complex<double>(0.0,0.0));
                
                const double r2=r*r;
                const double s2=s*s;
                const double u2=u*u;
                const double t2=t*t;

                //The equations in these loops are from Table 10.
                mf=0.0;
                for(m=0;m<=M;m++) {
                    fill(A1,A1+numSets,complex<double>(0.0,0.0));
                    fill(A2,A2+numSets,complex<double>(0.0,0.0));
                    fill(A3,A3+numSets,complex<double>(0.0,0.0));
                    fill(A4,A4+numSets,complex<double>(0.0,0.0));
                    fill(A5,A5+numSets,complex<double>(0.0,0.0));
                    fill(A6,A6+numSets,complex<double>(0.0,0.0));
                    fill(B1,B1+numSets,complex<double>(0.0,0.0));
                    fill(B2,B2+numSets,complex<double>(0.0,0.0));
                    fill(B3,B3+numSets,complex<double>(0.0,0.0));
                    fill(B4,B4+numSets,complex<double>(0.0,0.0));
                    fill(B5,B5+numSets,complex<double>(0.0,0.0));
                    fill(B6,B6+numSets,complex<double>(0.0,0.0));
                    
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
                        
                        for(curSet=0;curSet<numSets;curSet++) {
                            const complex <double> rhoC=nCoeff[n]*C(n,m,curSet);
                            const complex <double> rhoS=nCoeff[n]*S(n,m,curSet);

                            A1[curSet]+=rhoC*HVal;
                            A2[curSet]+=rhoC*dHVal;
                            A3[curSet]+=rhoC*Lmn;
                            A4[curSet]+=rhoC*d2HVal;
                            A5[curSet]+=rhoC*dLmn;
                            A6[curSet]+=rhoC*Omn;

                            B1[curSet]+=rhoS*HVal;
                            B2[curSet]+=rhoS*dHVal;
                            B3[curSet]+=rhoS*Lmn;
                            B4[curSet]+=rhoS*d2HVal;
                            B5[curSet]+=rhoS*dLmn;
                            B6[curSet]+=rhoS*Omn;
                        }
                        
                        nf++;
                    }
                    
                    for(curSet=0;curSet<numSets;curSet++) {
                        if(m>=1) {
                            a1[curSet]+=mf*(A1[curSet]*rm[m-1]+B1[curSet]*im[m-1]);
                            a2[curSet]+=mf*(B1[curSet]*rm[m-1]-A1[curSet]*im[m-1]);

                            a13[curSet]+=mf*(A2[curSet]*rm[m-1]+B2[curSet]*im[m-1]);
                            a14[curSet]-=mf*(A3[curSet]*rm[m-1]+B3[curSet]*im[m-1]);
                            a23[curSet]+=mf*(B2[curSet]*rm[m-1]-A2[curSet]*im[m-1]);
                            a24[curSet]-=mf*(B3[curSet]*rm[m-1]-A3[curSet]*im[m-1]);

                            if(m>=2) {
                                a11[curSet]+=mf*(mf-1)*(A1[curSet]*rm[m-2]+B1[curSet]*im[m-2]);
                                a12[curSet]+=mf*(mf-1)*(B1[curSet]*rm[m-2]-A1[curSet]*im[m-2]);
                            }
                        }
                        a3[curSet]+=(A2[curSet]*rm[m]+B2[curSet]*im[m]);
                        a4[curSet]-=(A3[curSet]*rm[m]+B3[curSet]*im[m]);
                        
                        a33[curSet]+=(A4[curSet]*rm[m]+B4[curSet]*im[m]);
                        a34[curSet]-=(A5[curSet]*rm[m]+B5[curSet]*im[m]);
                        a44[curSet]+=(A6[curSet]*rm[m]+B6[curSet]*im[m]);
                    }
                    mf++;
                }
                
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
                    double d2xdrdr;
                    double d2xdAzdAz;
                    double d2xdEldEl;
                    double d2xdrdAz;
                    double d2xdrdEl;
                    double d2xdAzdEl;
                    double d2ydrdr;
                    double d2ydAzdAz;
                    double d2ydEldEl;
                    double d2ydrdAz;
                    double d2ydrdEl;
                    double d2ydAzdEl;
                    double d2zdrdr;
                    double d2zdAzdAz;
                    double d2zdEldEl;
                    double d2zdrdAz;
                    double d2zdrdEl;
                    double d2zdAzdEl;
                    
                    if(spherDerivs) {
                        double J[9];
                        double H[27];

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
                        
                        calcSpherInvHessianCPP(H,pointCur,0);

                        d2xdrdr=H[0];//H(1,1,1);
                        d2xdAzdAz=H[4];//H(2,2,1);
                        d2xdEldEl=H[8];//H(3,3,1);
                        d2xdrdAz=H[3];//H(1,2,1);
                        d2xdrdEl=H[6];//H(1,3,1);
                        d2xdAzdEl=H[7];//H(2,3,1);

                        d2ydrdr=H[9];//H(1,1,2);
                        d2ydAzdAz=H[13];//H(2,2,2);
                        d2ydEldEl=H[17];//H(3,3,2);
                        d2ydrdAz=H[12];//H(1,2,2);
                        d2ydrdEl=H[15];//H(1,3,2);
                        d2ydAzdEl=H[16];//H(2,3,2);

                        d2zdrdr=H[18];//H(1,1,3);
                        d2zdAzdAz=H[22];//H(2,2,3);
                        d2zdEldEl=H[26];//H(3,3,3);
                        d2zdrdAz=H[21];//H(1,2,3);
                        d2zdrdEl=H[24];//H(1,3,3);
                        d2zdAzdEl=H[25];//H(2,3,3);
                    }

                    for(curSet=0;curSet<numSets;curSet++) {               
                        a1[curSet]/=r;
                        a2[curSet]/=r;
                        a3[curSet]/=r;
                        a4[curSet]/=r;
                        a11[curSet]/=r2;
                        a12[curSet]/=r2;
                        a13[curSet]/=r2;
                        a14[curSet]/=r2;
                        a23[curSet]/=r2;
                        a24[curSet]/=r2;
                        a33[curSet]/=r2;
                        a34[curSet]/=r2;
                        a44[curSet]/=r2;

                        const complex<double> dVdx=crScal*(a1[curSet]+s*a4[curSet]);
                        const complex<double> dVdy=crScal*(a2[curSet]+t*a4[curSet]);
                        const complex<double> dVdz=crScal*(a3[curSet]+u*a4[curSet]);
                        const complex<double> d2Vdxdx=crScal*(a11[curSet]+2.0*s*a14[curSet]+a4[curSet]/r+s2*a44[curSet]-s2*a4[curSet]/r);
                        const complex<double> d2Vdydy=crScal*(-a11[curSet]+2.0*t*a24[curSet]+a4[curSet]/r+t2*a44[curSet]-t2*a4[curSet]/r);
                        const complex<double> d2Vdzdz=crScal*(a33[curSet]+2.0*u*a34[curSet]+a4[curSet]/r+u2*a44[curSet]-u2*a4[curSet]/r);
                        const complex<double> d2Vdxdy=crScal*(a12[curSet]+s*t*a44[curSet]+s*a24[curSet]+t*a14[curSet]-s*t*a4[curSet]/r);
                        const complex<double> d2Vdxdz=crScal*(a13[curSet]+s*u*a44[curSet]+s*a34[curSet]+u*a14[curSet]-s*u*a4[curSet]/r);
                        const complex<double> d2Vdydz=crScal*(a23[curSet]+t*u*a44[curSet]+t*a34[curSet]+u*a24[curSet]-t*u*a4[curSet]/r);

                        //Compute the gradient, if desired.
                        if(gradVRet!=NULL) {
                            const size_t idx=3*curSet+3*numSets*curPoint;

                            if(spherDerivs) {
                                //gradV(:,curPoint)=J'*[dVdx;dVdy;dVdz];
                                gradVRet[0+idx]=dVdx*dxdr+dVdy*dydr+dVdz*dzdr;
                                gradVRet[1+idx]=dVdx*dxdAz+dVdy*dydAz+dVdz*dzdAz;
                                gradVRet[2+idx]=dVdx*dxdEl+dVdy*dydEl+dVdz*dzdEl;
                            } else {
                                gradVRet[0+idx]=dVdx;
                                gradVRet[1+idx]=dVdy;
                                gradVRet[2+idx]=dVdz;
                            }
                        }

                        //Compute the Hessian.
                        if(spherDerivs) {
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
                            const size_t idx=9*curSet+9*numSets*curPoint;

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
        }
        
        if(spherDerivs&&systemType==2) {
            //Flip signs of the elevation terms reflecting the
            //difference in definition between systems 0 and 2.
            
            if(gradVRet!=NULL) {
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t idx=3*curSet+3*numSets*curPoint;
                    
                    gradVRet[2+idx]=-gradVRet[2+idx];
                }
            }

            if(HessianVRet!=NULL) {
                for(curSet=0;curSet<numSets;curSet++) {
                    const size_t idx=9*curSet+9*numSets*curPoint;
                    
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
