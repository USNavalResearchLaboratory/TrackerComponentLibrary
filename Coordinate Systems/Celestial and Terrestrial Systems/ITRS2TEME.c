/*ITRS2TEME Convert from the International Terrestrial Reference
 *          System (ITRS) into the True Equator Mean Equinox (TEME) of date 
 *          coordinate system. The TEME system is non-standard and is
 *          generally only used in the Specialized General Perturbations 4
 *          (SGP4) orbital propagation algorithm. Note that the velocity
 *          conversion does not include the (small) centrifugal effect of
 *          polar motion.
 *
 *INPUTS: x The NXnumVec collection of vectors to convert. N can be 3, or
 *          6. If the vectors are 3D, then they are position. 6D vectors
 *          are assumed to be position and velocity, whereby the angular
 *          velocity of the Earth's rotation is taken into account using a
 *          non-relativistic formula.
 * Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
 *          The units of the date are days. The full date is the sum of
 *          both terms. The date is broken into two parts to provide more
 *          bits of precision. It does not matter how the date is split.
 *  deltaTTUT1 An optional parameter specifying the difference between TT
 *          and UT1 in seconds. This information can be obtained from
 * http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
 *          or 
 * http://www.usno.navy.mil/USNO/earth-orientation/eo-products
 *          If this parameter is omitted or if an empty matrix is passed,
 *          then the value provided by the function getEOP will be used
 *          instead.
 *     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
 *          including the effects of tides and librations. If this
 *          parameter is omitted or an empty matrix is passed the value
 *          from the function getEOP will be used.
 *      LOD The difference between the length of the day using terrestrial
 *          time, international atomic time, or UTC without leap seconds
 *          and the length of the day in UT1. This is an instantaneous
 *          parameter (in seconds) proportional to the rotation rate of the
 *          Earth. This is only needed if more than just position
 *          components are being converted.
 *
 *The conversion from the TEME to the pseudo-Earth-Fixed (PEF) coordinate
 *system is described in [1] and the relationship between the ITRS and the
 *PEF is described in [2].
 *
 *The velocity transformation deals with the instantaneous rotational
 *velocity of the Earth using a simple Newtonian velocity addition.
 *Basically, the axis of rotation in the Pseudo-Ears-Fixed (PEF) frame is
 *the z-axis (The PEF is akin to a less-accurate version of the TIRS).
 *The rotation rate in that system is Constants.IERSMeanEarthRotationRate
 *adjusted using the Length-of-Day (LOD) Earth Orientation Parameter (EOP).
 *Thus, in the PEF, the angular velocity vector is [0;0;omega], where omega
 *is the angular velocity accounting for the LOD EOP. Consequently, one
 *account for rotation by transforming from the ITRS to the PEF,
 *adding the cross product of Omega with the position in the PEF, and
 *then converting to the TEME. This is a simple Newtonian conversion.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[vec,rotMat]=ITRS2TEME(x,Jul1,Jul2);
 *or if more parameters are known,
 *[vec,rotMat]=ITRS2TEME(x,Jul1,Jul2,deltaTTUT1,xpyp,LOD);
 *
 *REFERENCES:
 *[1] D. A. Vallado, P. Crawford, R. Hujsak, and T. Kelso, "Implementing
 *    the revised SGP4 in STK," in Proceedings of the AGI User Exchange,
 *    Washington, DC, 17?18 Oct. 2006, slides. [Online].
 *    Available: http://www.agi.com/downloads/events/2006-agi-user-exchange/8_revised_sgp4_vallado2.pdf
 *[2] D. A. Vallado, J. H. Seago, and P. K. Seidelmann, "Implementation
 *    issues surrounding the new IAU reference systems for astrodynamics,"
 *    in Proceedings of the 16th AAS/AIAA Space Flight Mechanics
 *    Conference, Tampa, FL, 22?26 Jan. 2006. [Online].
 *    Available: http://www.centerforspace.com/downloads/files/pubs/AAS-06-134.pdf
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numRow,numVec;
    mxArray *retMat;
    double *xVec, *retData;
    double TT1, TT2, UT11, UT12;
    //The if-statements below should properly initialize all of the EOP.
    //The following initializations to zero are to suppress warnings when
    //compiling with -Wconditional-uninitialized.
    double xp=0;
    double yp=0;
    double deltaT=0;
    double LOD=0;
    double ITRS2TEME[3][3];
    double PEF2TEME[3][3];
    double WInv[3][3];//The inverse polar motion matrix to go from ITRS to PEF.
    double Omega[3];//The angular velocity vector for the Earth's rotation.

    if(nrhs<3||nrhs>6){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
 
    checkRealDoubleArray(prhs[0]);
    
    numRow=mxGetM(prhs[0]);
    numVec=mxGetN(prhs[0]);
    
    if(!(numRow==3||numRow==6)) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }
    
    xVec=mxGetDoubles(prhs[0]);
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);
    
    //If some values from the function getEOP will be needed
    if(nrhs<6||mxIsEmpty(prhs[3])||mxIsEmpty(prhs[4])||mxIsEmpty(prhs[5])) {
        mxArray *retVals[5];
        double *xpyp;
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        int retVal;
        
        //Get the time in UTC to look up the parameters by going to TAI and
        //then UTC.
        retVal=iauTttai(TT1, TT2, &JulUTC[0], &JulUTC[1]);
        if(retVal!=0) {
            mexErrMsgTxt("An error occurred computing TAI.");
        }
        retVal=iauTaiutc(JulUTC[0], JulUTC[1], &JulUTC[0], &JulUTC[1]);
        switch(retVal){
            case 1:
                mexWarnMsgTxt("Dubious Date entered.");
                break;
            case -1:
                mexErrMsgTxt("Unacceptable date entered");
                break;
            default:
                break;
        }
        
        JulUTCMATLAB[0]=doubleMat2Matlab(&JulUTC[0],1,1);
        JulUTCMATLAB[1]=doubleMat2Matlab(&JulUTC[1],1,1);

        //Get the Earth orientation parameters for the given date.
        mexCallMATLAB(5,retVals,2,JulUTCMATLAB,"getEOP");
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
        
        checkRealDoubleArray(retVals[0]);
        checkRealDoubleArray(retVals[1]);
        if(mxGetM(retVals[0])!=2||mxGetN(retVals[0])!=1||mxGetM(retVals[1])!=2||mxGetN(retVals[1])!=1) {
            mxDestroyArray(retVals[0]);
            mxDestroyArray(retVals[1]);
            mxDestroyArray(retVals[2]);
            mxDestroyArray(retVals[3]);
            mxDestroyArray(retVals[4]);
            mexErrMsgTxt("Error using the getEOP function.");
            return;
        }
        
        xpyp=mxGetDoubles(retVals[0]);
        xp=xpyp[0];
        yp=xpyp[1];
        //The celestial pole offsets are not used.
        
        //This is TT-UT1
        deltaT=getDoubleFromMatlab(retVals[3]);
        LOD=getDoubleFromMatlab(retVals[4]);
        //Free the returned arrays.
        mxDestroyArray(retVals[0]);
        mxDestroyArray(retVals[1]);
        mxDestroyArray(retVals[2]);
        mxDestroyArray(retVals[3]);
        mxDestroyArray(retVals[4]);
    }
    
    //If deltaT=TT-UT1 is given
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        deltaT=getDoubleFromMatlab(prhs[3]);
    }
    
    //Obtain UT1 from terestrial time and deltaT.
    iauTtut1(TT1, TT2, deltaT, &UT11, &UT12);
    
    //Get polar motion values, if given.
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        size_t dim1, dim2;
        
        checkRealDoubleArray(prhs[4]);
        dim1=mxGetM(prhs[4]);
        dim2=mxGetN(prhs[4]);
        
        if((dim1==2&&dim2==1)||(dim1==1&&dim2==2)) {
            double *xpyp=mxGetDoubles(prhs[4]);
        
            xp=xpyp[0];
            yp=xpyp[1];
        } else {
            mexErrMsgTxt("The celestial pole offsets have the wrong dimensionality.");
            return;
        }
    }
    
    //If LOD is given
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        LOD=getDoubleFromMatlab(prhs[5]);
    }

    {
     double GMST2000;
     double TEME2PEF[3][3];
     double TEME2ITRS[3][3];
     double W[3][3];
     double omega;
    
     //Get Greenwhich mean sidereal time under the IAU's 2000 model. This
     //is given in radians and will be used to build a rotation matrix to
     //rotate into the PEF system.
     GMST2000=iauGmst00(UT11, UT12,TT1, TT2);
     {
         double cosGMST,sinGMST;
         cosGMST=cos(GMST2000);
         sinGMST=sin(GMST2000);
         //Build the rotation matrix to rotate by GMST about the z-axis. This
         //will put the position vector in the PEF system.
         TEME2PEF[0][0]=cosGMST;
         TEME2PEF[0][1]=sinGMST;
         TEME2PEF[0][2]=0;
         TEME2PEF[1][0]=-sinGMST;
         TEME2PEF[1][1]=cosGMST;
         TEME2PEF[1][2]=0;
         TEME2PEF[2][0]=0;
         TEME2PEF[2][1]=0;
         TEME2PEF[2][2]=1.0;
     }
     //The inverse rotation is just the transpose
     iauTr(TEME2PEF, PEF2TEME);
     //To go from PEF to ITRS, we need to build the polar motion matrix
     //using the IAU's 2000 conventions.
     {
         double cosXp,sinXp,cosYp,sinYp;
         cosXp=cos(xp);
         sinXp=sin(xp);
         cosYp=cos(yp);
         sinYp=sin(yp);
         W[0][0]=cosXp;
         W[0][1]=sinXp*sinYp;
         W[0][2]=sinXp*cosYp;
         W[1][0]=0;
         W[1][1]=cosYp;
         W[1][2]=-sinYp;
         W[2][0]=-sinXp;
         W[2][1]=cosXp*sinXp;
         W[2][2]=cosXp*cosYp;
     }
     //The inverse rotation is just the transpose
     iauTr(W, WInv);
     
     //The total rotation matrix is thus the product of the two rotations.
     //TEME2ITRS=W*TEME2PEF;
     iauRxr(W, TEME2PEF, TEME2ITRS);
     //We want the inverse rotation
     iauTr(TEME2ITRS, ITRS2TEME);
     //The angular velocity vector of the Earth in the TIRS in radians.
     omega=getScalarMatlabClassConst("Constants","IERSMeanEarthRotationRate");
     //Adjust for LOD
     omega=omega*(1-LOD/86400.0);//86400.0 seconds are in a TT day.
     Omega[0]=0;
     Omega[1]=0;
     Omega[2]=omega;     
    }
    
    //Allocate space for the return vectors.
    retMat=mxCreateDoubleMatrix(numRow,numVec,mxREAL);
    retData=mxGetDoubles(retMat);
    
    {
        size_t curVec;
        
        for(curVec=0;curVec<numVec;curVec++) {
            //Multiply the position vector with the rotation matrix.
            iauRxp(ITRS2TEME, xVec+numRow*curVec, retData+numRow*curVec);
            //If a velocity vector was given.
            if(numRow>3) {
                double *posITRS=xVec+numRow*curVec;
                double *velITRS=xVec+numRow*curVec+3;//Velocity in TEME
                double posPEF[3];
                double velPEF[3];
                double *retDataVel=retData+numRow*curVec+3;
                double rotVel[3];
                //If a velocity was provided with the position, first
                //convert to PEF coordinates, then account for the rotation
                //of the Earth, then rotate into TEME coordinates.
                
                //Convert velocity from ITRS to PEF.
                iauRxp(WInv, velITRS, velPEF);
                //Convert position from ITRS to PEF
                iauRxp(WInv, posITRS, posPEF);

                //Evaluate the cross product for the angular velocity due
                //to the Earth's rotation.
                iauPxp(Omega, posPEF, rotVel);

                //Add the instantaneous velocity due to rotation.
                iauPpp(velPEF, rotVel, retDataVel);

                //Rotate from the PEF into the TEME
                iauRxp(PEF2TEME, retDataVel, retDataVel);
            }
        }
    }
    
    plhs[0]=retMat;
    
    if(nlhs>1) {
        double *elPtr;
        size_t i,j;
        
        plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
        elPtr=mxGetDoubles(plhs[1]);
        
        for (i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                elPtr[i+3*j]=ITRS2TEME[i][j];
            }
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
