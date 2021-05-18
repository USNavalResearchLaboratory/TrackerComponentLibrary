/**CIRS2TIRS Convert vectors of position and possibly velocity from the
 *           Celestial Intermediate Reference System (CIRS) to the 
 *           Terrestrial Intermediate Reference System (TIRS).
 *
 *INPUTS: x The NXnumVec collection of vectors to convert. N can be 3, or
 *          6. If the vectors are 3D, then they can be either position or
 *          velocity. 6D vectors are assumed to be position and velocity,
 *          whereby the angular velocity of the Earth's rotation is taken
 *          into account using a non-relativistic formula.
 * Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
 *          The units of the date are days. The full date is the sum of
 *          both terms. The date is broken into two parts to provide more
 *          bits of precision. It does not matter how the date is split.
 * deltaTTUT1 An optional parameter specifying the difference between TT
 *          and UT1 in seconds. This information can be obtained from
 * http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
 *          or 
 * http://www.usno.navy.mil/USNO/earth-orientation/eo-products
 *          If this parameter is omitted or if an empty matrix is passed,
 *          then the value provided by the function getEOP will be used
 *          instead.
 *      LOD The difference between the length of the day using terrestrial
 *          time, international atomic time, or UTC without leap seconds
 *          and the length of the day in UT1. This is an instantaneous
 *          parameter (in seconds) proportional to the rotation rate of the
 *          Earth. This is only needed if more than just position
 *          components are being converted.
 *
 *OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from CIRS
 *             coordinates to TIRS coordinates.
 *      rotMat The 3X3 rotation matrix used for the conversion of the
 *             positions.
 *
 *The conversion functions from the International Astronomical Union's
 *(IAU) Standard's of Fundamental Astronomy library are put together to get
 *the necessary rotation matrix for the position.
 *
 *The velocity transformation deals with the instantaneous rotational
 *velocity of the Earth using a simple Newtonian velocity addition.
 *Basically, the axis of rotation in the Terrestrial Intermediate Reference
 *System TIRS is the z-axis. The rotation rate in that system is
 *Constants.IERSMeanEarthRotationRate adjusted using the Length-of-Day
 *(LOD) Earth Orientation Parameter (EOP). Thus, in the TIRS, the angular
 *velocity vector is [0;0;omega], where omega is the angular velocity
 *accounting for the LOD EOP. Consequently, one accounts for rotation by
 *transforming from the GCRS to the TIRS, and subtracting the cross product
 *of Omega with the position in the TIRS.
 *This is a simple Newtonian conversion.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[vec,rotMat]=CIRS2TIRS(x,Jul1,Jul2);
 *or if more parameters are known,
 *[vec,rotMat]=CIRS2TIRS(x,Jul1,Jul2,deltaTTUT1,LOD);
 *
 *Different celestial coordinate systems are compared in [1].
 *
 *REFERENCES:
 *[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double TT1,TT2,*xVec;
    double deltaT=0;
    double LOD=0;
    size_t numRow,numVec;
    double CIRS2TIRS[3][3];
    double Omega[3];//The rotation vector in the TIRS
    mxArray *retMat;
    double *retData;

    if(nrhs<3||nrhs>5){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    checkRealDoubleArray(prhs[0]);
    
    numRow = mxGetM(prhs[0]);
    numVec = mxGetN(prhs[0]);
    
    if(!(numRow==3||numRow==6)) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }
    
    xVec=mxGetDoubles(prhs[0]);
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);
        
    //If some values from the function getEOP will be needed
    if(nrhs<=4||mxIsEmpty(prhs[3])||mxIsEmpty(prhs[4])) {
        mxArray *retVals[5];
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
        
        //We do not need the polar motion coordinates.
        mxDestroyArray(retVals[0]);
        //We do not need the celestial pole offsets.
        mxDestroyArray(retVals[1]);

        //This is TT-UT1
        deltaT=getDoubleFromMatlab(retVals[3]);
        LOD=getDoubleFromMatlab(retVals[4]);
        //Free the returned arrays.
        mxDestroyArray(retVals[2]);
        mxDestroyArray(retVals[3]);
        mxDestroyArray(retVals[4]);
    }

    //If deltaT=TT-UT1 is given
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        deltaT=getDoubleFromMatlab(prhs[3]);
    }
    //If LOD is given
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        LOD=getDoubleFromMatlab(prhs[4]);
    }
    
    //Compute the rotation matrix for going from CIRS to TIRS as well as
    //the instantaneous vector angular momentum due to the Earth's rotation
    //in GCRS coordinates.
    {
        double UT11, UT12;
        double era, omega;
        //Obtain UT1 from terestrial time and deltaT.
        iauTtut1(TT1, TT2, deltaT, &UT11, &UT12);
 
        //Find the Earth rotation angle for the given UT1 time. 
        era = iauEra00(UT11, UT12);
        
        //Construct the rotation matrix.
        CIRS2TIRS[0][0]=1;
        CIRS2TIRS[0][1]=0;
        CIRS2TIRS[0][2]=0;
        CIRS2TIRS[1][0]=0;
        CIRS2TIRS[1][1]=1;
        CIRS2TIRS[1][2]=0;
        CIRS2TIRS[2][0]=0;
        CIRS2TIRS[2][1]=0;
        CIRS2TIRS[2][2]=1;     
        iauRz(era, CIRS2TIRS);
        
        //Next, to be able to transform the velocity, the rotation of the Earth
        //has to be taken into account. 

        //The angular velocity vector of the Earth in the TIRS in radians.
        omega=getScalarMatlabClassConst("Constants","IERSMeanEarthRotationRate");
        //Adjust for LOD
        omega=omega*(1-LOD/86400.0);//86400.0 is the number of seconds in a TT
                                    //day.
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
            iauRxp(CIRS2TIRS, xVec+numRow*curVec, retData+numRow*curVec);
            
            //If a velocity vector was given.
            if(numRow>3) {
                double *posCIRS=xVec+numRow*curVec;
                double posTIRS[3];
                double *velCIRS=xVec+numRow*curVec+3;//Velocity in GCRS
                double velTIRS[3];
                double *retDataVel=retData+numRow*curVec+3;
                double rotVel[3];
                //If a velocity was provided with the position, first
                //convert to TIRS coordinates, then account for the
                //rotation of the Earth.
                
                //Convert velocity from CIRS to TIRS.
                iauRxp(CIRS2TIRS, velCIRS, velTIRS);
                //Convert position from CIRS to TIRS
                iauRxp(CIRS2TIRS, posCIRS, posTIRS);
                                
                //Evaluate the cross product for the angular velocity due
                //to the Earth's rotation.
                iauPxp(Omega, posTIRS, rotVel);
                
                //Subtract out the instantaneous velocity due to rotation.
                iauPmp(velTIRS, rotVel, retDataVel);
            }
        }
    }
    plhs[0]=retMat;
    
    //If the rotation matrix is desired on the output.
    if(nlhs>1) {
        double *elPtr;
        size_t i,j;
        
        plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
        elPtr=mxGetDoubles(plhs[1]);
        
        for (i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                elPtr[i+3*j]=CIRS2TIRS[i][j];
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
