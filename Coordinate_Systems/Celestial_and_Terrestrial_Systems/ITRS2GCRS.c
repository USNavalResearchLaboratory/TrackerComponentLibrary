/**ITRS2GCRS Convert vectors of position and possibly velocity from the
 *           International Terrestrial Reference System  (ITRS), a type of
 *           Earth-Centered Earth-Fixed (ECEF)  coordinate system, to the
 *           Geocentric Celestrial Reference System (GCRS), a type of
 *           Earth-Centered Inertial (ECI) coordinate system. Note
 *           that the velocity correction includes the centrifugal effects
 *           of the conversion from the ITRS into the terrestrial 
 *           intermediate reference system (TIRS), but omits the effects
 *           of the conversion from the celestial intermediate reference
 *           system (CIRS) into the GCRS. The period of the Celestial
 *           Intermediate Pole (CIP) motion in the GCRS is on the order
 *           of 14 months and thus is significantly smaller than the
 *           rotation effects of the Earth in the TIRS. The velocity
 *           conversion also does not include the (small) centrifugal
 *           effect of polar motion.
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
 * deltaTTUT1 An optional parameter specifying the difference between TT
 *          and UT1 in seconds. This information can be obtained from
 * http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
 *          or 
 * http://www.usno.navy.mil/USNO/earth-orientation/eo-products
 *          If this parameter is omitted or if an empty matrix is passed,
 *          then the value provided by the function getEOP will be used
 *          instead.
 *     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
 *          including the effects of tides and librations. If this
 *          parameter is omitted or if an empty matrix is passed, the value
 *          from the function getEOP will be used.
 *     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
 *          the IAU 2006/2000A precession/nutation model in radians If this
 *          parameter is omitted or if an empty matrix is passed, the value
 *          from the function getEOP will be used.
 *      LOD The difference between the length of the day using terrestrial
 *          time, international atomic time, or UTC without leap seconds
 *          and the length of the day in UT1. This is an instantaneous
 *          parameter (in seconds) proportional to the rotation rate of the
 *          Earth. This is only needed if more than just position
 *          components are being converted.
 *
 *OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from ITRS
 *             coordinates to GCRS coordinates.
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
 *accounting for the LOD EOP. Consequently, one account for rotation by
 *transforming from the ITRS to the TIRS, adding the cross product of
 *Omega with the position in the TIRS, and then converting to the GCRS.
 *This is a simple Newtonian conversion.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[vec,rotMat]=ITRS2GCRS(x,Jul1,Jul2);
 *or if more parameters are known,
 *[vec,rotMat]=ITRS2GCRS(x,Jul1,Jul2,deltaTTUT1,xpyp,dXdY,LOD);
 *
 *Different celestial coordinate systems are compared in [1].
 *
 *REFERENCES:
 *[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 *March 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"
//For sqrt
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numRow,numVec;
    mxArray *retMat;
    double *xVec, *retData;
    double TT1, TT2, UT11, UT12;
    //The if-statements below should properly initialize all of the EOP.
    //The following initializations to zero are to suppress warnings when
    //compiling with -Wconditional-uninitialized.
    double dX=0;
    double dY=0;
    double xp=0;
    double yp=0;
    double deltaT=0;
    double LOD=0;
    double ITRS2GCRS[3][3];
    double TIRS2GCRS[3][3];
    double invrPom[3][3];//Inverse polar motion matrix. TIRS=IPOM*ITRS.
    double Omega[3];

    if(nrhs<3||nrhs>7){
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
    
    //If some values from the function getEOP will be needed.
    if(nrhs<=6||mxIsEmpty(prhs[3])||mxIsEmpty(prhs[4])||mxIsEmpty(prhs[5])||mxIsEmpty(prhs[6])) {
        mxArray *retVals[5];
        double *xpyp, *dXdY;
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
        dXdY=mxGetDoubles(retVals[1]);
        xp=xpyp[0];
        yp=xpyp[1];
        dX=dXdY[0];
        dY=dXdY[1];
        
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
    
    //If no values for the polar motion coordinates are given, then use
    //zeros.
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
    
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        size_t dim1, dim2;
        
        checkRealDoubleArray(prhs[5]);
        dim1=mxGetM(prhs[5]);
        dim2=mxGetN(prhs[5]);
        
        if((dim1==2&&dim2==1)||(dim1==1&&dim2==2)) {
            double *dXdY=mxGetDoubles(prhs[5]);
        
            dX=dXdY[0];
            dY=dXdY[1];
        } else {
            mexErrMsgTxt("The polar motion coordinates have the wrong dimensionality.");
            return;
        }
    }
    
    //If LOD is given
    if(nrhs>6&&!mxIsEmpty(prhs[6])) {
        LOD=getDoubleFromMatlab(prhs[6]);
    }
    
    //Compute the rotation matrix for going from ITRS to GCRS as well as
    //the instantaneous vector angular momentum due to the Earth's rotation
    //in TIRS coordinates.
    {
    double x, y, s, era, sp;
    double rpom[3][3], rc2i[3][3];
    double GCRS2ITRS[3][3];
    double omega;
        
    //Get the X,Y coordinates of the Celestial Intermediate Pole (CIP) and
    //the Celestial Intermediate Origin (CIO) locator s, using the IAU 2006
    //precession and IAU 2000A nutation models.
    iauXys06a(TT1, TT2, &x, &y, &s);
    
    //Add the CIP offsets.
    x += dX;
    y += dY;
    
    //Get the GCRS-to-CIRS matrix
    iauC2ixys(x, y, s, rc2i);
    
    //Find the Earth rotation angle for the given UT1 time.
    era = iauEra00(UT11, UT12);
    
    //Get the Terrestrial Intermediate Origin (TIO) locator s' in radians
    sp=iauSp00(TT1,TT2);
    
    //Get the polar motion matrix
    iauPom00(xp,yp,sp,rpom);
    
    //Combine the GCRS-to-CIRS matrix, the Earth rotation angle, and the
    //polar motion matrix to get a transformation matrix to get the
    //rotation matrix to go from GCRS to ITRS.
    iauC2tcio(rc2i, era, rpom,GCRS2ITRS);
    
    //To go from the ITRS to the GCRS, we need to use the inverse rotation
    //matrix, which is just the transpose of the rotation matrix.
    iauTr(GCRS2ITRS, ITRS2GCRS);
    
    //Next, to be able to transform the velocity, the rotation of the Earth
    //has to be taken into account. This requires first transforming from
    //ITRS to TIRS coordinates, where the rotational axis is the z-axis.
    //That transformation requires the inverse polar motion matrix, which,
    //being a rotation matrix, is given by its transpose.
    iauTr(rpom, invrPom);
    //Then, one must transform from TIRS to GCRS, which can be done by
    //taking the inverse of the rotation matrix from GCRS to ITRS computed
    //using the identity matrix for polar motion (i.e. no polar motion
    //means leaving it in the TIRS.
    {
        double rident[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        double GCRS2TIRS[3][3];
        iauC2tcio(rc2i, era, rident,GCRS2TIRS);
        iauTr(GCRS2TIRS, TIRS2GCRS);
    }
    
    //The angular velocity vector of the Earth in the TIRS in radians.
    omega=getScalarMatlabClassConst("Constants","IERSMeanEarthRotationRate");
    //Adjust for LOD
    omega=omega*(1-LOD/86400.0);//86400.0 is the number of seconds in a TT day.
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
            iauRxp(ITRS2GCRS, xVec+numRow*curVec, retData+numRow*curVec);
            
            //If a velocity vector was given.
            if(numRow>3) {
                double posTIRS[3];
                double velTIRS[3];
                double *posITRS=xVec+numRow*curVec;
                double *velITRS=xVec+numRow*curVec+3;//Velocity in GCRS
                double *retDataVel=retData+numRow*curVec+3;
                double rotVel[3];
                //If a velocity was provided with the position, then first
                //convert into the TIRS, then account for the rotation of
                //the Earth, then rotate into the GCRS.
                
                //Convert velocity from ITRS to TIRS.
                iauRxp(invrPom, velITRS, velTIRS);
                //Convert position from ITRS to TIRS
                iauRxp(invrPom, posITRS, posTIRS);
                
                //Evaluate the cross product for the angular velocity due
                //to the Earth's rotation.
                iauPxp(Omega, posTIRS, rotVel);
                
                //Add the instantaneous velocity due to rotation.
                iauPpp(velTIRS, rotVel, retDataVel);
                
                //Rotate from TIRS to GCRS
                iauRxp(TIRS2GCRS, retDataVel, retDataVel);
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
                elPtr[i+3*j]=ITRS2GCRS[i][j];
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
