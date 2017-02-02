/**STARCAT2OBS  Convert data for the location of stars as typically
 *              supplied by a star catalog, such as the Hipparcos catalog,
 *              to local ENU observed coordinates at the receiver,
 *              including corrections for parallax, aberation,
 *              gravitational deflection by the sun and a low-fidelity
 *              refraction model. The catalog data is assumed to be at the
 *              J2000 epoch.
 *
 *catData   catData is a matrix of stars (one per row) that are to be
 *          converted, where each row has the following format:
 *catData(:,1) RArad    Right Ascension in (rad) ICRS at the J2000.0 epoch 
 *catData(:,2) DErad    Declination in (rad) ICRS at the J2000.0 epoch
 *catData(:,3) Plx      Parallax (rad)
 *catData(:,4) pmRA     Proper motion in Right Ascension (rad/yr) in the
 *                      form dRA/dt and not cos(Dec)*dRA/dt.
 *catData(:,5) pmDE     Proper motion in Declination (rad/yr)
 *catData(:,6) vRad     Radial velocity of the star in meters per second
 *                      with a positive sign meaning that a star is
 *                      receding (going away from) Earth.
 *Jul1,Jul2 Two parts of a pseudo-Julian date given in UTC for when the
 *          observation is made. The units of the date are days. The full
 *          date is the sum of both terms. The date is broken into two
 *          parts to provide more bits of precision. It does not matter how
 *          the date is split.
 *zObs      zObs=[lat;lon;h], the longitude, geodetic latitude and height
 *          above the reference ellipsoid of the observer using the WGS-84
 *          reference ellipsoid. The units of lat and lon are radians and
 *          the height is in meters. East and North are the positive
 *          directions.
 *R         The relative humidity at the observer (between 0 and 1). If
 *          this parameter is omitted or an empty matrix is passed, then
 *          Constants.standardRelHumid is used.
 *P         The atmospheric pressure at the observer in Pascals (N/m^2). If
 *          this parameter is omitted or an empty matrix is passed, then
 *          Constants.standardAtmosphericPressure is used.
 *T         The air temperature at the observer in degrees Kelvin. If this
 *          parameter is omitted or an empty matrix is passed, then
 *          Constants.standardTemp is used.
 *wl        The wavelength at which the observation is made in units of
 *          meters. If this parameter is omitted or an empty matrix is
 *          passed, then a wavelength of 0.574 micrometers is used, which
 *          is in the visible spectrum (a rather yellow color).
 *deltaT    The difference between UTC and UT1 in seconds. This
 *          information can be obtained from
 *          http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
 *          or 
 *          http://www.usno.navy.mil/USNO/earth-orientation/eo-products
 *          If this parameter is omitted or if an empty matrix is passed,
 *          then the value provided by the function getEOP will be used
 *          instead.
 *xpyp      xpyp=[xp,yp], the polar motion coordinates of the respect
 *          to the International Terrestrial Reference System. As
 *          described in Section 5.1 of the IERS Conventions 2010,
 *          values are published by the IERS and should have been
 *          updated to account for the additional temporal effects of
 *          ocean tides and librations. If this parameter is omitted,
 *          the value provided by the function getEOP will be used instead.
 *
 *OUTPUTS: zSpher For N stars, this is a 2XN matrix with each colum being
 *                the location of a star in [azimuth;elevation] in radians
 *                taken with respect to a local East-North-Up coordinate
 *                system defined on the WGS-84 ellipsoid. Azimuth is
 *                measured in radians North of East.
 *         uObs   For N stars, this is a 3XN matrix of unit vectors in
 *                WGS-84 ENU coordinates pointing toward the stars. 
 *
 *This is a mex wrapper for the function iauAtco13 in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[zObs,uObs]=starCat2Obs(catData,Jul1,Jul2,zObs,P,T,R,wl,dut1,xpyp);
 *or
 *[zObs,uObs]=starCat2Obs(catData,Jul1,Jul2,zObs);
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"
#include "CoordFuncs.hpp"

const double halfPi=1.5707963267948966192313216916398;
const double pi=3.1415926535897932384626433832795;
//The constant to convert arcseconds to radians. The multiplications go
//as->arcminutes->deg->rad
const double as2Rad=(1.0/60.0)*(1.0/60.0)*(pi/180.0);
const double rad2as=1.0/as2Rad;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *catData;
    mxArray *zSpherMATLAB;
    //This is only used if nlhs>1. It is initialized here to  avoid a
    //warning if compiled with -Wconditional-uninitialized
    mxArray *uObsMATLAB=NULL;
    size_t numStars,i;
    double Jul1,Jul2;
    double uENU[9];//To hold WGS-84 East-North-Up axes at the observer.
    //The if-statements below should properly initialize all of the EOP.
    //The following initializations to zero are to suppress warnings when
    //compiling with -Wconditional-uninitialized.
    double *zObs, P,T,R,wl;
    double xp=0;
    double yp=0;
    double deltaT=0;
    //Pointers to the rows of catData
    double *RArad, *DErad, *Plx, *pmRA, *pmDE,*vRad;

    if(nrhs<4||nrhs>20) {
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    //Get the catalog
    numStars=mxGetM(prhs[0]);
    if(numStars==0||mxGetN(prhs[0])!=6) {
        mexErrMsgTxt("The catalog data has the wrong dimensionality.");
        return;
    }
    checkRealDoubleArray(prhs[0]);
    catData=(double*)mxGetData(prhs[0]);
    
    //Get the UTC time
    Jul1=getDoubleFromMatlab(prhs[1]);
    Jul2=getDoubleFromMatlab(prhs[2]);
    
    //Get the location
    if(mxGetM(prhs[3])!=3||mxGetN(prhs[3])!=1) {
        mexErrMsgTxt("The observer location has the wrong dimensionality.");
        return;
    }
    checkRealDoubleArray(prhs[3]);
    zObs=(double*)mxGetData(prhs[3]);
    
    //Get the relative humidity
    if(nrhs>4&&mxGetM(prhs[4])!=0) {
       R=getDoubleFromMatlab(prhs[4]);
    } else {//Otherwise get the value from the Constants class.
       R=getScalarMatlabClassConst("Constants","standardRelHumid");
    }
    
    //Get the pressure
    if(nrhs>5&&mxGetM(prhs[5])!=0) {
        P=getDoubleFromMatlab(prhs[5]);
    } else {//Otherwise, get the value from the Constants class.
        P=getScalarMatlabClassConst("Constants","standardAtmosphericPressure");
    }
    
    //Get the temperature
    if(nrhs>6&&mxGetM(prhs[6])!=0) {
        T=getDoubleFromMatlab(prhs[6]);
    } else {//Otherwise get the value from the Constants class.
        T=getScalarMatlabClassConst("Constants","standardTemp");
    }
    
    if(nrhs>7&&mxGetM(prhs[7])!=0) {
       wl= getDoubleFromMatlab(prhs[7]);
    } else {
       wl=0.574e-6;
    }
    
    //Convert the meteorlogical units for use in the iauAtco13 function.
    //Convert from Pascals to millibars.
    P*=0.01;
    //Convert from degrees Kelvin to degrees Centigrade.
    T-=273.15;
    //Convert from meters to micrometers.
    wl*=1e6;
    
    //If any default values will be needed, load them.
    if(nrhs<=9||mxGetM(prhs[8])==0||mxGetM(prhs[9])==0){
        mxArray *retVals[3];
        double *xpyp,*dXdY;
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        int retVal;
        
        //Get the time in UTC to look up the parameters by going to TAI and
        //then UTC.
        retVal=iauTttai(Jul1, Jul2, &JulUTC[0], &JulUTC[1]);
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
        mexCallMATLAB(3,retVals,2,JulUTCMATLAB,"getEOP");
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
        
        checkRealDoubleArray(retVals[0]);
        checkRealDoubleArray(retVals[1]);
        if(mxGetM(retVals[0])!=2||mxGetN(retVals[0])!=1||mxGetM(retVals[1])!=2||mxGetN(retVals[1])!=1) {
            mxDestroyArray(retVals[0]);
            mxDestroyArray(retVals[1]);
            mxDestroyArray(retVals[2]);
            mexErrMsgTxt("Error using the getEOP function.");
            return;
        }
        xpyp=(double*)mxGetData(retVals[0]);
        dXdY=(double*)mxGetData(retVals[1]);//The celestial pole offsets are not used.
        xp=xpyp[0];
        yp=xpyp[1];
        deltaT=getDoubleFromMatlab(retVals[2]);
        
        //Free the returned arrays.
        mxDestroyArray(retVals[0]);
        mxDestroyArray(retVals[1]);
        mxDestroyArray(retVals[2]);
    }

    //Get the UTC UT1 offset, if provided.
    if(nrhs>8&&mxGetM(prhs[8])!=0) {
        deltaT=getDoubleFromMatlab(prhs[8]);
    }
    
    //Get the components of the polar motion coordinates.
    if(nrhs>9&&mxGetM(prhs[9])!=0) {
        double *XpYp;
        if(mxGetM(prhs[9])!=2&&mxGetN(prhs[9])!=1) {
            mexErrMsgTxt("The polar motion coordinate vector has the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[9]);
        XpYp=(double*)mxGetData(prhs[9]);
        
        xp=XpYp[0];
        yp=XpYp[1];
    }
    
    //Get the axes at the local observer.
    {
        double c[3];
        double a=getScalarMatlabClassConst("Constants","WGS84SemiMajorAxis");
        double f=getScalarMatlabClassConst("Constants","WGS84Flattening");

        getENUAxesCPP(uENU, c,zObs,false,a,f);
    }

    //Set the pointers for each column of data in the catalog.
    RArad=catData;
    DErad=RArad+numStars;
    Plx=DErad+numStars;
    pmRA=Plx+numStars;
    pmDE=pmRA+numStars;
    vRad=pmDE+numStars;
    
    //Allocate space for the return values
    zSpherMATLAB=mxCreateDoubleMatrix(2,numStars,mxREAL);
    plhs[0]=zSpherMATLAB;
    
    if(nlhs>1) {
        uObsMATLAB=mxCreateDoubleMatrix(3,numStars,mxREAL);
        plhs[1]=uObsMATLAB;
    }
    
    {
        double *uEast,*uNorth,*uUp;
        double *zRet;
        //This is only used if nlhs>1. It is initialized here to  avoid a
        //warning if compiled with -Wconditional-uninitialized
        double *uRet=NULL;
        uEast=uENU;
        uNorth=uENU+3;
        uUp=uENU+6;
        
        zRet=(double*)mxGetData(zSpherMATLAB);
        if(nlhs>1) {
            uRet=(double*)mxGetData(uObsMATLAB);
        }
        
        for(i=0;i<numStars;i++) {
            //Variables to hold the result of the iauAtco13 function.
            double aob,zob,hob,dob,rob,eo;
            int retVal;
            
            retVal=iauAtco13(RArad[i],//ICRS right ascension at J2000.0, radians.
                          DErad[i],//ICRS declination at J2000.0, radians.
                          pmRA[i],//RA proper motion, radians/year (in the form dRA/dt and not cos(Dec)*dRA/dt).
                          pmDE[i],//Dec proper motion (radians/year).
                          Plx[i]*rad2as,//parallax (arcseconds).
                          vRad[i]/1000.0,//radial velocity (km/s, +ve if receding).
                          Jul1, Jul2,//Quasi-Julian UTC date.
                          -deltaT,//UT1-UTC in seconds.
                          zObs[1],//WGS-84Longitude, radians East.
                          zObs[0],//WGS-84 geodetic latitude (radians North).
                          zObs[2],//WGS-84 Ellipsoidal height in meters.
                          xp,yp,//polar motion coordinates (radians)
                          P,//Pressure at the observer in millibars (hectoPascals).
                          T,//Temperature at the observer (deg C).
                          R,//Relative humidity at the observer (0-1).
                          wl,//Wavelength (micrometers).
                          &aob,//Observed azimuth (radians East of North).
                          &zob,//Observed zenith distance (radians).
                          &hob,//Observed hour angle (radians).
                          &dob,//Observed declination (radians North).
                          &rob,//Observed CIO-based right ascension (radians).
                          &eo);//Equation of the Origins (ERA-GST).

            if(retVal<0) {
                mxDestroyArray(zSpherMATLAB);
                if(nlhs>1) {
                    mxDestroyArray(uObsMATLAB);
                }
                
                mexErrMsgTxt("An error occurred during the transformation to local coordinates.");
                return;
            }

            zRet[0]=halfPi-aob;//Convert radians East of North to North of East
            zRet[1]=halfPi-zob;
            
            if(nlhs>1) {
                //Convert the  azimith and zenith distance (zenith
                //distance=pi/2-elevation) to a unit vector in the local
                //ENU coordinate system of the observer.
                uRet[0]=cos(aob)*sin(zob);
                uRet[1]=sin(aob)*sin(zob);
                uRet[2]=cos(zob);
                
                uRet+=3;
            }
            zRet+=2;
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
