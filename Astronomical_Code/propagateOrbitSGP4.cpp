/**PROPAGATEORBITSGP4 Use the SGP4/SDP4 propagator to obtain the position
 *                    and velocity of a non-maneuvering Earth-orbiting
 *                    satellite after a given time offset. The SGP4
 *                    elements are given in the obsolete True Equator Mean
 *                    Equinox (TEME) of date coordinate system as is the
 *                    result. The propagator is of interest,
 *                    because two-line element (TLE) sets, which are sets
 *                    of satellite ephemerides published by the U. S. Air
 *                    Force are given in terms of SGP4 orbital elements.
 *                    The function TLE2SGP4OrbEls  can be used to
 *                    get orbital elements for this function from a TLE.
 *                    The propagator is actually a combination of the
 *                    Simplified General Perturbations 4 (SGP4) dynamic
 *                    model for satellites that remain close to the Earth
 *                    and the Simplified Deep Space Perturbations 4 (SDP4)
 *                    model for satellites that go far enough from the
 *                    Earth that the gravitational effects of the Sun and
 *                    Moon begin to matter. The propagator only works with
 *                    satellites in orbit, not those on an escape
 *                    trajectory and is inaccurate with decaying
 *                    trajectories. The accuracy of the deep space
 *                    propagator, which is used when
 *                    2*pi/SGP4Elements(6)>=225, is not very good.
 * 
 *INPUTS: SGP4Elements A 7X1 or 1X7 vector of SGP4 orbital elements for the
 *                     SGP4 propagator that are close to but not the same
 *                     as Keplerian orbital elements having the same names.
 *                     These are:
 *                     1) eccentricity (0-1)
 *                     2) inclination  (radians)
 *                     3) argument of perigee/periapsis (radians)
 *                     4) right ascension of the ascending node (radians)
 *                     5) mean anomaly (radians)
 *                     6) mean motion (radians per second [TT])
 *                     7) BSTAR drag term. This pseudo ballistic
 *                       coefficient has units of inverse Earth radii.
 *                       Normally, a ballistic coefficient BC is mass per
 *                       unit area divided by a drag coefficient. The BSTAR
 *                       drag term is a scaled version of the normal
 *                       ballistic coefficient. That is,
 *                       BSTAR=BC*rho0/2 where rho0=2.461e-5 kg/m^2. The
 *                       ballistic coefficient BC is Cd*Area/mass where Cd
 *                       is the drag coefficient of the object.
 *              deltaT An NX1 or 1XN vector of the time offsets in seconds
 *                     (TT) from the epoch time at which one wishes to
 *                     determine the state of the target. Note that there
 *                     are exactly 86400 TT seconds in a TT Julian day.
 *   TTEpoch1, TTEpoch2 The epoch time of the SGP4Elements given in
 *                     terrestrial time (TT) as a two-part Julian date (as
 *                     a fractional day count). These are only required if
 *                     2*pi/SGP4Elements(6)>=225, which is when the deep
 *                     space propagator is used. Otherwise, these inputs
 *                     can be omitted or empty matrices can be passed. The
 *                     added precision of using a split date does not
 *                     matter in this instance, since Vallado's SGP4
 *                     library routine just uses the summed value.  
 *             opsMode An optional parameter specifying the orbital
 *                     propagation mode. Possible values are
 *                     0 (The default if omitted) Use a model that is
 *                       supposed to be similar to what the Air Force Space
 *                       Command (AFSPC) uses when they publish orbital
 *                       elements.
 *                     1 Use a model that is supposed to be a slight
 *                       improvement over what the AFSPC uses.
 *        gravityModel An optional parameter specifying which gravitational
 *                     model should be used. Since the gravitational models
 *                     stem from ellipsoidal Earth models, the choice of
 *                     the model also affects the radius of the Earth term
 *                     used in the algorithm. Possible values are
 *                     0 (The default if omitted). Use a gravitational
 *                       model based on the WGS-72 reference ellipsoid.
 *                       This appears to be what the AFSPC uses in their
 *                       TLE sets.
 *                     1 Use a gravitational model based on the WGS-84
 *                       reference ellipsoid.
 *        
 *OUTPUTS:      xState The 6XN target state at deltaT offset from the epoch
 *                     time. The state consists of position and velocity
 *                     components in the order [x;y;z;xDot;yDot;zDot] and
 *                     is in the obsolete TEME coordinate system having
 *                     units of meters (for position) and meters per
 *                     second (for velocity).
 *          errorState This indicates whether any errors occurred during
 *                     the propagation. Possible values are
 *                      0: There is no problem.
 *                      1: Mean element problem. Eccentricity >= 1.0 or
 *                              eccentricity < -0.001 or semimajor axis
 *                              implied by the elements < 0.95.
 *                      2: Mean motion less than 0.
 *                      3: Partially osculating element problem (these are
 *                         derived from the SGP4Elements), eccentricity <
 *                         0.0  or  eccentricity > 1.0
 *                      4: Semi-latus rectum of the orbit implied by the
 *                             elements is < 0.
 *                      5: Epoch elements are sub-orbital. This error is
 *                         no longer used.
 *                      6: Satellite has decayed.
 *
 *The SGP4 propagator is described in [1] and [2]. The implementation here
 *uses the code that Vallado released into the public domain with his
 *paper, and which was downloaded from
 *http://www.centerforspace.com/downloads/
 *as an external library for performing the key propagation step.
 *
 *Note that this is NOT the official SGP4 orbital propagator used by the
 *U.S. Air Force and cannot be assumed to be as reliable or produce
 *identical results to the official propagator. Information on obtaining
 *the U.S. Air Force's official propagator is given at
 *http://www.afspc.af.mil/units/ASDA/
 *
 *The algorithm can be compiled for use in Matlab  using the
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT,TTEpoch1,TTEpoch2,opsMode,gravityModel);
 *or as
 *[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT,TTEpoch1,TTEpoch2);
 *or if no deep-space model is needed, as
 *[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT);
 *
 *REFERENCES:
 *[1] F. R. Hoots and R. L. Roehrich, "Spacetrack report no. 3: Models for
 *    propagation of NORAD element sets," Department of Defense, Tech.
 *    Rep., 31 Dec. 1988. [Online].
 *    Available: http://www.amsat.org/amsat/ftp/docs/spacetrk.pdf
 *[2] D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, "Revisiting
 *    spacetrack report # 3: Rev 2," in Proceedings of the AIAA/AAS
 *    Astrodynamics Specialist Conference and Exhibit, Keystone, CO, 21-24
 *    Aug. 2006. [Online].
 *    Available: http://celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.pdf
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is for the core SGP4 propagation routine.*/
#include "SGP4.h"
/*This header is required by Matlab.*/
#include "mex.h"
/*This is for input validation*/
#include "MexValidation.h"
//For max
#include <algorithm>

#define pi 3.14159265358979323846264338327950288419716939937510582097494459

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //The Julan date in TT at the epoch used in the orbital propagation
    //algorithm: 0:00 January 1 1950
    const double epochDate=2433281.5;
    double *deltaT;
    size_t numTimes;
    double SGP4Elements[7],elementEpochTime;
    bool opsMode=0;
    char opsChar;
    bool gravityModel=0;
    gravconsttype gravConstType;
    mxArray *retMat;
    double *retVec;
    mxArray *errorMat;
    double *errorVals;

    if(nrhs<2||nrhs>6||nrhs==3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    {
        const size_t M=mxGetM(prhs[0]);
        const size_t N=mxGetN(prhs[0]);
        if(!((M==7&&N==1)||(N==7&&M==1))) {
            mexErrMsgTxt("The SGP4 elements have the wrong dimensionality.");
            return; 
        }
    }
    
    checkRealDoubleArray(prhs[0]);
    {
        size_t i;
        double *theEls;
        theEls=mxGetDoubles(prhs[0]);
        
        //Copy the elements into SGP4Elements, adjusting the units from SI
        //to the values used here.
        for(i=0;i<5;i++) {
            SGP4Elements[i]=theEls[i];
        }
        //The multipliation by 60 changes the value from radians per second
        //to radians per minute as desired by the SGP4 code.
        SGP4Elements[5]=theEls[5]*60.0;
        SGP4Elements[6]=theEls[6];
    }
    
    {
        const size_t M=mxGetM(prhs[1]);
        const size_t N=mxGetN(prhs[1]);
        
        if((M!=1&&N!=1)||mxIsEmpty(prhs[1])) {
            mexErrMsgTxt("The times have the wrong dimensionality.");
            return;
        }
        numTimes=std::max(M,N);
    }
    checkRealDoubleArray(prhs[1]);
    deltaT=mxGetDoubles(prhs[1]);

    if(nrhs>2&&!mxIsEmpty(prhs[2])&&!mxIsEmpty(prhs[3])) {
        const double TT1=getDoubleFromMatlab(prhs[2]);
        const double TT2=getDoubleFromMatlab(prhs[3]);
        
        if(TT1>TT2) {
            elementEpochTime=(TT1-epochDate)+TT2;
        } else {
            elementEpochTime=(TT2-epochDate)+TT1;
        }
    } else {
        //No time is given. Make sure that the deep space propagator is not
        //going to be used. Otherwise, the time value does not matter.
        elementEpochTime=0;
        
        if(2*pi/SGP4Elements[5]>=225) {
            mexErrMsgTxt("The elements imply the use of the deep space propagator, but no time was given.");
            return;  
        }
    }
    
    if(nrhs>4) {
        opsMode=getBoolFromMatlab(prhs[4]);
    }
    
    if(opsMode==0) {
        opsChar='a';
    } else {
        opsChar='i';
    }
    
    if(nrhs>5) {
        gravityModel=getBoolFromMatlab(prhs[5]);
    }
    
    if(gravityModel==0) {
        gravConstType=wgs72;
    } else {
        gravConstType=wgs84;
    }
    
    //Allocate space for the return values
    retMat=mxCreateDoubleMatrix(6,numTimes,mxREAL);
    retVec=mxGetDoubles(retMat);
    errorMat=mxCreateDoubleMatrix(numTimes,1,mxREAL);
    errorVals=mxGetDoubles(errorMat);
    {
        //This will be filled with the ephemeris information needed for
        //propagating the satellite.
        elsetrec satRec;
        size_t i;

        SGP4Funcs::sgp4init(gravConstType,//Specified WGS-84 or WGS-72
                 opsChar,//Mode of operation, 'a' or 'i'.
                 0,//Satellite number; unused in the sgp4 function.
                 elementEpochTime,//Epoch time of the orbital elements.
                 0,//xndot; unused in the sgp4 function.
                 0,//xnddot; unused in the sgp4 function.
                 SGP4Elements[6],//BSTAR drag term.
                 SGP4Elements[0],//Eccentricity
                 SGP4Elements[2],//Argument of perigee
                 SGP4Elements[1],//Inclination
                 SGP4Elements[4],//Mean anomaly
                 SGP4Elements[5],//Mean motion
                 SGP4Elements[3],//Righ ascension of the ascending node.
                 satRec);//Gets filled by the function.
                
        //Do the actual propagation.
        //The division by 60 transforms the time value from seconds to
        //minutes, as the SGP4 code wants the result in minutes.
        for(i=0;i<numTimes;i++) {
            double *r=retVec+6*i;
            double *v=retVec+6*i+3;
            int j;
            
            SGP4Funcs::sgp4(satRec,  deltaT[i]/60.0, r,  v);
            
            //Convert the units from kilometers and kilometers per second
            //to meters and meters per second.
            for(j=0;j<3;j++) {
                r[j]=1000*r[j];
                v[j]=1000*v[j];
            }
            
            errorVals[i]=(double)satRec.error;
        }
    }
    
    //Save the result
    plhs[0]=retMat;
    
    //If the error value is desired.
    if(nlhs>1) {
        plhs[1]=errorMat;
    } else{
        mxDestroyArray(errorMat);
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
