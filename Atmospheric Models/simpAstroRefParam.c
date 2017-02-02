/**SIMPASTROREFPARAM Obtain the parameters A and B for a very simple
*                    astronomical refraction model of the form:
*                    DeltaZ=A*tan(zetaObs)+B*tan(zetaObs)^3, where zetasObs
*                    is the observed zenith angle of a celestial body and
*                    DeltaZ is the bias due to refraction that should be
*                    subtracted from zetaObs to get the true direction of
*                    arrival. Such a model is invalid for observed zenith
*                    angles near 90 degrees. However, the function
*                    simpAstroRefrac, which uses these constants, uses a
*                    modified formula that is more robust.
*
*INPUTS:R  The relative humidity at the observer (between 0 and 1). If
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
*
*OUTPUTS: A The A value for the model in radians.
*         B The B value for the model in radians.
*
*This function is just a wrapper for the function iauRefco in the
*International Astronomical Union's (IAU) Standard's of Fundamental
*Astronomy library.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[A,B]=simpAstroRefParam(R,P,T,wl);
*or
*[A,B]=simpAstroRefParam();
*
*April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double R, P, T, wl;
    double A,B;
    
    if(nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
        
    //Get the relative humidity
    if(nrhs>0&&mxGetM(prhs[0])!=0) {
       R=getDoubleFromMatlab(prhs[0]);
    } else {//Otherwise get the value from the Constants class.
       R=getScalarMatlabClassConst("Constants","standardRelHumid");
    }
    
    //Get the pressure
    if(nrhs>1&&mxGetM(prhs[1])!=0) {
        P=getDoubleFromMatlab(prhs[1]);
    } else {//Otherwise, get the value from the Constants class.
        P=getScalarMatlabClassConst("Constants","standardAtmosphericPressure");
    }
    
    //Get the temperature
    if(nrhs>2&&mxGetM(prhs[2])!=0) {
        T=getDoubleFromMatlab(prhs[2]);
    } else {//Otherwise get the value from the Constants class.
        T=getScalarMatlabClassConst("Constants","standardTemp");
    }

    //wl is inputted in meters, but needs to be in micrometers for the
    //function iauAtco13
    if(nrhs>3&&mxGetM(prhs[3])!=0) {
       wl= getDoubleFromMatlab(prhs[3]);
    } else {
       wl=0.574e-6;
    }
    
    //Convert from Pascals to millibars.
    P*=0.01;
    //Convert from degrees Kelvin to degrees Centigrade.
    T+=-273.15;
    //Convert from meters to micrometers.
    wl*=1e6;
    
    //Get the parameters for the simple refraction model.
    iauRefco(P, T, R, wl,&A, &B);
    
    //Allocate space for the return values.
    plhs[0]=doubleMat2Matlab(&A,1,1);
    if(nlhs>1) {
        plhs[1]=doubleMat2Matlab(&B,1,1);
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
