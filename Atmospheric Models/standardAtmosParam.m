function [N,rho,T,P]=standardAtmosParam(dayOfYear,secondOfTheDay,latLonAlt,rhow)
%%STANDARDATMOSPARAM  Get basic parameters for atmospheric refractivity,
%                     density, pressure and temperature from the
%                     NRLMSISE-00 atmospheric model. The model is best
%                     suited for altitudes below 90km as the anomalous
%                     oxygen parameter of the NRLMSISE-00 model is not
%                     used by this function.
%
%INPUTS: dayOfYear The integer day of the year in the Gregorian calendar in
%                  universal coordinated time (UTC). Counting starts at 1.
%                  The resolution of the model is not sufficient for it to
%                  matter whether 365 or 366 is given as the day if it
%                  isn't/is a leap year.
%   secondOfTheDay The second of the day. This starts at zero. The
%                  resolution of the model is not high enough for leap
%                  seconds to matter, so values above 86400.0 are just
%                  clipped to 86400.0.
%        latLonAlt The location under consideration given in WGS-84
%                  ellipsoidal coordinates of latitude and longitude in
%                  radians and ellipsoidal height in meters.
%             rhow An optional parameter specifying the mass density of
%                  water vapor at the point in question in units of
%                  kilograms per cubic meter. If omitted, the air is
%                  assumed to be dry (rhow=0). The total density of the air
%                  is assumed to be the sum of the dry air density and
%                  rhow.
%
%OUTPUTS: N The refractivity of the atmosphere. In this model, N is always
%           real. N=10^6*(n-1) where n is the index of refraction. This is
%           generally valid for frequencies from L-band (1GHz) to 10 GHz
%           (the middle of X-band).
%       rho The atmospheric density at the point in question in units of
%           kilograms per cubic meter. 
%         T The temperature at the point in question with units of degrees
%           Kelvin.
%         P The atmospheric pressure at the point in question in units of
%           Newtons per square meter (Pascals). It assumes that the gasses
%           can be treated as ideal gasses.
%
%The dry air density and temperature is obtained from the NRLMSISE-00
%atmospheric model using the default parameters for magnetic and solar
%parameters, which are generally valid below 90km, using the function
%NRLMSISE00GasTemp. The reslting standard atmospheric parameters are then
%passed to the atmosParam4GasTemp function.
%
%The atmospheric pressure is obtained using the Ideal Gas Law as described
%in [1]. The properties of atmospheric constituents used in
%atmosParam4GasTemp are more recent than those used in the NRLMSISE-00
%standard. Consequently, pressure obtained for the associated altitude is
%not completely numerically consistent with the function
%NRLMSISE00Alt4Pres.
%
%REFERENCES:
%[1] D. P. Drob, THE NRLMSISE-00 AND HWM-93 USERS GUIDE: Version 1.50,
%    Nov. 2003.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
   rhow=0; 
end

[gasTable,t]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt);
T=t(2);%The temperature in Kelving at the point.

[N,rho,P]=atmosParam4GasTemp(gasTable,T,rhow);
end

%LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
