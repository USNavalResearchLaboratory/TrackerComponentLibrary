function relHumid=absHumid2RelHumid(absHumid,T,algChoice)
%%ABSHUMID2RELHUMID Convert relative humidity to absolute humidity under
%                   the assumption of the validity of the Ideal Gas Law and
%                   Dalton's Law of Partial Pressures.
%
%INPUTS: absHumid The absolute humidity with SI units of kilograms of water
%                 per cubic meter.
%               T The temperature in degrees Kelvin.
%       algChoice An optional parameter specifying the algorithm used for
%                 the dew point computation. The choices are
%                 0 The corrected version of the Clausius-Clapeyron
%                   equation for use over land or in the upper air.
%                 1 The empirical Magnus-type equation for use over water
%                   or in the upper air (-40C to 50C temperature).
%                 2 The empirical Magnus-type equation for use over ice
%                   (-80C to 0C temperature)
%                 If algChoice is omitted, then the default value of 0,
%                 the corrected Clausius-Clapeyron equation is used.
%
%OUTPUTS: relHumid  The relative humidity as a fraction from 0 to 1.
%
%The function relHumid2AbsHumid describes the relationship between relative
%and absolute humidities. This function simply performs the inverse
%operation.
%
%More information on the algorithms for the dew point calculation is given
%in the comments to the function dewPointPres4Temp.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    algChoice=0;
end

R=Constants.molarGasConstant;
%Water is H2O, so its molar mass should be 
M=2*Constants.elementAMU(1)+Constants.elementAMU(8);

%The extra factor of 1000 is to convert the units to grams per cubic
%meter.
PH2O=1000*absHumid*R*T/M;

%The partial pressure of water at the dew point for the given temperature.
PSat=dewPointPres4Temp(T,algChoice);

%The partial pressure of water in the air based on the relative humidity.
relHumid=PH2O/PSat;
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
