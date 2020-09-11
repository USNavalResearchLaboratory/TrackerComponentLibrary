function specHumid=relHumid2SpecHumid(relHumid,T,mpVDryAir,defChoice,algChoice)
%RELHUMID2SPECHUMID  Convert a relative humidity value to a speciic
%                    humidity value, for either definition of specific
%                    humidity, assuming the validity of the Ideal Gas Law.
%                    Specific humidity can be defined as either the ratio
%                    of the density(mass per unit volume) of water in the
%                    air to the mass density of the other gasses, or it can
%                    be defined as the mass density of the water in the air
%                    to the total mass density of the air including the
%                    water.
%
%INPUTS: relHumid The relative humidity as a fraction from 0 to 1.
%               T The temperature in degrees Kelvin.
%       mpVDryAir The mass density (mass in kilograms per unit volume in
%                 cubic meters) of dry air (the air not counting the
%                 water).
%       defChoice An optimal parameter specifying the definition of
%                 specific humidity to use. The choices are
%                 0 Define specific humidity as the mass density of water
%                   over the mass density of dry air.
%                 1 Define specific humidity as the mass density of water
%                   over the total mass density of the air.
%                 If defChoice is omitted, then the default value of 0 is
%                 used.
%       algChoice An optional parameter specifying the algorithm used for
%                 the dew point computation. The choices are
%                 0 The corrected version of the Clausius-Clapeyron
%                   equation for use over land or in the upper air.
%                 1 The empirical Magnus-type equation for use over water
%                   or in the upper air (-40C to 50C temperature).
%                 2 The empirical Magnus-type equation for use over ice
%                  (-80C to 0C temperature)
%                 If algChoice is omitted, then the default value of 0,
%                 the corrected Clausius-Clapeyron equation is used.
%
%OUTPUTS: specHumid The specific humidity using the given definition. The
%                   specific humidity is unitless.
%
%Given the temperature, the partial pressure of water at saturation (the
%dew point) can be found (pSat). The actual partial pressure of water in
%the atmosphere is then just pH2O=relHumid*pSat. Next, the Ideal Gas Law,
%which one can find in most introductory textbooks on physics is used. The
%law says that
%P*V=n*R*T
%where P= pressure, V is volume, n is the number of moles (a count of
%atoms) of the substance in the volume, R is the universal gas constant and
%T is the temperature in degrees Kelvin. n can be written as m/M where m is
%the mass of the substance in the volume and M is the molar mass of the
%substance (often m is in grams and M is in atomic mass units a.k.a. grams
%per mol).
%
%When considering only the water in the air, one can use the Ideal Gas Law
%to write:
%mH2O/V=pH2O*MH2O/(R*T)
%where mH2O is the mass of water in volume V and MH2O is the molar mass of
%the water. Thus, if one knows the mass per unit volume of the remaining
%gasses in the air, one can find the specific humidity. The mass per unit
%volume of dry air and a temperature can be found in a standard atmosphere
%at a particular time and place using
%[N,mpVDryAir,T,P]=standardAtmosParam(Jul1,Jul2,point);
%
%More information on the algorithms for the dew point calculation is given
%in the comments to the function dewPointPres4Temp.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    defChoice=0;
end

if(nargin<5)
    algChoice=0;
end

%The partial pressure of water at the dew point for the given temperature.
PSat=dewPointPres4Temp(T,algChoice);

%The partial pressure of water in the air based on the relative humidity.
pH2O=relHumid*PSat;
R=Constants.molarGasConstant;
%Water is H2O, so its molar mass should be 
MH2O=2*Constants.elementAMU(1)+Constants.elementAMU(8);

%The factor of 1/1000 converts from grams per cubic meter to kilograms per
%cubic meter.
mpVH2O=(1/1000)*pH2O*MH2O/(R*T);

if(defChoice~=0)
    specHumid=mpVH2O/(mpVDryAir+mpVH2O);
else
    specHumid=mpVH2O/mpVDryAir;
end

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
