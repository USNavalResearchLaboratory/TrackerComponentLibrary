function relHumid=specHumid2RelHumid(specHumid,T,mpVDryAir,defChoice,algChoice)
%SPECHUMID2RELHUMID  Convert a specific humidity to a relative humidity
%                    assuming the validity of the Ideal Gas Law.
%
%INPUTS: specHumid The specific humidity which is a nonnegative number.
%                T The temperature in degrees Kelvin.
%        mpVDryAir The mass density (mass in kilograms per unit volume in
%                  cubic meters) of dry air (the air not counting the
%                  water).
%        defChoice An optimal parameter specifying the definition of
%                  specific humidity to use. The choices are
%                  0 Define specific humidity as the mass density of water
%                    over the mass density of dry air.
%                  1 Define specific humidity as the mass density of water
%                    over the total mass density of the air.
%                  If defChoice is omitted, then the default value of 0 is
%                  used.
%        algChoice An optional parameter specifying the algorithm used for
%                  the dew point computation. The choices are
%                  0 The corrected version of the Clausius-Clapeyron
%                    equation for use over land or in the upper air.
%                  1 The empirical Magnus-type equation for use over water
%                    or in the upper air (-40C to 50C temperature).
%                  2 The empirical Magnus-type equation for use over ice
%                   (-80C to 0C temperature)
%                  If algChoice is omitted, then the default value of 0,
%                  the corrected Clausius-Clapeyron equation is used.
%
%OUTPUTS: relHumid The relative humidity as a fractional value from 0 to 1.
%
%The function relHumid2SpecHumid describes how the specific humidity is
%related to the relative humidity. This function just invers the procedure.
%
%More information on the algorithms for the dew point calculation is given
%in the comments to the function dewPointPres4Temp.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(defChoice))
    defChoice=0;
end

if(nargin<5||isempty(algChoice))
    algChoice=0;
end

if(defChoice~=0)
    mpVH2O=mpVDryAir*specHumid/(1-specHumid);
else
    mpVH2O=mpVDryAir*specHumid;
end

R=Constants.molarGasConstant;
%Water is H2O, so its molar mass should be 
MH2O=2*Constants.elementAMU(1)+Constants.elementAMU(8);

%The factor of 1000 converts from kilograms per cubic meter to grams per
%cubic meter.
pH2O=1000*mpVH2O*R*T/MH2O;

%The partial pressure of water at the dew point for the given temperature.
PSat=dewPointPres4Temp(T,algChoice);

%The partial pressure of water in the air based on the relative humidity.
relHumid=pH2O/PSat;

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
