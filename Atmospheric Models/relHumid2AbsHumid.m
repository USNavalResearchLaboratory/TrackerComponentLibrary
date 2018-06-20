function absHumid=relHumid2AbsHumid(relHumid,T,algChoice)
%%RELHUMID2ABSHUMID Convert relative humidity to absolute humidity under
%                   the assumption of the validity of the Ideal Gas Law and
%                   Dalton's Law of Partial Pressures.
%
%INPUTS: relHumid  The relative humidity as a fraction from 0 to 1.
%                T The temperature in degrees Kelvin.
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
%OUTPUTS: absHumid The absolute humidity with SI units of kilograms of 
%                  water per cubic meter.
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
%per mol). Thus, after some algebra, one obtains
%(m/V)=P*M/(R*T)
%where the ratio m/V is the absolute humidity. Dalton's law of partial
%pressues came into play in that we assumed that the ideal gas law could be
%applied to a single constituent of the gas mixture in the atmosphere
%independently of all others.
%
%More information on the algorithms for the dew point calculation is given
%in the comments to the function dewPointPres4Temp.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    algChoice=0;
end

%The partial pressure of water at the dew point for the given temperature.
PSat=dewPointPres4Temp(T,algChoice);

%The partial pressure of water in the air based on the relative humidity.
PH2O=relHumid*PSat;
R=Constants.molarGasConstant;
%Water is H2O, so its molar mass should be 
M=2*Constants.elementAMU(1)+Constants.elementAMU(8);

%The extra factor of 1/1000 is to convert the units to kilograms per cubic
%meter.
absHumid=(1/1000)*PH2O*M/(R*T);

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

