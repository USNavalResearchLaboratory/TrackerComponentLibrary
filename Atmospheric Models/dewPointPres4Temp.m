function p=dewPointPres4Temp(T,algChoice)
%%DEWPOINTPRES4TEMP Obtain the saturation pressure of water for a given
%                   temperature. This is the partial vapor pressure of
%                   water that (in equilibrium) can not be exceeded. In
%                   other words, this is the dew point pressure. Above this
%                   pressure, water would begin to come out of gaseous form
%                   at this temperature to restore equilibrium.
%
%INPUTS: T A vector or matrix of temperatures in degrees Kelvin.
% algChoice An optional parameter specifying the algorithm to use. The
%          choices are
%          0 The corrected version of the Clausius-Clapeyron equation for
%            use over land or in the upper air.
%          1 The empirical Magnus-type equation for use over water or in
%            the upper air (-40C to 50C temperature).
%          2 The empirical Magnus-type equation for use over ice (-80C to
%            0C temperature)
%          If algChoice is omitted, then the default value of 0, the
%          corrected Clausius-Clapeyron equation of [1,2,3] is used.
%          Choices 1 and 2 use the methods from [4].
%
%OUTPUTS: p The dew point pressure in units of Newtons per square meter
%           (Pascals) corresponding to each of the temperatures in T.
%
%Formulae 0 is from [1], where formula 1 is also mentioned.
%Two minor comments were printed on that paper, specifically [2] and [3].
%
%The paper demonstrates how the traditional form of the frequently used
%Clausius-Clapeyron equation for finding the saturation of pressure of air
%is based on the incorrect assumption that the latent heat of vaporization
%(also known as the anthalpy of vaporization) of water is constant as a
%function of temperature. This assumption leads to inconsistent results.
%The paper derives a more accurate form of the Clausius-Clapeyron equation
%by relaxing that assumption. The corrected version of the Clausius-
%Clapeyron equation is shown to be more accurate than the empirical Magnus-
%type equation (algorithm 1). Note that by switching the constants, the
%corrected Clausius-Clapeyron equation can be used to find the saturation
%pressure of other substances, such as acetone.
%
%The Magnus-type equations are from [4]. Note that an additional pressure-
%dependent Magnus type equation to deal with departures from the ideal gas
%law is also presented in that paper, but is not provided here.
%
%Note that the function dewPointTemp4Pres is the inverse of this function.
%
%REFERENCES:
%[1] D. Koutsoyiannis, "Clausius-Clapeyron equation and saturation vapour
%    pressure: simple theory reconcided with practice," European Journal of
%    Physics, vol. 33, no. 2, pp. 295-305, Mar. 2012.
%[2] T. López-Arias, "Comment on 'Clausius-Clapeyron equation and
%    saturation vapour pressure: simple theory reconcided with practice',"
%    European Journal of Physics, vol. 33, no. 3, pp. L11-L12, May 2012.
%[3] D. Koutsoyiannis, "Reply to 'Comment on "Clausius-Clapeyron equation
%    and saturation vapour pressure: simple theory reconcided with
%    practice"'," European Journal of Physics, vol. 33, no. 3, pp. L13-L14,
%    May. 2012.
%[4] O. A. Alduchov and R. E. Eskridge, "Improved Magnus Form Approximation
%    of Saturation Vapor Pressure," Journal of Applied Meteorology, vol.
%    35, no. 4, pp. 601-609, Apr. 1996.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    algChoice=0;
end

switch(algChoice)
    case 1
        %Convert the temperature to degrees Centigrade.
        T=T+Constants.absoluteZero;
        %Equation 21 in the Alduchov and Eskridge paper. The result is in
        %hectoPascals.
        p=6.1094*exp(17.625*T./(243.04+T));
    case 2
        %Convert the temperature to degrees Centigrade.
        T=T+Constants.absoluteZero;
        %Equation 23 in the Alduchov and Eskridge paper for the saturation
        %pressure over ice. The result is in hectoPascals.
        p=6.1121*exp(22.587*T./(273.86+T));
    otherwise
        %Equation 23 in the March 2012 Koutsoyiannis paper with the
        %parameters for water vapor from Section 5. The equation is 
        %p=p0*exp(alpha/(R*T0)*(1-T0/T))*(T0/T)^((cL-cP)/R);
        %where R is the specific gas constant of water vapor with units of
        %J/(kg*K), T0 is the temperature at the triple point of water, with
        %units of Kelvin, p0 is the pressure at the triple point of water
        %in hecoPascals, cP and cL are the specific heat of water vapor and
        %liquid water at constant pressure with units of J/(kg*K), and
        %alpha=L0+(cL-cP)*T0, where L0 is the latent heat of water for a
        %constant pressure at the triple point with units of J/kg.
        %
        %However, to make the data better match reality, the authors
        %tweaked alpha. Thus, the specific values of their constants are
        %used here rather than using the most accurate theoretical
        %formulation based upon available information regarding water at
        %the triple point.
        
        %The Koutsoyiannis paper's value of the pressure of H2O at the
        %triple point in hectoPascals.
        p0=6.11657;
        T0=273.16;%The temperature at the triple point of water in Kelvin.
        p=p0*exp(24.921*(1-T0./T)).*(T0./T).^5.06;
end
    %Convert from hectoPascals to Pascals.
    p=100*p;
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
