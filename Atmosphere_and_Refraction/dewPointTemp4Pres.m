function T=dewPointTemp4Pres(p,algChoice)
%%DEWPOINTTEMP4PRES For a given partial vapor pressure of water, find the
%                   dew point temperature. that is the temperature at which
%                   the air would be saturated with water, in equilibrium.
%
%INPUTS: p The (partial) pressure of the water vapor for which the dew
%          point temperature is desired in units of Newtons per square
%          meter (Pascals).
% algChoice An optional parameter specifying the algorithm to use. The
%          choices are
%          0 The corrected version of the Clausius-Clapeyron equation for
%            use over land or in the upper air.
%          1 The empirical Magnus-type equation for use over water or in
%            the upper air (-40C to 50C temperature).
%          2 The empirical Magnus-type equation for use over ice (-80C to
%            0C temperature)
%          If algChoice is omitted, then the default value of 0, the
%          corrected Clausius-Clapeyron equation is used. Choice 0 is taken
%          from [1]. Choices 1 and 2, are taken from [2].
%
%As is the case in the function dewPointPres4Temp, which is the inverse of
%this function, formulae 0 is from [1], where the numerical inversion
%method given in equations 44 and 45 in Section 5 of the paper are used.
%
%Formulas 1 and 2 are obtained by solving for the inverse of the Magnus
%approximations given in [2]. Simple, explicit solutions for the function
%inverses can be found.
%
%Note that the function dewPointPres4Temp is the inverse of this function.
%
%REFERENCES:
%[1] D. Koutsoyiannis, "Clausius-Clapeyron equation and saturation vapour
%    pressure: simple theory reconcided with practice," European Journal of
%    Physics, vol. 33, no. 2, pp. 295-305, Mar. 2012.
%[2] O. A. Alduchov and R. E. Eskridge, "Improved Magnus Form Approximation
%    of Saturation Vapor Pressure," Journal of Applied Meteorology, vol.
%    35, no. 4, pp. 601-609, Apr. 1996.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algChoice))
    algChoice=0;
end

switch(algChoice)
    case 1
        %This evaluates the inverse of Equation 21 in the Alduchov and
        %Eskridge paper.
        
        p=p/100;%Convert to hectoPascals.
        logRat=log(p/6.1094);
        
        T=-243.04*logRat./(logRat-17.625);
        %Convert the temperature to degrees Kelvin.
        T=T-Constants.absoluteZero;
    case 2
        %This evaluates the inverse of Equation 23 in the Alduchov and
        %Eskridge paper for the saturation temperature over ice.
        
        p=p/100;%Convert to hectoPascals.
        logRat=log(p/6.1121);
        
        T=-273.86*logRat./(logRat-22.587);      
        %Convert the temperature to degrees Kelvin.
        T=T-Constants.absoluteZero;
    otherwise
        %The inverse of Equation 23 in the Koutsoyiannis paper is not
        %simple to find. However, the paper provides iterative solutions in
        %the form of Equations 44 and 45 for specific values related to
        %water. For general values, the initial estimate in Equation 44 is
        %T0/T=1+1/(alpha/(R*T0)-(cL-cP)/R)*ln(p0/p);
        %and the iteration to refine the ratio is
        %T0/T_{new}=1+(R*T0/alpha)*ln(p0/p)+((cL-cP)/R)/(alpha/(R*T0))*ln(T0/T_{old})
        %where R is the specific gas constant of water vapor with units of
        %J/(kg*K), T0 is the temperature at the triple point of water, with
        %units of Kelvin, p0 is the pressure at the triple point of water
        %in hecoPascals, cP and cL are the specific heat of water vapor and
        %liquid water at constant pressure with units of J/(kg*K), and
        %alpha=L0+(cL-cP)*T0, where L0 is the latent heat of water for a
        %constant pressure at the triple point with units of J/kg.
        %
        %In the implementation here, the specific numerical from the paper
        %are used. This means that alpha has been tweaked a little bit to
        %better match their real data.

        numIter=27;

        p0=6.11657*100;%The pressure at the triple point of water in Pascals.
        T0=273.16;%The temperature at the triple point of water in Kelvin.
        
        LPRat=log(p0./p);
        
        TRat=1+1/(24.921-5.06)*LPRat;
        for curIter=1:numIter
            TRat=1+(1/(24.921))*LPRat+(5.06/24.921).*log(TRat);
        end
        T=T0./TRat;
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
