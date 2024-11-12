function N=approxRefractivity(T,P,Pw)
%%APPROXREFRACTIVITY Given the temperature, pressure and water vapor
%   partial pressure, approximate the refractivity of air using the
%   standard approximation from the International Telecommunications Union
%   (ITU). Compare this function to atmosParam4GasTemp.
%
%INPUTS: T The temperature in degrees Kelvin.
%        P The atmospheric pressure in Pascals (Total pressure: dry
%          pressure plus the partial pressure of water vapor).
%       Pw The partial pressure of water vapor in the atmosphere in
%          Pascals. 
%
%OUTPUTS: N The refractivity of the atmosphere.
%
%This function implements the approximation given in Equation 6 in Annex I
%of [1]. It does not depend on frequency.
%
%REFERENCES:
%[1] International Telecommunication Union, "Recommendation ITU-R
%    P.453-11: The radio refractive index: Its formula and refractivity
%    data," Tech. Rep., Jul. 2015.
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Convert from Pascals to hectopascals.
P=P/100;
Pw=Pw/100;

%Equation 6 in Annex I.
N=77.6*(P/T)-5.6*(Pw/T)+3.75e5*(Pw/T^2);

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
