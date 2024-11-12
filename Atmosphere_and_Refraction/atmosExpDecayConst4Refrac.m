function [ce,deltaN]=atmosExpDecayConst4Refrac(Ns)
%%ATMOSEXPDECAYCONST4REFRAC Given the refractivity of air at sea level,
%   obtain the approximate decay constant for the refractivity as a
%   function of height as given in Appendix A of [1]. Also obtain the
%   change in refractivity 1km above sea level. This can be used as
%   the ce input in stdRefracBiasApprox. The refractivity is taken to
%   change as a function of height above sea level (a reference sphere) as
%   N=Ns*exp(-ce*h).
%
%INPUTS: Ns A matrix of refractivity of air at sea level.
%
%OUTPUTS: ce The exponential decay constant of the refractivity with units
%            of inverse meters. This has the same dimensions as Ns.
%     deltaN The change in atmospheric refractivity going 1km up  from sea
%            level. This has the same dimensions as Ns.
%
%This function implements the approximation of Appendix A of [1].
%
%REFERENCES:
%[1] B. R. Bean and G. D. Thayer, CRPL Exponential Reference Atmosphere.
%    Washington, D.C.: U. S. Department of Commerce, National Bureau of
%    Standards, Oct. 1959. [Online]. Available:
%    https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nbsmonograph4.pdf
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

deltaN=-7.32*exp(0.005577*Ns);
ce=log(Ns./(Ns+deltaN))/1e3;

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
