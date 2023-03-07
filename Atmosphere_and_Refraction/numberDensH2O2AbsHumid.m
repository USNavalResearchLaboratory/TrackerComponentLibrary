function absHumid=numberDensH2O2AbsHumid(numberDensity)
%%NUMBERDENSH2O2ABSHUMID Given number of water molecules per cubic meter of
%                   air, determine the absolute humidity of the air.
%
%INPUT: numberDensity The number of water moleculer per cubic meter of
%                     air.
%
%OUTPUTS: absHumid The absolute humidity with units of kilograms of
%                  water per cubic meter.
%
%The number density has units of particles per cubic meter. Multiplied by
%the atomic mass of the gas, we have atomic mass units (AMU) per cubic
%meter. Multiplied by the value of the atomic mass unit in kilograms, we
%get kilograms per cubic meter.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The atomic weight of water.
AMUH20=Constants.gasProp('H2O');

absHumid=numberDensity*Constants.atomicMassUnit*AMUH20;

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
