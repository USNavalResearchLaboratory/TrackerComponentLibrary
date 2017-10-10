function specHumid=absHumid2SpecHumid(absHumid,mpVDryAir,defChoice)
%%ABSHUMID2SPECHUMID Convert from absolute humidity (the mass of water in
%                    air per unit volume) to specific humidity for either
%                    definition of specific humidity, assuming the validity
%                    of the Ideal Gas Law. Specific humidity can be defined
%                    as either the ratio of the density(mass per unit
%                    volume) of water in the air to the mass density of the
%                    other gasses, or it can be defined as the mass density
%                    of the water in the air to the total mass density of
%                    the air including the water.
%
%INPUTS: absHumid The absolute humidity with SI units of kilograms of water
%                 per cubic meter.
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
%
%OUTPUTS: specHumid The specific humidity using the given definition. The
%                   specific humidity is unitless.
%
%The specific humidity is just the absolute humidity divided by the
%mass density of dry air under the first definition and it is the absolute
%humidity divided by the sum of the absolute humidity and the mass density
%of dry air under the second definition.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    defChoice=0;
end

if(defChoice~=0)
    specHumid=absHumid/(absHumid+mpVDryAir);
else
    specHumid=absHumid/mpVDryAir;
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
