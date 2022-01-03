function latLonEnd=directRhumbSpherProblem(latLonStart,azimuth,dist,r)
%%DIRECTRHUMBSPHERPROBLEM Given a starting location on the surface of a
%               reference sphere as well as a constant heading a distance
%               to travel, determine the location of a vessel traveling at
%               the constant heading after traveling the distance. A
%               constant heading trajectory is a rhumb line (loxodrome) and
%               is usually not the shortest path between two points. This
%               function uses a spherical Earth approximation; the functon
%               directRhumbProblem uses an ellipsoidal Earth
%               approximation.
%
%INPUTS: latLonStart The 2X1 initial point given in latitude and longitude
%                    in radians in the format [latitude;longitude]
%                    (on a reference sphere, latitude is spherical
%                    elevation; longitude is spherical azimuth). This
%                    cannot be a pole.
%            azimuth The constant heading that is to be traveled in
%                    radians East of North.
%               dist The NX1 or 1XN set of distances on the sphere-
%                    approximated Earth starting from latLonStart where one
%                    wishes to find the stopping point. The distances must
%                    be less than the rhumb distance to reach the pole or
%                    invalid results will be obtained.
%                  r The assumed radius of the spherical Earth model. If
%                    omitted or an empty matrix is passed, the default of
%                    osculatingSpher4LatLon(latLonStart) is used.
%
%OUTPUTS: latLonPoints A2XN matrix of geocentric latitude and longitudes of
%                      the final points of the spherical geodesic
%                      trajectory given in radians as [latitude;longitude].
%
%This function implements Equation 8 and 17 in [1]. Equation 17 is used for
%the case where the azimuth is very close to pi/2;
%
%EXAMPLE:
%We show that after solving the indirect rhumb problem, the trajectory
%obtained leads to the desired endpoint when given to
%directRhumbSpherProblem when using the same Earth radius for both.
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% r=osculatingSpher4LatLon(latLonStart);
% [azStart,dist]=indirectRhumbSpherProblem(latLonStart,latLonEnd,r);
% latLonPoint=directRhumbSpherProblem(latLonStart,azStart,dist,r);
% max(abs(wrapRange(latLonPoint-latLonEnd,-pi,pi)))
%The maximum difference is within expected finite precision limits.
%
%REFERENCES:
%[1] G. H. Kaplan, "Practical Sailing Formulas for Rhumb-Line Tracks on an
%    Oblate Earth," Navigation: Journal of the Institute of Navigation,
%    vol. 42, no. 2, pp. 313-326 Summer 1995.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(r))
    r=osculatingSpher4LatLon(latLonStart);
end
dist=dist(:).';

phi0=latLonStart(1);
lambda0=latLonStart(2);

%If the azimuth is too close to due East/West, then use the solution for
%due East/West to avoid numerical precision issues. Otherwise, use the
%solution for other headings in general.
if(abs(abs(azimuth)-pi/2)>1e-8)
    %Equation 8 in [1].
    phi=phi0+(dist/r)*cos(azimuth);
    lambda=lambda0+tan(azimuth)*(log(tan(pi/4+phi/2))-log(tan(pi/4+phi0/2)));
else
    if(azimuth==pi/2)
        phi=phi0*ones(size(dist));
    else
        phi=phi0+(dist/r)*cos(azimuth);
    end
    %Equation 17 in [1] for nearly due East/West, adapted to the sphere. It
    %differs from Equation 9 just in the fact that phi0 and phi are
    %averaged.
    lambda=lambda0+dist*sin(azimuth)./(r*cos((1/2)*(phi0+phi)));
end

latLonEnd=[phi;lambda];
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
