function [azimuth,dist]=indirectRhumbSpherProblem(latLonStart,latLonEnd,r)
%%INDIRECTRHUMBSPHERPROBLEM Given a starting and ending latitude and
%               longitude on a reference sphere, determine the heading (in
%               radians East of North) and the distance one must travel on
%               the shortest constant-heading course to go from the
%               starting point to the stopping point. A constant heading
%               course follows a rhumb line (loxodrome) and is usually not
%               the shortest path between two points. This function uses a
%               spherical Earth approximation; the functon
%               indirectRhumbProblem uses an ellipsoidal Earth
%               approximation.
%
%INPUTS: latLonStart The 2X1 initial point given in latitude and longitude
%                    in radians in the format [latitude;longitude]
%                    (on a reference sphere, latitude is spherical
%                    elevation; longitude is spherical azimuth). This
%                    cannot be a pole.
%          latLonEnd The 2X1 final point given in latitude and longitude in
%                    radians in the format [latitude;longitude].
%                  r The assumed radius of the spherical Earth model. If
%                    omitted or an empty matrix is passed, the default of
%                    osculatingSpher4LatLon(latLonStart) is used.
%
%OUTPUTS: azimuth The constant heading in radians East of North that one
%                 must travel to go on a constant-heading course from
%                 latLonStart to latLonEnd. The azimuth value will not be
%                 accurate if latLonStart is at a geographic pole.
%            dist The distance that one must travel on a constant-heading
%                 course to go from latLonStart to latLonEnd.
%
%This function implements Equation 2 and 5 in [1]. For the case where
%abs(azimuth) is very close to pi/2, a limit in Equation 5 is taken so that
%the distance can still be obtained.
%
%A spherical approximation of the Earth is worse than an ellipsoidal
%approximation. However, spherical rhumb lines are faster to compute,
%because explicit expressions are available.
%
%EXAMPLE:
%This computes a trajectory that crosses the international date and and
%goes from the Northern hemisphere to the southern hemisphere. We show that
%solving the direct problem with the computed heading and distance provides
%the original point.  We then plot the trajectory on an image of the
%spherical Earth. For better plotting, the radius of the Earth has been
%normalized to 1.
% N=500;
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% 
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% [azimuth,dist]=indirectRhumbSpherProblem(latLonStart,latLonEnd,1);
% 
% distVals=linspace(0,dist,N);
% latLonWayPoints=directRhumbSpherProblem(latLonStart,azimuth,distVals,1);
% 
% %Show that the direct algorithm reaches the same endpoint as the indirect
% %algorithm.
% max(abs(wrapRange(latLonWayPoints(:,N)-latLonEnd,-pi,pi)))
% 
% xStartCart=ellips2Cart([latLonStart;0],1,0);
% xEndCart=ellips2Cart([latLonEnd;0],1,0);
% %The path is displayed slightly above the Earth's surface to make it
% %easier to see.
% pathPoints=ellips2Cart([latLonWayPoints;0.02*ones(1,N)],1,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],1,0);
% scatter3(xStartCart(1),xStartCart(2),xStartCart(3),100,'filled')
% scatter3(xEndCart(1),xEndCart(2),xEndCart(3),100,'filled')
% plot3(pathPoints(1,:),pathPoints(2,:),pathPoints(3,:),'-r','linewidth',4)
%
%REFERENCES:
%[1] J. Alexander, "Loxodromes: A rhumb way to go," Mathematics Magazine,
%    vol. 77, no. 5, pp. 349-356, Dec. 2004.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(r))
    r=osculatingSpher4LatLon(latLonStart);
end

lat1=latLonStart(1);
lat2=latLonEnd(1);
lon1=latLonStart(2);
lon2=latLonEnd(2);

%Expressions for Sigma are from Equation 6 in [1].
Sigma1=log(tan((1/2)*((pi/2)+lat1)));
Sigma2=log(tan((1/2)*((pi/2)+lat2)));
SigmaDiff=wrapRange(Sigma2-Sigma1,-pi,pi);

latDiff=wrapRange(lat2-lat1,-pi/2,pi/2);
lonDiff=wrapRange(lon2-lon1,-pi,pi);
%Equation 2 in [1], getting rid of the ambiguity by using atan2.
azimuth=atan2(lonDiff,SigmaDiff);

c=r*sqrt(SigmaDiff^2+lonDiff^2);
%Equation 5 of [1] for the distance along the rhumb line. This is not
%good if the azimuth is near pi/2.
dist=c/abs(SigmaDiff/latDiff);

%If abs(azimuth) is very close to pi/2, then dist will be a NaN due to
%SigmaDiff and latDiff both being zero. In such an instance, we put in
%the asumptotic limit of the ratio as the numerator and denominator go
%to zero.
sel=~isfinite(dist);
dist(sel)=c*cos(lat1);
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
