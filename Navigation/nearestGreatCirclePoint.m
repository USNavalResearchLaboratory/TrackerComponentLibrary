function [latLonX,s1X,s2X]=nearestGreatCirclePoint(az1,latLon1,latLon2,r,getBoth)
%%NEARESTGRETCIRCLEPOINT Given a point at latLon1, a heading East of North
%           az1, and a second point at latLon2, find the point on on the
%           trajectory that is closest to latLon2 when considering a
%           spherical Earth. The closest point is such that a line from
%           latLon2 to the point intersects the trajectory at a right
%           angle. This point might be backwards on the trajectory (in the
%           az+pi direction). A second solution where the lines intersect
%           at a right angle (on the sphere) can also be obtained, if
%           desired. The second solution is always in the opposite
%           direction along the path (forward or reverse az direction) and
%           is halfway around the world from the first solution.
%
%INPUTS: az1 The azimuthal direction of the trajectory at point latLon1.
%            This is given in radians East of North.
%    latLon1 The reference point defining the trajectory. This is
%            [latitude;longitude] on the reference sphere given in radians.
%            the longitude is East longitude (West is negative).
%    latLon2 A second point whereby the shortest distance to the trajectory
%            of the first point is desired.
%          r The assumed radius of the spherical Earth model. If omitted or
%            an empty matrix is passed, the default of
%            osculatingSpher4LatLon(latLon1) is used.
%    getBoth Specified whether both the near and far solutions are desired.
%            The default if omitted or an empty matrix is passed is false.
%
%OUTPUTS: latLonX A 2X1 [latitude;longitude] point on the trajectory that
%                 is closest to latLon2. If getBoth is true, then this is
%                 2X2 and the second solution is the far point.
%             s1X The 1X1 or 1X2 great circle distances  distance between
%                 latLon1 and latLonX.
%             s2X The 1X1 or 1X2 great circle distances  distance between
%                 latLon2 and latLonX.
%
%This function implements the algorithm described at the beginning of
%Section 3 of [1], which plays a role in a more sophisticated approach for
%ellipsoidal models.
%
%EXAMPLE:
%Here, we have a point A with an azimuthal motion direction azA. We also
%have a point P. We want to find the point X on the trajectory described by
%A that is closest to the point P. The solution is a point PX that is at a
%90 degree angle to AX. We display both solutions (the closest and the
%other solution that forms a 90 degree angle. The second solution is always
%halfway around the world from the first solution. Note that the nearest
%solution is not always in the direction of azA -- it could be in the
%opposite direction (backwards on the path). Locations and distances are
%normalized to the unit circle.
% N=500;
% latLonA=[34.685169;139.443632]*(pi/180);
% azA=54*(pi/180);
% latLonP=[-33.8617;-151.2117]*(pi/180);
% [~,~,~,WPAP]=indirectGreatCircleProb(latLonA,latLonP,1,N);
% 
% distVals=linspace(0,2*pi,N+2);
% aDirPts=directGreatCircleProb(latLonA,azA,distVals,1);
% 
% latLonX=nearestGreatCirclePoint(azA,latLonA,latLonP,1,true);
% [azXP1,~,~,WPXP1]=indirectGreatCircleProb(latLonX(:,1),latLonP,1,N);
% [azXA1]=indirectGreatCircleProb(latLonX(:,1),latLonA,1);
% 
% [azXP2,~,~,WPXP2]=indirectGreatCircleProb(latLonX(:,2),latLonP,1,N);
% [azXA2]=indirectGreatCircleProb(latLonX(:,2),latLonA,1);
% 
% %These quantities should both be 90 -showing that the intersection angles
% %of both solutions with the trajectory of A are 90 degrees.
% abs(wrapRange((azXP1-azXA1),-pi,pi))*(180/pi)
% abs(wrapRange((azXP2-azXA2),-pi,pi))*(180/pi)
% 
% %Convert the waypoints to Cartesian to plot. Plot slightly above the
% %ellipsoid so that the line can be seen better.
% WPAP=ellips2Cart([WPAP;0.02*ones(1,N+2)],1,0);
% aDirPts=ellips2Cart([aDirPts;0.02*ones(1,N+2)],1,0);
% WPXP1=ellips2Cart([WPXP1;0.02*ones(1,N+2)],1,0);
% WPXP2=ellips2Cart([WPXP2;0.02*ones(1,N+2)],1,0);
% 
% figure(1)
% clf
% hold on
% h=plotMapOnEllipsoid([],1,0);
% plot3(WPAP(1,:),WPAP(2,:),WPAP(3,:),'-r','linewidth',4)
% plot3(aDirPts(1,:),aDirPts(2,:),aDirPts(3,:),'-g','linewidth',4)
% plot3(WPXP1(1,:),WPXP1(2,:),WPXP1(3,:),'-b','linewidth',4)
% plot3(WPXP2(1,:),WPXP2(2,:),WPXP2(3,:),'-c','linewidth',4)
% view(-60,0)
% set(h,'HandleVisibility','off');%Omit from the legend.
% legend('Point to Trajectory Start','Trajectory','Closest Solution','Other Solution')
%
%REFERENCES:
%[1] S. Baselga, J. C. Martinez-Llario, "Intersection and point-to-line
%    solutions for geodesics on the ellipsoid," Studi Geophysica et
%    Geodaetica, vol. 62, no. 3, pp. 353-363, Jul. 2018.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(getBoth))
    getBoth=false;
end

if(nargin<4||isempty(r))
    r=osculatingSpher4LatLon(latLon1);
end

%Everything is normalized to the unit sphere until the end.
[az12,s12]=indirectGreatCircleProb(latLon1,latLon2,1);

%The inner angle of the triangle.
az1=wrapRange(az1+pi,-pi,pi);

A=abs(az12-az1);
if(A>pi)
    A=2*pi-A;
end

if(getBoth)
    s12=[s12,s12+pi];
    latLonX=zeros(2,2);
else
    latLonX=zeros(2,1);
end

%Get the closest solution.
%Equation 8 in [1].
s2X=asin(sin(s12)*sin(A));
%Equation 10 in [1], simplified to turn the sine ratio into a tangent.
s1X=2*atan(tan((1/4)*(2*A+pi))*tan((s12-s2X)/2));
latLonX(:,1)=directGreatCircleProb(latLon1,az1,s1X(1),1);
if(getBoth)
    latLonX(:,2)=directGreatCircleProb(latLon1,az1,s1X(2),1);
end

s1X=abs(r*s1X);
s2X=abs(r*s2X);

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
