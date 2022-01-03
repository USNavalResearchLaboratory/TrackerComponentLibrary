function [latLonPoint1,latLonPoint2]=greatCircleTDOALoc(TDOA1,TDOA2,latLonRef,latLon1,latLon2,r,c)
%%GREATCIRCLETDOALOC Given a reference emitter and two other emitters being
%            observed by an object wishing to localize itself, OR given a
%            reference sensor and two other sensors wishing to localize an
%            emitter, determine the location of the object/emitter using
%            time-delay-of-arrival (TDOA) measurements under the assumption
%            of surface-wave propagation on a spherical Earth. This is the
%            type of localizing formerly used in the Long Range Navigation
%            (LORAN) systems.
%
%INPUTS: TDOA1 The scalar time of arrival at sensor 1 minus the time of
%              arrival at the reference.
%        TDOA2 The scalar time of arrival at sensor 2 minus the time of
%              arrival at the reference.
% latLonRef,latLon1,latLon2 The 2X1 latitudes and longitudes in radians in
%              the format [latitude;longitude] of the reference, sensor 1
%              and sensor 2. These are with respect to the assumed
%              reference sphere.
%            r The assumed radius of the spherical Earth model. If omitted
%              or an empty matrix is passed, the default of
%              r=osculatingSpher4LatLon(latLonRef) is used.
%            c The assumed speed of light. If this value is omitted or an
%              empty matrix is passed, then the default of
%              Constants.speedOfLight is used. Note that the speeds of
%              light for surface wave propagation and in the air are less
%              than the speed of light in a vacuum.
%
%OUTPUTS: latLonPoint1,latLonPoint2 The 2X1 solutions in latitude and
%              longitude with respect to the reference sphere.
%
%This function implements the algorithm of [1]. Note that in [1], it is
%suggested that one use an osculating sphere in the vicinity of what is
%being observed. The osculatingSpher4LatLon function can be used to obtain
%such a sphere.
%
%EXAMPLE 1:
%This example illustrates the geometry of how there can be two solutions.
%The reference sensor is plotted in red and sensors 1 and two are plotted
%in green and blue. The distances from the reference sensor to the other
%sensors are given as a red dashed line. The "true" distances to the target
%(plotted in yellow) are drawn in black. The problem is solved. In this
%case, the first solution is correct and we verify that the difference in
%the latitude and lognitudes found is within finite precision errors. We
%then plot the second solution in cyan and draw magenta lines from the
%sensors to it. Though it is farther away, the TDOAs are the same, which is
%verified to be within finite precision bounds.
% N=200;
% r=Constants.WGS84MeanRadius;%Spherical Earth radius.
% c=Constants.speedOfLight;
% latLonR=[34.685169;139.443632]*(pi/180);
% latLon1=[37.7917;-122.4633]*(pi/180);
% latLon2=[-33.8617;151.2117]*(pi/180);
% %The sensor/target location.
% latLonS=[-40;-110]*(pi/180);
%      
% xyzR=ellips2Cart([latLonR;0],r,0);
% xyz1=ellips2Cart([latLon1;0],r,0);
% xyz2=ellips2Cart([latLon2;0],r,0);
% xyzS=ellips2Cart([latLonS;0],r,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],r,0);
% scatter3(xyzR(1),xyzR(2),xyzR(3),200,'r','filled')
% scatter3(xyz1(1),xyz1(2),xyz1(3),200,'g','filled')
% scatter3(xyz2(1),xyz2(2),xyz2(3),200,'b','filled')
% scatter3(xyzS(1),xyzS(2),xyzS(3),200,'y','filled')
% view(-90,0)
% 
% %Get the distances between the sensors and also the target. They are
% %angles scaled by r and are thus just labeled as angles.
% [~,~,~,latLonWaypointsm1]=indirectGreatCircleProb(latLonR,latLon1,r,N);
% [~,~,~,latLonWaypointsm2]=indirectGreatCircleProb(latLonR,latLon2,r,N);
% [~,thetams,~,latLonWaypointsms]=indirectGreatCircleProb(latLonR,latLonS,r,N);
% [~,theta1s,~,latLonWaypoints1s]=indirectGreatCircleProb(latLon1,latLonS,r,N);
% [~,theta2s,~,latLonWaypoints2s]=indirectGreatCircleProb(latLon2,latLonS,r,N);
% 
% %Convert the waypoints to Cartesian to plot.
% latLonWaypointsm1=ellips2Cart([latLonWaypointsm1;zeros(1,N+2)],r,0);
% latLonWaypointsm2=ellips2Cart([latLonWaypointsm2;zeros(1,N+2)],r,0);
% latLonWaypointsms=ellips2Cart([latLonWaypointsms;zeros(1,N+2)],r,0);
% latLonWaypoints1s=ellips2Cart([latLonWaypoints1s;zeros(1,N+2)],r,0);
% latLonWaypoints2s=ellips2Cart([latLonWaypoints2s;zeros(1,N+2)],r,0);
% 
% plot3(latLonWaypointsm1(1,:),latLonWaypointsm1(2,:),latLonWaypointsm1(3,:),'--r','linewidth',4)
% plot3(latLonWaypointsm2(1,:),latLonWaypointsm2(2,:),latLonWaypointsm2(3,:),'--r','linewidth',4)
% plot3(latLonWaypointsms(1,:),latLonWaypointsms(2,:),latLonWaypointsms(3,:),'-k','linewidth',4)
% plot3(latLonWaypoints1s(1,:),latLonWaypoints1s(2,:),latLonWaypoints1s(3,:),'-k','linewidth',4)
% plot3(latLonWaypoints2s(1,:),latLonWaypoints2s(2,:),latLonWaypoints2s(3,:),'-k','linewidth',4)
% 
% %The error-free measurements.
% TDOA1=(theta1s-thetams)/c;
% TDOA2=(theta2s-thetams)/c;
% [sol1,sol2]=greatCircleTDOALoc(TDOA1,TDOA2,latLonR,latLon1,latLon2,r,c);
% errCur=min(max(abs(sol1-latLonS)),max(abs(sol2-latLonS)))
% %errCur will be within finite precision limits.
% xyzSAlt=ellips2Cart([sol2;0],r,0);
% %To see where the second solution is, one will have to rotate the globe.
% scatter3(xyzSAlt(1),xyzSAlt(2),xyzSAlt(3),200,'c','filled')
% 
% [~,thetamsAlt,~,latLonWaypointsmsAlt]=indirectGreatCircleProb(latLonR,sol2,r,N);
% [~,theta1sAlt,~,latLonWaypoints1sAlt]=indirectGreatCircleProb(latLon1,sol2,r,N);
% [~,theta2sAlt,~,latLonWaypoints2sAlt]=indirectGreatCircleProb(latLon2,sol2,r,N);
% 
% TDOA1Alt=(theta1sAlt-thetamsAlt)/c;
% TDOA2Alt=(theta2sAlt-thetamsAlt)/c;
% 
% TDOA1-TDOA1Alt
% TDOA2-TDOA2Alt
% %The differences in the TDOAs between solutions will be within finite
% %precision limits.
% %Convert the waypoints to Cartesian to plot.
% latLonWaypointsmsAlt=ellips2Cart([latLonWaypointsmsAlt;zeros(1,N+2)],r,0);
% latLonWaypoints1sAlt=ellips2Cart([latLonWaypoints1sAlt;zeros(1,N+2)],r,0);
% latLonWaypoints2sAlt=ellips2Cart([latLonWaypoints2sAlt;zeros(1,N+2)],r,0);
% plot3(latLonWaypointsmsAlt(1,:),latLonWaypointsmsAlt(2,:),latLonWaypointsmsAlt(3,:),'--m','linewidth',4)
% plot3(latLonWaypoints1sAlt(1,:),latLonWaypoints1sAlt(2,:),latLonWaypoints1sAlt(3,:),'--m','linewidth',4)
% plot3(latLonWaypoints2sAlt(1,:),latLonWaypoints2sAlt(2,:),latLonWaypoints2sAlt(3,:),'--m','linewidth',4)
%
%EXAMPLE 2:
%Figure 5 in [1] shows a particular geometry of the sensors and the target.
%What if the geometry changes? Consider loc1, loc2, and loc3 to be the
%reference and sensors 1 and 2. loc4List list holds locations of the target
%that are the same as in Fig. 5 and that differ (for example, behind all of
%the sensors). We go through all sensor locations in loc4List AND also go
%through all permutations of the points (e.g. switch sensors 1 and 2 or the
%reference, or switch sensor 1 with the target). Taking both solutions, we
%look at the minimum error solution, taking the maximum over all
%permutations of points. We display this value, which is seen to be within
%finite precision limits. That shows that this function works with many
%geometries beyond the one shown in Fig. 5 in [1].
% r=Constants.WGS84MeanRadius;%Spherical Earth radius.
% c=Constants.speedOfLight;
% loc1=[34.685169;139.443632]*(pi/180);
% loc2=[37.7917;-122.4633]*(pi/180);
% loc3=[-33.8617;151.2117]*(pi/180);
% 
% %For the first scenario.
% loc4List=[-40,   1,  50,  50, -20;
%          -110, 165, 110, 170, 120]*(pi/180);
% numS=size(loc4List,2);
%      
% idxPerms=genAllPermutations(4);
% numPerms=size(idxPerms,2);
% 
% maxErr=0;
% for curLast=1:numS
%     latLonList=[loc1,loc2,loc3,loc4List(:,curLast)];
% 
%     for curPerm=1:numPerms
%         thePerm=idxPerms(:,curPerm);
%         latLonListCur=latLonList(:,thePerm);
% 
%         latLonR=latLonListCur(:,1);
%         latLon1=latLonListCur(:,2);
%         latLon2=latLonListCur(:,3);
%         latLonS=latLonListCur(:,4);
% 
%         %Get the distances between the target and the sensors. They are
%         %angles scaled by r and are thus just labeled as angles.
%         [~,thetams]=indirectGreatCircleProb(latLonR,latLonS,r);
%         [~,theta1s]=indirectGreatCircleProb(latLon1,latLonS,r);
%         [~,theta2s]=indirectGreatCircleProb(latLon2,latLonS,r);
% 
%         %The error-free measurements.
%         TDOA1=(theta1s-thetams)/c;
%         TDOA2=(theta2s-thetams)/c;
% 
%         [sol1,sol2]=greatCircleTDOALoc(TDOA1,TDOA2,latLonR,latLon1,latLon2,r,c);
%         errCur=min(max(abs(sol1-latLonS)),max(abs(sol2-latLonS)));
%         maxErr=max(maxErr,errCur);
%     end
% end
% maxErr
%maxErr will be within finite precision limits.
%
%REFERENCES:
%[1] P. Williams and D. Last, "On Loran-C time-difference to co-ordinate
%    converters," in Proceedings of the 32nd Annual Convention & Technical
%    Symposium of the International Loran Association, Boulder, CO, 3-7
%    Nov. 2003.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<6||isempty(r))
    r=osculatingSpher4LatLon(latLonRef); 
end

%Equation 3.21 in [1].
P1=(c/r)*TDOA1;
P2=(c/r)*TDOA2;

cosP1=cos(P1);
sinP1=sin(P1);
cosP2=cos(P2);
sinP2=sin(P2);

%Get the normalized great circle distance and heading from the reference to
%each of the other known locations.
[azStart1,thetam1]=indirectGreatCircleProb(latLonRef,latLon1,1);
[azStart2,thetam2]=indirectGreatCircleProb(latLonRef,latLon2,1);

cosThetam1=cos(thetam1);
sinThetam1=sin(thetam1);
cosThetam2=cos(thetam2);
sinThetam2=sin(thetam2);

%Inferred from Fig. 5 in [1].
K=abs(wrapRange(azStart2-azStart1,-pi,pi));
cosK=cos(K);
sinK=sin(K);

a1=(cosP1-cosThetam1)/sinThetam1;
a2=(cosP2-cosThetam2)/sinThetam2;

u1=a1*cosK-a2;
u2=a1*sinK;
u3=a2*(sinP1/sinThetam1)-a1*(sinP2/sinThetam2);

%Equation 3.26 in [1].
val1=u3*u1;
denom=u1^2+u2^2;
%The max function should not be needed, but it is there in case some finite
%precision issue arises.
val2=u2*sqrt(max(0,denom-u3^2));

cosBeta1=(val1+val2)/denom;
cosBeta1Alt=(val1-val2)/(u1^2+u2^2);

%Equation 3.24 in [1], taking the inverse.
thetams=atan((cosP1-cosThetam1)/(sinP1+sinThetam1*cosBeta1));
thetams=thetams+(thetams<0)*pi;
thetamsAlt=atan((cosP1-cosThetam1)/(sinP1+sinThetam1*cosBeta1Alt));
thetamsAlt=thetamsAlt+(thetamsAlt<0)*pi;

xyzRef=ellips2Cart([latLonRef;0],r,0);
xyz1=ellips2Cart([latLon1;0],r,0);
xyz2=ellips2Cart([latLon2;0],r,0);

%Equation 3.29 in [1].
A=[xyz1.';
   xyz2.';
   xyzRef.'];

latLonElPoint1=getSol4Thetams(thetams,cosP1,sinP1,cosP2,sinP2,A,r);
latLonElPoint2=getSol4Thetams(thetamsAlt,cosP1,sinP1,cosP2,sinP2,A,r);
latLonPoint1=latLonElPoint1(1:2);
latLonPoint2=latLonElPoint2(1:2);
end

function latLonElPoint=getSol4Thetams(thetams,cosP1,sinP1,cosP2,sinP2,A,r)
%%GETSOL4THETAMS Get a solution for a particular thetams. This is broken
%                out as a separate function, since the work has to be done
%                twice, once for each solution.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

cosThetams=cos(thetams);
sinThetams=sin(thetams);

%Equation 3.23 in [1].
cosTheta1s=cosThetams*cosP1-sinThetams*sinP1;
cosTheta2s=cosThetams*cosP2-sinThetams*sinP2;

%Get the Cartesian location of the item to be localized.
%Equation 3.30 in [1].
b=r^2*[cosTheta1s;
       cosTheta2s;
       cosThetams];
%Equation 3.31 in [1].
sLoc=A\b;

latLonElPoint=Cart2Ellipse(sLoc,[],r,0);

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
