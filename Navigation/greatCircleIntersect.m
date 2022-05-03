function [latLonPoint,sAX,sCX]=greatCircleIntersect(az1,az2,latLon1,latLon2,r,getBoth)
%%GREATCIRCLEINTERSECT Earthbound sensors 1 and 2 located at latLon1 and
%           latLon2 and measure the azimuthal directions of a target given
%           in radians East of North, az1 and az2. This function determines
%           the latitude and longitude of the target on a spherical Earth.
%           The sensors cannot be placed at the North or South poles. Note
%           that the location of a target on a line running directly
%           between the sensors is unobservable. This type of localization
%           is also known as a "cross fix" or a "crossfix".
%
%INPUTS: az1 The direction of the target in radians East of North as
%            measured by the first sensor.
%        az2 The direction of the target in radians East of North as
%            measured by the second sensor.
%    latLon1 The 2X1 location of the first sensor given in latitude and
%            longitude in radians in the format [latitude;longitude].
%            This should not be the North or South pole. This is with
%            respect to the reference sphere.
%    latLon2 The 2X1 location of the second sensor given in latitude and
%            longitude in radians in the format [latitude;longitude].
%            This should not be the North or South pole.
%          r The assumed radius of the spherical Earth model. This is
%            needed for the sAX and sCX outputs to have the proper scale.
%            If omitted or an empty matrix is passed, the default of
%            osculatingSpher4LatLon(latLon1) is used.
%    getBoth By default, the solution returned is the one that is closest
%            to latLon1. However, an additional solution on the other side
%            of the Earth exists. If this is true, then both solutions will
%            be returned. The default if omitted or an empty matrix is
%            passed is false.
%
%OUTPUTS: latLonPoint The 2X1 [latitude;longitude] location in radians of
%                     the target (the closest solution to latLon1). If
%                     getBoth is true, then this is a 2X2 set of points,
%                     one for the near and the far solution.
%                 sAX The distance from sensor 1 to the target. If getBoth
%                     is true, this is 1X2 for both solutions.
%                 sCX The distance from sensor 2 to the target. If getBoth
%                     is true, this is 1X2 for both solutions.
%
%Great circles are geodesics on a sphere. A geodesic path between two
%points on a sphere is much simpler to determine as compared to an
%ellipsoid (which is what one gets using the indirectGeodeticProb
%function). For the spherical model, the latitude is taken as the elevation
%and the longitude as the azimuth on the sphere (we assume coordinate
%system 0 in the spher2Cart function).
%
%This function implements the algorithm used for initialization given in
%[1]. The A and C angles in the triangle in [1] are computed based on
%differencing the azimuthal angles toward the target and the azimuthal
%angles of the great circle route between the sensors. Additionally, under
%the assumptions inherent in how A and C are obtained, when sAX or sCX are
%negative, we have to add r*pi to get the correct positive values to use.
%
%EXAMPLE:
%This has the two sensors across the Pacific Ocean on an ideal spherical
%Earth. The ideal spherical model is used to demonstrate that in the ideal,
%noise-free case, the algorithm produces the correctly localized result.
%The spherical triangle being solves is also plotted on a spherical Earth.
%Note that using the same radius for each case is important to get
%consistent results.
% N=100;
% latLonA=[34.685169;139.443632]*(pi/180);
% latLonC=[-33.8617;151.2117]*(pi/180);
% latLonX=[37.7917;-122.4633]*(pi/180);
% 
% %Get the azimuth
% [azAX,~,~,WPAX]=indirectGreatCircleProb(latLonA,latLonX,1,N);
% [azCX,~,~,WPCX]=indirectGreatCircleProb(latLonC,latLonX,1,N);
% [~,~,~,WPAC]=indirectGreatCircleProb(latLonA,latLonC,1,N);
% 
% %Convert the waypoints to Cartesian to plot.
% WPAX=ellips2Cart([WPAX;zeros(1,N+2)],1,0);
% WPAC=ellips2Cart([WPAC;zeros(1,N+2)],1,0);
% WPCX=ellips2Cart([WPCX;zeros(1,N+2)],1,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],1,0);
% plot3(WPAX(1,:),WPAX(2,:),WPAX(3,:),'-r','linewidth',4)
% plot3(WPAC(1,:),WPAC(2,:),WPAC(3,:),'-g','linewidth',4)
% plot3(WPCX(1,:),WPCX(2,:),WPCX(3,:),'-b','linewidth',4)
% view(-90,0);
% 
% latLonPointX=greatCircleIntersect(azAX,azCX,latLonA,latLonC,1);
% max(abs(latLonPointX-latLonX))
%One will see that the points agree within finite precision limits.
%
%REFERENCES:
%[1] S. Baselga, J. C. Martinez-Llario, "Intersection and point-to-line
%    solutions for geodesics on the ellipsoid," Studi Geophysica et
%    Geodaetica, vol. 62, no. 3, pp. 353-363, Jul. 2018.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(getBoth))
    getBoth=false;
end

if(nargin<5||isempty(r))
    r=osculatingSpher4LatLon(latLon1);
end

%The great circle from latLon1 to latLon2.
[azStart,~,azEnd]=indirectGreatCircleProb(latLon1,latLon2,r);
%If azStart is the angle to the inside of the triangle, we want azEnd to be
%the angle to the inside of the triangle too (not the outside), so it has
%to go in the opposite direction along A-C. azStart is already in the range
%of -pi to pi.

azEnd=wrapRange(azEnd+pi,-pi,pi);
%azStart and azEnd are in the range of -pi to pi.
az1=wrapRange(az1,-pi,pi);
az2=wrapRange(az2,-pi,pi);

A=abs(azStart-az1);
if(A>pi)
    A=2*pi-A;
end

C=abs(azEnd-az2);
if(C>pi)
    C=2*pi-C;
end

%Get unit vectors to the points on the unit sphere.
uA=spher2Cart([latLon1(2);latLon1(1)]);
uC=spher2Cart([latLon2(2);latLon2(1)]);

%This is the spherical distance across the surface of the unit sphere
%between points 1 and 2. Essentially, get the angle (in radians) between
%the points and then multiply by the radius (here 1) to get the arclength.
sAC=angBetweenVecs(uA,uC);

sinsAC=sin(sAC);
cossAC=cos(sAC);

%Equation 3 in [1].
sAX=r*atan(sinsAC/(cossAC*cos(A)+sin(A)/tan(C)));

%Go in the direction of observation, not reverse to it.
if(sAX<0)
    sAXAlt=sAX;
    sAX=sAX+r*pi;
else
    sAXAlt=sAX-r*pi;
end

%Equation 4 in [1].
sCX=r*atan(sinsAC/(cossAC*cos(C)+sin(C)/tan(A)));

%Both latLon points obtained should be the same, except for finite
%precision errors.
latLonPoint=directGreatCircleProb(latLon1,az1,sAX,r);

if(getBoth)
    sAX=[sAX,sAXAlt];
    latLonPoint(:,2)=directGreatCircleProb(latLon1,az1,sAXAlt,r);
end

if(nargout>2)
    %If sCX is desired, we need to make sure that we get the correct sCX.
    %Over long distances, the first intersection point chosen by sCX might
    %not be the same as that chosen by sAX.

    %This should be close to latLonPoint, unless it is the wrong solution and
    %is on the other side of the Earth.
    latLonPointB=directGreatCircleProb(latLon2,az2,sCX,r);
    
    if(max(abs(latLonPointB-latLonPoint))<0.1)
        %If the ordering agrees.
        if(getBoth)
            if(sCX<0)
                sCXAlt=sCX+r*pi;
            else
                sCXAlt=sCX-r*pi;
            end
            sCX=[sCX,sCXAlt];
        end
    else
        if(sCX<0)
            sCXAlt=sCX;
            sCX=sCX+r*pi;
        else
            sCXAlt=sCX-r*pi;
        end
        if(getBoth)
           sCX=[sCX,sCXAlt]; 
        end
    end
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
