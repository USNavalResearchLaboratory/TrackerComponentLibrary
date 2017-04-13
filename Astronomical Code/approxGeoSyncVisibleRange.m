function [lonRange,maxEl]=approxGeoSyncVisibleRange(latlonObs,a,f,rG)
%%APPROXGEOSYNCVISIBLERANGE At a particular location on the Earth,
%               a certain range of geosynchronous equatorial satellites
%               (geostationary statellites) might be visible. Here, we use
%               the approximation that the Earth is an ellipsoid and that
%               geosynchronous satellites are in a perfectly circular orbit
%               around the equator of the Earth to quickly determine what
%               range of longitudes of satellites are visible to the
%               observer. The solution is only geometric and does not
%               consider that atmospheric refraction allows one to see
%               below the geometric horizon.
%
%INPUTS: latlonObs A 2X1 point given in ellipsoidal latitude and longitude
%                  coordinates in radians. The point will be taken to be on
%                  the surface of a reference ellipsoid (height=0).
%                a The semi-major axis of the reference ellipsoid. If this
%                  argument is omitted or an empty matrix is passed, the
%                  value in Constants.WGS84SemiMajorAxis is used.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted or an empty matrix is passed,
%                  the value in Constants.WGS84Flattening is used.
%               rG The orbital radius of geosynchronous satellites. If this
%                  parameter is omitted or an empty matrix is passed, then
%                  a standard value (derived below) is used having units of
%                  meters (which should match a).
%
%OUTPUTS: lonRange A 2X1 vector where the longitudes (in radians) over
%                  which gepsynchronous satellites are visible. If none are
%                  visible, then an empty matrix is returned. Near the
%                  poles, none will be visible.
%            maxEl The maximum elevation angle in radians (with respect to
%                  the local tangent plane of the reference ellipsoid at
%                  point latlonObs) that a geosynchronous satellite will be
%                  visible. If non are visible, then an empty matrix is
%                  returned. If any satellites are visible, then the
%                  minimum elevation will be zero.
%
%This function is derived by taking a geosynchronous orbit as a circular
%orbit with the same period as the rotation of the Earth, aligned with the
%equator. The Earth is taken to be an ellipsoid. Satellites are visible if
%they are not below the local tangent plane to the ellipsoid. The tangent
%plane is defined as  the plane orthogonal to a local normal to the
%ellipsoid.
%
%To begin, we will determine the value rG, which is the radius of a
%geostationary orbit. Approximating the gravitation force of the Earth as a
%point with mass M, Newton's Universal Law of Gravitation tells us that the
%magnitude of acceleration due to gravity is
%norm(g)=G*M/r^2
%where r is the difference from the center of mass and G is the universal
%gravitational constant. On the other hand, the magnitude of centripetal
%acceleration of a rotating body is 
%norm(a)=omega^2*r
%where omega is the radis of rotation. it might seem odd to equate a
%gravitational acceleration with a centripetal acceleration, but that is
%exactly one method of deriving Newton's law of gravitation from Kepler's
%equations. Thus, here, we set the two equations equal and solve for the
%radius to get
%r^3=G*M/omega^2
%This is how one can getermine the geostationary orbit radius. For the
%parameters, if rG is omitted, we use G*M=Constants.WGS84GMWithAtmosphere
%and omega=Constants.WGS84EarthRotationRate.
%
%The reference ellipsoid of the Earth has the equator in the x-y plane and
%the z plane going approximately through the rotational axis of the Earth.
%The Cartesian location xyzO of the point latlonObs on the surface of the
%reference ellipsoid can be determined using the ellips2Cart function. A
%normal vector u to the surface of the reference ellipsoid can be obtained
%using the getENUAxes function. A point along the circle of possible
%geosynchronous satellites is xyzG=[cos(lambda);sin(lambda);0]. A
%geostationary satellite is visible if the angle between the unit vector u
%and the vector xyzG-xyzO is less than pi/2. That is, the satellite is
%above the horizon. We can determine the angle between the vectors using a
%dot product. For two vectors a and b it is known that
%dot(a,b)=norm(a)*norm(b)*cos(theta)
%where theta is the angle separating vectors a and b. when the dot product
%is zero, the vectors are orthogonal. Thus, solving dot(u,xyzG-xyzO)==0 for
%lambda gives the two solutions to the problem. That is what is done here
%to determine lonRange.
%
%On the other hand, to determine the maximum elevation above the horizon,
%we note that for a target at zero longitude, the maximum elevation is
%obtained with respect to a geosynchronous satellite at zero longitude. Due
%to symmetry, changing the longitude of latlonObs does not change the
%maximum. Thus, we compute the maximum elevation as
% maxEl=pi/2-angBetweenVecs([rG;0;0]-xyzAlt,uAlt)
%where xyzAlt is the location at the specified lattude and zero longitude
%on the reference ellipsoid and ualt is the normal to the reference
%ellipsoid at that point. [rG;0;0] is the location of the highest elevation
%angle geostationary satellite at that angle.
%
%EXAMPLE 1:
%Here, we verify that on the equator, the maximum elevation angle is
%directly overhead and one can see a bit less than 180 degrees.
% [lonRange,maxEl]=approxGeoSyncVisibleRange([0;0])
%We now consider how little ones sees at a high latitude
% [lonRange,maxEl]=approxGeoSyncVisibleRange([81*(pi/180);0])
%At 81 degrees latitude, the maximum elevation angle is about 0.328 degrees
%(0.00573 radians) and one can only see a span of latitudes of 0.5396
%radians, which is about 30.9 degrees.
%Finally, at 82 degrees, empty matrices are returned, because geosynchrnous
%satellites are no longer visible.
%[lonRange,maxEl]=approxGeoSyncVisibleRange([82*(pi/180);0])
%
%EXAMPLE 2:
%Here, we plot maxEl as a function of latitude. Points below the horizon
%are simply set to zero elevation.
% numPoints=500;
% lat=linspace(-pi/2,pi/2,numPoints);
% maxEl=zeros(numPoints,1);
% for curPoint=1:numPoints
%     [~,retVal]=approxGeoSyncVisibleRange([lat(curPoint);0]);
%     if(~isempty(retVal))
%         maxEl(curPoint)=retVal;
%     end
% end
% figure(1)
% clf
% plot(lat*(180/pi),maxEl*(180/pi),'linewidth',2)
% axis([-90,90,0,90])
% h1=xlabel('Latitude (degrees)');
% h2=ylabel('Maximum GEO elevation (degrees)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 3:
%Here, we plot the width of the longitudinal angular viewing region of
%geosynchronous satellites. We center the viewing region at zero lattitude
%so that a simple diff function gives us the span without worrying about
%the circular nature of longitude.
% numPoints=500;
% lat=linspace(-pi/2,pi/2,numPoints);
% spanVals=zeros(numPoints,1);
% for curPoint=1:numPoints
%     lonRange=approxGeoSyncVisibleRange([lat(curPoint);0]);
%     if(~isempty(lonRange))
%         spanVals(curPoint)=diff(lonRange);
%     end
% end
% figure(2)
% clf
% plot(lat*(180/pi),spanVals*(180/pi),'linewidth',2)
% axis([-90,90,0,180])
% h1=xlabel('Latitude (degrees)');
% h2=ylabel('Visible longitudinal span (degrees)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<2||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(rG))
    GM=Constants.WGS84GMWithAtmosphere;
    omega=Constants.WGS84EarthRotationRate;
    rG=nthroot(GM/omega^2,3);
end

u=getENUAxes(latlonObs,true,a,f);
xyzCart=ellips2Cart([latlonObs;0],a,f);

c=dot(u,xyzCart)/rG;

%There are issues if y is too small. here, we use an asymptotic
%approximation for y=0 when y becomes very small.
if(abs(u(2))<2^(-1024))
    %Use the asymptotic solution for y=0.
    lambda1=-acos(c/u(1));
    lambda2=-lambda1;
    
    if(~isreal(lambda1))%Mark as having no real solution.
        lambda1=[];
        lambda2=[];
    end
else%No asymptotics needed.
    uxyMag2=sum(u(1:2).^2);
    radical=abs(u(2))*sqrt(uxyMag2-c^2);
    
    if(~isreal(radical))
        %Mark as having no real solution.
        lambda1=[];
        lambda2=[];
    else%Compute the real solutions.
        denom2=u(1)^2*u(2)+u(2)^3;

        val1=(c*u(2)^2+u(1)*radical)/denom2;
        lambda1=atan2(val1,(c*u(1)-radical)/uxyMag2);

        val1=(c*u(2)^2-u(1)*radical)/denom2;
        lambda2=atan2(val1,(c*u(1)+radical)/uxyMag2);
    end
end

if(~isempty(lambda1))%Real solutions exist.
    lonRange=[lambda1;lambda2];
else%No real solutions exist.
    lonRange=[];
end

maxEl=[];
if(nargout>1&&~isempty(lonRange))
    %If the maximum elevation of a geosynchronous satellite above the local
    %tangent plane is desired.
    shiftedLatLonAlt=[latlonObs(1);0;0];
    
    xyzCartAlt=ellips2Cart(shiftedLatLonAlt,a,f);
    uAlt=getENUAxes(shiftedLatLonAlt,true,a,f);

    maxEl=pi/2-angBetweenVecs([rG;0;0]-xyzCartAlt,uAlt);
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
