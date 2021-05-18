function xECEF=NED2ECEF(plhOrigin,xNED,a,f)
%%NED2ECEF Convert from a local North-East-Down (NED) Cartesian cooridnate
%          system to an Earth-centered Earth-fixed (ECEF) Cartesian
%          coordinate system. The alignment of the coordinate axes with the
%          reference ellipsoid is the one used in common standards such as
%          that used by the DoD's WGS-84 standard and the International
%          Earth Rotation and Reference Systems Service's (IERS)
%          international terrestrial reference frame (ITRF). The global
%          ECEF z-axis is North.
%
%INPUTS: plhOrigin A 3X1 point in [latitude;longitude; ellipsoidal height]
%             coordinates with respect to a reference ellipsoid
%             parameterized by a and f, derving as the origin of the local
%             NED coordinate system. Latitude and longitude are given in
%             radians North and East.
%        xNED A 3XnumPts collection of 3D Cartesian points in the local NED
%             coordinate system that should be converted to a global ECEF
%             coordinate system.
%           a The semi-major axis of the reference ellipsoid. If this
%             argument is omitted, the value in
%             Constants.WGS84SemiMajorAxis is used.
%           f The flattening factor of the reference ellipsoid. If this
%             argument is omitted, the value in Constants.WGS84Flattening
%             is used.
%
%OUTPUTS: xECEF A 3XnumPts matrix of the points in xNED converted to ECEF
%              coordinates.
%
%The NED coordinate system is a common local tangent plane cooridnate
%system. Here, we take "up" to be that described by a reference ellipsoid
%and not true gravitational "up". This limits the accuracy of the local
%coordinate system.
%
%The conversions are mentioned in [1].
%
%EXAMPLE:
%Here, we convert 3 points in a local NED tangent-plane coordinate system
%to ECEF and then to ellipsoidal coordinates. In the NED coordinate system,
%the points are offset in the  North, East, and up directions. One will see
%that in ellipsoidal coordinates, as compared to the origin, there are
%respective offsets in latitude (with no change in longitude), longitude
%(with no change in latitude) and height (with no change in latitude or
%longitude).
% lat=19.7241*(pi/180);
% lon=-155.0868*(pi/180);
% height=0;
% plhPoint=[lat;lon;height];%Around Hilo, Hawaii.
% xNED=[1e3,   0,    0;
%         0, 1e3,    0;
%         0,   0, -1e3];
% xECEF=NED2ECEF(plhPoint,xNED);
% latLonHeight=Cart2Ellipse(xECEF);
% diffVals=bsxfun(@minus,latLonHeight,plhPoint);
% %Convert angles to degrees
% diffVals(1:2,:)=diffVals(1:2,:)*(180/pi);
% %One can see the change in latitude on the first point, the change in
% %longitude on the second point, and no change in laittude or longitude in
% %the third point.
% diffVals(1:2,:)
% %The change in height with the third point is the full kilometer. The
% %changes with the other points is just due to the tangent point
% %approximation and are much smaller.
% diffVals(3,:)
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

uENU=getENUAxes(plhOrigin,false,a,f);
originCart=ellips2Cart(plhOrigin,a,f);

numPts=size(xNED,2);
xECEF=zeros(3,numPts);
for k=1:numPts
    xECEF(:,k)=xNED(1,k)*uENU(:,2)+xNED(2,k)*uENU(:,1)-xNED(3,k)*uENU(:,3);
end

%Correct for the origin shift.
xECEF=bsxfun(@plus,xECEF,originCart);

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
