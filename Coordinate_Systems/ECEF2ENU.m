function [xENU,M]=ECEF2ENU(plhOrigin,xECEF,a,f)
%%ECEF2ENU Convert from an Earth-centered Earth-fixed (ECEF) Cartesian
%          coordinate system to a local East-North-Up (ENU) coordinate
%          system. The alignment of the coordinate axes with the reference
%          ellipsoid is the one used in common standards such as that used
%          by the DoD's WGS-84 standard and the International Earth
%          Rotation and Reference Systems Service's (IERS) international
%          terrestrial reference frame (ITRF). The global ECEF z-axis is
%          North.
%
%INPUTS: plhOrigin A 3X1 point in [latitude;longitude; ellipsoidal height]
%             coordinates with respect to a reference ellipsoid
%             parameterized by a and f, derving as the origin of the local
%             ENU coordinate system. Latitude and longitude are given in
%             radians North and East.
%       xECEF A 3XnumPts collection of 3D Cartesian points in an ECEF
%             coordinate system that should be converted to a local ENU
%             coordinate system.
%           a The semi-major axis of the reference ellipsoid. If this
%             argument is omitted, the value in
%             Constants.WGS84SemiMajorAxis is used.
%           f The flattening factor of the reference ellipsoid. If this
%             argument is omitted, the value in Constants.WGS84Flattening
%             is used.
%
%OUTPUTS: xENU A 3XnumPts matrix of the points in xECEF converted to ENU
%              coordinates.
%            M A rotation matrix such that rotates a point from ECEF to ENU
%              coordinates.
%
%The ENU coordinate system is a common local tangent plane cooridnate
%system. Here, we take "up" to be that described by a reference ellipsoid
%and not true gravitational "up". This limits the accuracy of the local
%coordinate system.
%
%The conversions are mentioned in [1].
%
%EXAMPLE:
%In this example, we convert a point from the local ENU coordinate system
%to ECEF and then back to the local ENU cordinate system. The relative
%error of the converted-back point is considered.
% lat=19.7241*(pi/180);
% lon=-155.0868*(pi/180);
% height=0;
% plhPoint=[lat;lon;height];%Around Hilo, Hawaii.
% 
% xENU=[100;200;300];
% RelErr=abs(ECEF2ENU(plhPoint,ENU2ECEF(plhPoint,xENU))-xENU)./xENU
%One will find the relative error to be under 1e-11, which is on the order
%of finite-precision limits.
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

%Make the reference point the origin.
origShiftPts=bsxfun(@minus,xECEF,originCart);

numPts=size(xECEF,2);
xENU=zeros(3,numPts);
%Extract the components of the points in ENU coordinates with the reference
%point at the new origin.
xENU(1,:)=sum(bsxfun(@times,uENU(:,1),origShiftPts),1);
xENU(2,:)=sum(bsxfun(@times,uENU(:,2),origShiftPts),1);
xENU(3,:)=sum(bsxfun(@times,uENU(:,3),origShiftPts),1);

if(nargout>1)
    M=uENU';
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
