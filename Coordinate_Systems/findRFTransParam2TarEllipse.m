function M=findRFTransParam2TarEllipse(plhRx,plhTar,zRot,a,f)
%%FINDRFTRANSPARAM2TARELLIPS Find the transformation matrix needed to
%            rotate a global Cartesian vector to local radar-facing
%            coordinates such that a specified target is on the boresight
%            of the radar. This makes the local z-axis of the radar point
%            at the target when the radar and target locations are
%            specified in ellipsoidal [latitude;longitude;height]
%            coordinates on the Earth.
%
%INPUTS: plhRx The [latitude; longitude; height] of the receiver with
%              respect to the reference ellipsoid. Latitude and longitude
%              are given in radians.
%       plhTar The [latitude; longitude; height] of the target with
%              respect to the reference ellipsoid. Latitude and longitude
%              are given in radians.
%         zRot This is an additional (counterclockwise) rotation in
%              radians about the pointing direction of the radar (the
%              z-axis). The default if omitted or an empty matrix is
%              passed is zero. If the receiver and transmitter are in the
%              same plane, then zRot=0 means that the local y axis points
%              up.
%            a The semi-major axis of the reference ellipsoid. If this
%              argument is omitted, the value in
%              Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted, the value in Constants.WGS84Flattening
%              is used.
%
%OUTPUTS: M A rotation matrix for the transformation from global Cartesian
%           coordinates to local radar-facing coordinates. The z axis
%           represents the pointing direction of the radar. This can be
%           directly fed into the RUV coordinate transform functions.
%
%Rotations are made with respect to the local tangent plane of the
%reference ellipsoid at the receiver location using the findRFTransParam
%function.
%
%This function differs from findRFTransParam2Tar not only in that the
%target and receiver locations are specified in Cartesian coordiante, but
%since the rotation are taken with respect to a local East-North-Up plane,
%the orientation of the non-z axes with zRot=0 differs as a function of
%position differently from what findRFTransParam2Tar does.
%
%EXAMPLE:
%We point a radar at a target and show that the direction components of
%r-u-v coordinates and spherical coordinates (systemType=1) are zero.
% %Receiver and target locations, latitude, longitude, ellipsoidal height.
% llRRx=[deg2rad([20.735566;-155.943136]);100];
% llRTar=[deg2rad([21.991950;-151.998054]);50e3];
% %Convert to Cartesian.
% lRx=ellips2Cart(llRRx);
% lTar=ellips2Cart(llRTar);
% 
% %Get a rotation matrix to point the radar at the 
% M=findRFTransParam2TarEllipse(llRRx,llRTar);
% 
% systemType=1;%Note this is required for spherical measurements.
% spherLoc=Cart2Sphere(lTar,systemType,false,lRx,lRx,M);
% ruvLoc=Cart2Ruv(lTar,false,lRx,lRx,M);
% %The spherical angles and u and v in ruv coordinates are all zero, meaning
% %the radar is pointing straight at the target.
% spherLoc(2:3)
% ruvLoc(2:3)
%
%February 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<4||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<3||isempty(zRot))
    zRot=0;
end

lTar=ellips2Cart(plhTar,a,f);

xTarENU=ECEF2ENU(plhRx,lTar,a,f);
xTarENUSpher=Cart2Sphere(xTarENU,0);

az=pi/2-xTarENUSpher(2);
el=xTarENUSpher(3);
M=findRFTransParam(plhRx,az,el,zRot,a,f);

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
