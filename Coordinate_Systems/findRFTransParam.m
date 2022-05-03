function M=findRFTransParam(plhPoint,az,el,zRot,a,f,useEUN)
%%FINDRFTRANSPARAMS Find the transformation matrix needed to rotate a
%            global Cartesian vector to local radar-facing coordinates.
%            With az=0, el=0 and zRot=0, the transformation matrix returned
%            by this function makes the x-axis point West, the y-axis point
%            up and the z-axis point North (unless useEUN is true), which
%            is a right-handed coordinate system. With el=0 and zRot=0, az
%            being nonzero performs a clockwise rotation of the z and x
%            axes in the local tangent place of the receiver. el being
%            nonzero then raises the z axis above the local tangent plane
%            (after the azimuthal rotation), and zRot being nonzero then
%            rotates the x-y axes about the z axis counterclockwise.
%
%INPUTS: plhPoint The point at which the axes are to be found given in
%               terms of [latitude;longitude] with the geodetic latitude
%               and longitude. The latitude should be between -pi/2 and
%               pi/2. A height term can also be included, but it will be
%               ignored. This point defines the local East-North-Up
%               coordinate system with respect to a reference ellipsoid.
%            az The radar azimuth, in radians. If not provided, a value of
%               zero is assumed. The azimuth is the radians East of true
%               North (clockwise), as defined by the reference ellipsoid,
%               of the pointing direction of the radar. The default if
%               omitted or an empty matrix is passed is zero.
%            el The radar elevation, in radians. If not provided, a value
%               of zero is assumed. This is the angle above the horizon
%               (local level --the tangent to the reference ellipsoid) of
%               the pointing direction of the radar. The default is omitted
%               or an empty matrix is passed is zero.
%          zRot This is an additional (counterclockwise) rotation in
%               radians about the pointing direction of the radar (the
%               z-axis). The default if omitted or an empty matrix is
%               passed is zero.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in Constants.WGS84Flattening
%               is used.
%        useEUN If this is true, then when all rotations are zero, the axes
%               will point East-Up-North instead of West-Up-North. The WUN
%               coordinate system is right-handed. The EUN coordinate
%               system is is left-handed. The default if omitted or an
%               empty matrix is passed is false.
%
%OUTPUTS: M A rotation matrix for the transformation from global Cartesian
%           coordinates to local radar-facing coordinates. The z axis
%           represents the pointing direction of the radar. This can be
%           directly fed into the RUV coordinate transform functions.
%
%Details of this transformation can be found in Chapter 5 of [1].
%
%EXAMPLE:
%Confusion sometimes arises regarding the  direction of the x and y axes
%when pointing a radar. Here, we try to show this through an example.
% at the San Juan Airport (SJU; ICAO: TJSJ)
% radarLocDegs=[ 18.437082;-66.001724; 0];%Latitude, longitude, altitude.
% radarLocRad=radarLocDegs*(pi/180);%Convert to radians.
% 
% targetLocDeg=radarLocDegs+[-0.2;0.2;1000];%To the South-East and 1km up.
% targetLocRad=[targetLocDeg(1:2)*pi/180;targetLocDeg(3)];%Convert to
%                                                         %radians.
% 
% %The radar is pointing East, not elevated at all.
% az=90*(pi/180);
% el=0;
% M=findRFTransParam(radarLocRad,az,el);%Rotation matrix, global to local
% 
% %Convert to Cartesian
% radarCart = ellips2Cart(radarLocRad);
% targetCart = ellips2Cart(targetLocRad);
% targetRuvw=Cart2Ruv(targetCart,false,radarCart,radarCart,M,true)
%The last three element of targetRuvw are a unit vector in the local
%coordinate system of the radar. The local coordinate system of the radar
%has the z-axis pointing East, the y-axis Up, and the x-axis North. The
%target is South-East. Thus, the unit vector has a negative x-component and
%positive y and z components. Many folks expect the x component to be
%positive, because, when standing on the ground facing in the direction of
%the radar, the target will be forward, up and to the right. However, the
%x-axis is to one's left. To make the x-axis to the right, one should set
%zRot to pi to make the y axis point down. Thus, using
%M=findRFTransParam(radarLocRad,az,el,pi);
%one will get a conversion with a positive x component in the local
%coordinate system (and now a negative y component). For both x and y
%components to be positive, use a -90 degree rotation (-pi/2 radians) to
%make the x axis point up and the y axis point right.
%
%REFERENCES:
%[1] B. L. Diamond and M. E. Austin, "The Aries Program: Coordinates,
%    Transformations, Trajectories and Tracking",  Massachusetts Institute
%    of Technology Lincoln Laboratory, no. 1975-15, 5 Sept. 1975.
%
%August 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(useEUN))
    useEUN=false;
end

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(zRot))
    zRot=0;
end

if(nargin<3||isempty(el))
    el=0;
end

if(nargin<2||isempty(az))
    az=0;
end

%Cartesian ECEF to ENU.
u=getENUAxes(plhPoint,false,a,f);

%Cartesian ENU to radar-face XYZ
A = [-cos(az),                  sin(az),          0;
     -sin(az)*sin(el), -cos(az)*sin(el),    cos(el);
      sin(az)*cos(el),	cos(az)*cos(el),    sin(el)];

%Rotate around the local z axis in radar-face XYZ
Mz=Euler1Ang2RotMat(zRot,'z','right');

M=Mz*A*u.';

if(useEUN)
    %If a left-handed East-Up-North coordinate system should be used.
    M(1,:)=-M(1,:); 
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
