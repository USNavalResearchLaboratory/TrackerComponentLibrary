function M=findRFTransParam2Tar(lRx,lTar,zRot)
%%FINDRFTRANSPARAM2TAR Find the rotation matrix needed to rotate a
%            Cartesian vector to local radar-facing coordinates such
%            that a specified target is on the boresight of the radar. This
%            makes the local z-axis of the radar point at the target.
%
%INPUTS: lRx The 3X1 Cartesian location of the receiver.
%       lTar The 3X1 Cartesian location of the target.
%       zRot This is an additional (counterclockwise) rotation in
%            radians about the pointing direction of the radar (the
%            z-axis). The default if omitted or an empty matrix is passed
%            is zero. If the receiver and transmitter are in the same
%            plane, then zRot=0 means that the local y-axis points in the
%            direction of the global z-axis.
%
%OUTPUTS: M A 3X3 rotation matrix for the transformation from global 
%           Cartesian coordinates to local radar-facing coordinates. The z-
%           axis represents the pointing direction of the radar. This can
%           be directly fed into the RUV coordinate transform functions.
%
%EXAMPLE:
% lRx=[10e3;-15e3;4e3];
% lTar=[20e3;5e3;4e3];
% M=findRFTransParam2Tar(lRx,lTar);
% %One sees that the uv components are zero.
% uv=getUVDirection(lTar,lRx,M)
% %The rotated local y axis points in the direction of the global z axis.
% M'*[0;1;0]
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(zRot))
    zRot=0;
end

xTarSpher=Cart2Sphere(lTar-lRx,0);

az=pi/2-xTarSpher(2);
el=xTarSpher(3);

%Cartesian XYZ to radar-face XYZ
A = [-cos(az),                  sin(az),          0;
     -sin(az)*sin(el), -cos(az)*sin(el),    cos(el);
      sin(az)*cos(el),	cos(az)*cos(el),    sin(el)];

%Rotate around the local z axis in radar-face XYZ
Mz=Euler1Ang2RotMat(zRot,'z','right');

M=Mz*A;
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
