function uv=spherAng2Uv(azEl,systemType,includeW,Ms,Muv)
%%SPHERANG2UV Convert azimuth and elevation in spherical coordinates into a
%             [u;v] or [u;v;w] direction cosine vector in 3D. Direction
%             cosines u and v are the first two elements of a unit vector
%             in 3D.
%
%INPUTS: azEl A 2XN set of N points in [azimuth;elevation] in radians to
%           convert to direction cosines.
%  systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. 
%             Elevation is measured up from the x-y plane (towards the
%             z-axis). This is consistent with common spherical
%             coordinate systems for specifying longitude (azimuth) and
%             geocentric latitude (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%  includeW An optional boolean value indicating whether a third direction
%           cosine component should be included. The u and v direction
%           cosines are two parts of a 3D unit vector. Generally, one might
%           assume that the target is in front of the sensor, so the third
%           component would be positive and is not needed. However, the
%           third component can be included if ambiguity exists. The
%           default if this parameter is omitted or an empty matrix is
%           passed is false.
%    Ms,Muv If either the spherical coordinate system or the u-v coordinate
%           system is rotated compared to the global Cartesian coordinate
%           system, these optional 3X3 matrices provide the rotations. Ms
%           is a 3X3 matrix to go from the alignment of a global
%           Cartesian coordinate system to that in which the spherical
%           coordinates are computed. Similarly, Muv is a rotation matrix
%           to go from the alignment of a global Cartesian cordinate system
%           to that in which the u-v(-w) coordinates are computed. If
%           either of these in omitted or an empty matrix is passed, then
%           the missing one is replaced with the identity matrix.
%
%OUTPUTS: uv A 2XN (without w) or 3XN (with w) set of direction cosines
%            values corresponding to the specified angles.
%
%Direction cosines and spherical coordinate systems are discussed in [1].
%
%EXAMPLE:
%In this example, a spherical value is converted to uv coordiantes and then
%back, demonostrating the consistency of the functions spherAng2Uv and
%uv2SpherAng. The relative error should be about zero.
% zSpher=[0.4;0.7];
% systemType=3;
% includeW=true;
% Ms=Euler2Ang2RotMat(0.1,0.2,'xy');
% Muv=Euler1Ang2RotMat(0.25,'z');
% convBack=uv2SpherAng(spherAng2Uv(zSpher,systemType,includeW,Ms,Muv),systemType,Ms,Muv);
% relErr=max(abs((convBack-zSpher)./zSpher))
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(Muv))
    Muv=eye(3,3);
end

if(nargin<4||isempty(Ms))
    Ms=eye(3,3);
end

if(nargin<3||isempty(includeW))
    includeW=false; 
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

azimuth=azEl(1,:);
elevation=azEl(2,:);

if(systemType==2)
    elevation=pi/2-elevation;
    systemType=0;
elseif(systemType==3)
    azimuth=pi/2-azimuth;
    systemType=0;
end

switch(systemType)
    case 0
        uv=Muv*Ms'*[cos(azimuth).*cos(elevation);
                    sin(azimuth).*cos(elevation);
                    sin(elevation)];
    case 1
        uv=Muv*Ms'*[sin(azimuth).*cos(elevation);
                    sin(elevation);
                    cos(azimuth).*cos(elevation)];
    otherwise
        error('Invalid system type specified.')
end

if(includeW==false)
    uv=uv(1:2,:);
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
