function angDiff=angBetweenAzElDirs(z1,z2,systemType,angType)
%%ANGBETWEENAZELDIRS Given spherical azimuth and elevation values, find the
%      angle between the directions implied by the values. This function
%      converts the directions into unit vectors and then calls
%      angBetweenVecs.
%
%INPUTS: z1, z2 The 1XN set of [azimuth;elevation] in radians (unless
%            angType=1)of the first and second vectors in each of N pairs.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are:
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
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
%   angType An optional parameter specifying the units of the angles in z1
%           and z2 and the units desired for the output. Possible values
%           are:
%           0 (The default if omitted or an empty matrix is passed) The
%             values are given in radians.
%           1 The values are given in degrees.
%
%OUTPUTS: angDiff An 1XN matrix of the angular differences (in radians
%                 unless angType=1) between the vectors in v1 and the
%                 corresponding vectors in v2. The distances can range from
%                 0 to pi radians (or 0 to 180 in degrees).
%
%EXAMPLE:
%This is an example where the directions represented by the spherical
%coordinates are two degrees apart.
% z1=[45;89];
% z2=[-135;89];
% systemType=0;
% angType=1;%Degrees
% angDiff=angBetweenAzElDirs(z1,z2,systemType,angType)
%
%March 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(angType))
    angType=0;
end

if(nargin<3||isempty(systemType))
    systemType=0;
end

if(angType==1)
    angDiff=rad2deg(angBetweenVecs(spher2Cart(deg2rad(z1),systemType),spher2Cart(deg2rad(z2),systemType)));
else
    angDiff=angBetweenVecs(spher2Cart(z1,systemType),spher2Cart(z2,systemType));
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
