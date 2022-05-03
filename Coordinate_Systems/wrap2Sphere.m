function wrapVals=wrap2Sphere(azEl,systemType)
%%WRAP2SPHERE Given azimuth-elevation directions in spherical coordinates
%             that might be outside of the standard range of -pi to pi for
%             azimuth and -pi/2 to pi/2 for evevation, wrap the points to
%             the sphere. This function is useful when offsets (like
%             simulated noise) are added to spherical coordinates or if
%             one is performing Newton's method in spherical coordinates
%             and is likely to pass a pole. For example, if using
%             systemType=0 and elevation exceeds pi/2 by a small epsilon,
%             then it must be reduced to pi/2-epsilon and azimuth offset by
%             pi.
%
%INPUTS: azEl A 2XN set of N spherical coordinates in radians that should
%             be wrapped to the normal sphere. Each point is of the format
%             [azimuth;elevation] and the angles are measured according to
%             systemType.
%  systemType An optional parameter specifying the axis from which the
%             angles are measured in radians. Possible values are
%             0 (The default if omitted) Azimuth is measured 
%               counterclockwise from the x-axis in the x-y plane.
%               Elevation is measured up from the x-y plane (towards the
%               z-axis). This is consistent with common spherical
%               coordinate systems for specifying longitude (azimuth) and
%               geocentric latitude (elevation).
%             1 Azimuth is measured counterclockwise from the z-axis in the
%               z-x plane. Elevation is measured up from the z-x plane
%               (towards the y-axis). This is consistent with some
%               spherical coordinate systems that use the z-axis as the
%               boresight direction of the radar.
%             2 This is the same as 0 except instead of being given
%               elevation, one is given the angle away from the z-axis,
%               which is (pi/2-elevation).
%             3 This is the same as 0 except azimuth is measured clockwise
%               from the y-axis in the x-y plane instead of
%               counterclockwise from the x-axis. This coordinate system
%               often arises when given "bearings" in a local East-North-Up
%               coordinate system, where the bearing directions are
%               measured East of North.
%
%OUTPUTS: wrapVals Values of azEl that have been wrapped to the sphere in
%                  radians. azEl(1,:) ranges from -pi to pi (pi not
%                  included) and azEl(2,:) is from -pi/2 to pi/2.
%
%This function works by first converting the direction into a unit vector
%and then by converting the unit vector back into azimuth and elevation.
%The use of trigonometric functions for the first step assures the smooth
%aliasing of points over the poles and the use of inverse trigonometric
%functions for the second step places the results in the correct standard
%region.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0;
end

azimuth=azEl(1,:);
elevation=azEl(2,:);

%Get unit vectors pointing in the direction of the target.
switch(systemType)
    case 0
        x=cos(azimuth).*cos(elevation);
        y=sin(azimuth).*cos(elevation);
        z=sin(elevation);
        
        azimuth=atan2(y,x);
        elevation=asin(z);
    case 1
        x=sin(azimuth).*cos(elevation);
        y=sin(elevation);
        z=cos(azimuth).*cos(elevation);
        
        azimuth=atan2(x,z);
        elevation=asin(y);
    case 2
        x=cos(azimuth).*sin(elevation);
        y=sin(azimuth).*sin(elevation);
        z=cos(elevation);

        azimuth=atan2(y,x);
        elevation=pi/2-asin(z);
    case 3
        x=cos(pi/2-azimuth).*cos(elevation);
        y=sin(pi/2-azimuth).*cos(elevation);
        z=sin(elevation);
        
        azimuth=pi/2-atan2(y,x);
        elevation=asin(z);
    otherwise
        error('Invalid system type specified.')
end
wrapVals=[azimuth;elevation];
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
