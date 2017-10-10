function uv=spherAng2Uv(azEl,systemType,includeW)
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
% includeW An optional boolean value indicating whether a third direction
%          cosine component should be included. The u and v direction
%          cosines are two parts of a 3D unit vector. Generally, one might
%          assume that the target is in front of the sensor, so the third
%          component would be positive and is not needed. However, the
%          third component can be included if ambiguity exists. The default
%          if this parameter is omitted or an empty matrix is passed is
%          false.
%
%OUTPUTS: uv A 2XN (without w) or 3XN (with w) set of direction cosines
%            values corresponding to the specified angles.
%
%Direction cosines and spherical coordinate systems are discussed in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(includeW))
   includeW=false; 
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

azimuth=azEl(1,:);
elevation=azEl(2,:);
numPoints=size(azEl,2);

if(includeW)
    uv=zeros(3,numPoints);
else
    uv=zeros(2,numPoints);
end

if(systemType==2)
    elevation=pi/2-elevation;
    systemType=0;
end

switch(systemType)
    case 0
        uv(1,:)=cos(azimuth).*cos(elevation);
        uv(2,:)=sin(azimuth).*cos(elevation);
        if(includeW)
            uv(3,:)=sin(elevation);
        end
    case 1
        uv(1,:)=sin(azimuth).*cos(elevation);
        uv(2,:)=sin(elevation);
        if(includeW)
            uv(3,:)=cos(azimuth).*cos(elevation);
        end
    otherwise
        error('Invalid system type specified.')
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
