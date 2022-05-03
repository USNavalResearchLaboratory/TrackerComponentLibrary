function zMean=meanDirectionAzEl(zAzEl,systemType,w)
%%MEANDIRECTIONAZEL Given multiple direction specified by spherical
%       coordinates (azimuths and elevations in radians), find the mean
%       direction from the values. Optionally, weights can be provided for
%       a weighed mean.
%
%INPUTS: zAzEl A 2XnumSample Array of [azimuth;elevation] pairs
% systemType An optional parameter specifying the axes from which
%            the angles are measured. Possible values are:
%            0 (The default if omitted) Azimuth is measured
%              counterclockwise from the x-axis in the x-y plane.
%              Elevation is measured up from the x-y plane (towards the
%              z-axis). This is consistent with common spherical
%              coordinate systems for specifying longitude (azimuth)
%              and geocentric latitude (elevation).
%            1 Azimuth is measured counterclockwise from the z-axis in
%              the z-x plane. Elevation is measured up from the z-x plane
%              (towards the y-axis). This is consistent with some
%              spherical coordinate systems that use the z axis as
%              the boresight direction of the radar.
%            2 This is the same as 0 except instead of being given
%              elevation, one desires the angle away from the z-axis,
%              which is (pi/2-elevation).
%            3 This is the same as 0 except azimuth is measured clockwise
%              from the y-axis in the x-y plane instead of counterclockwise
%              from the x-axis. This coordinate system often arises when
%              given "bearings" in a local East-North-Up coordinate system,
%              where the bearing directions are measured East of North.
%          w An optional vector of weights having the same dimensionality
%            as ang, if a weighted average is desired. The weights do not
%            have to sum to one. If omitted or an empty matrix is passed,
%            all measurements are treated as having the same weight.
%
%OUTPUT: zMean A 2X1 fused measurement in the format of
%              [azimuth;elevation].
%
%Due to the spherical nature of the measurements, one cannot simply average
%the values. Thus, they are first converted to Cartesian unit vectors, the
%unit vectors are averaged and then the result it converted back to
%spherical azimuth and elevation angles.
%
%EXAMPLE:
%We have two points that are 180 degrees apart in azimuth and the same
%positive elevation. The mean of those point is the vertical (90 degrees
%elevation and azny azimuthal value), which is what is obtained.
% azEl=zeros(2,2);
% azEl(:,1)=[50;89]*(pi/180);
% azEl(:,2)=[-130;89]*(pi/180);
% meanDirectionAzEl(azEl)*(180/pi)
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    w=[];
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Convert to Cartesian unit vectors.
uVecs=spher2Cart(zAzEl,systemType);

if(~isempty(w))
    meanVec=sum(bsxfun(@times,uVecs,w(:).'),2);
else
    numVec=size(uVecs,2);
    meanVec=sum(uVecs,2)/numVec;
end

zMeanRAzEl=Cart2Sphere(meanVec,systemType);
zMean=zMeanRAzEl(2:3);

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
