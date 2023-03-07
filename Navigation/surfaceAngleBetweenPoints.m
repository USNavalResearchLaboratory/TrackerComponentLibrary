function angle=surfaceAngleBetweenPoints(latLon1,latLonVertex,latLon2,useSmallestPosSol,a,f)
%%SURFACEANGLEBETWEENPOINTS Given an angle formed by the three points 
%               latLon1->latLonVertex->latLon2 in that order, determine the
%               angle between the geodesics latLonVertex->latLon1 and
%               latLonVertex->latLon2 across the ground as measured at
%               latLonVertex. If angle1 is the angle in radians East of
%               North in the direction of latLon1 as measured in the local
%               tangent plane at latLonVertex and angle2 is the same in the
%               direction of latLon2, then if useSmallestPosSol is false,
%               then angle2=angle1+angle. On the other hand, if
%               useSmallestPosSol is true, then angle is the positive angle
%               between the (unordered) latLonVertex->latLon1 and
%               latLonVertex->latLon2 rays going the shortest distance
%               around the circle (clockwise or counterclockwise).
%
%INPUTS: latLon1 The 2X1 [latitude;longitude] of the first endpoint in
%                radians.
%   latLonVertex The 2X1 [latitude;longitude] of the vertex in radians.
%        latLon2 The 2X1 [latitude;longitude] of the second endpoint in
%                radians.
% useSmallestPosSol If this input is false then angle is what would have to
%                be added to the angle from latLonVertex in radians East of
%                North to get the geodesic heading in radians East of North
%                to go to latLon2. Otherwise, if this input is true, then
%                this is the postivie angle between those two legs, the
%                shortest way around the circle. The default if omitted or
%                an empty matrix is passed is true.
%              a The semi-major axis of the reference ellipsoid (in
%                meters). If this argument is omitted or an empty matrix is
%                passed, the value in Constants.WGS84SemiMajorAxis is used.
%                If f=0, then this input is not used and it can be omitted
%                or an empty matrix can be passed.
%              f The flattening factor of the reference ellipsoid. If this
%                argument is omitted or an empty matrix is passed, the
%                value in Constants.WGS84Flattening is used.
%
%OUTPUTS: angle The latLon1->latLonVertex->latLon2 angle as specified by
%               useSmallestSol.
%
%This function just computed the angle to each point from latLonVertex
%using indirectGeodeticProb (if f is not 0) or greatCircleAzimuth (if f=0)
%and takes the differences, wrapping the result to +/-pi and taking the
%absolute value if useSmallestPosSol is true.
%
%EXAMPLE:
%This just shows that the function outputs the expected relations between
%the angles. That is, the output is the same regardingess of the ordering
%of latLon1 and latLon2 if useSmallestPosSol is true and if
%useSmallestPosSol is false, the angular relation between angle1 angle2 and
%angle is verified.
% latLon1=deg2rad([38.563737;-121.434812]);%Near Sacramento.
% latLonVertex=deg2rad([34.105551;-118.136759]);%Near LA.
% latLon2=deg2rad([33.337986;-111.980393]);%Near Phoenix.
% useSmallestPosSol=true;
% vertexAng1=surfaceAngleBetweenPoints(latLon1,latLonVertex,latLon2,useSmallestPosSol);
% vertexAng2=surfaceAngleBetweenPoints(latLon2,latLonVertex,latLon1,useSmallestPosSol);
% %Check that the ordering doesn't matter when useSmallestPosSol is true.
% vertexAng1==vertexAng2
% 
% %Now, compute it ordered.
% useSmallestPosSol=false;
% vertexAng1=surfaceAngleBetweenPoints(latLon1,latLonVertex,latLon2,useSmallestPosSol);
% vertexAng2=surfaceAngleBetweenPoints(latLon2,latLonVertex,latLon1,useSmallestPosSol);
% 
% angle1=indirectGeodeticProb(latLonVertex,latLon1);
% angle2=indirectGeodeticProb(latLonVertex,latLon2);
% %Get to the other angle with the first angle and vertexAng1 or vertexAng2.
% %The following values should be zero within finite precision limits.
% (angle1+vertexAng1)-angle2
% (angle2+vertexAng2)-angle1
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(useSmallestPosSol))
    useSmallestPosSol=true;
end

if(f~=0)
    alpha1=indirectGeodeticProb(latLonVertex,latLon1,a,f);
    alpha2=indirectGeodeticProb(latLonVertex,latLon2,a,f);
else
    alpha1=greatCircleAzimuth(latLonVertex,latLon1);
    alpha2=greatCircleAzimuth(latLonVertex,latLon2);
end

angle=wrapRange(alpha2-alpha1,-pi,pi);

if(useSmallestPosSol)
    angle=abs(angle);
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
