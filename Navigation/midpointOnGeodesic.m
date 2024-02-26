function latLonHalf=midpointOnGeodesic(latLonStart,latLonEnd,a,f)
%%MIDPOINTONGEODESIC Given two points on the surface of the
%    reference ellipsoid or sphere (f=0), find the points that is halfway
%    along the geodesic between the points. This can be useful for finding
%    a midway point to use as a reference point. For example, if one wishes
%    to use an osculating sphere that fits decently to both points on the
%    reference ellipsoid, choosing the midpoints between them as the
%    reference is one option.
%
%INPUTS: latLonStart The [latitude;longitude] in radians of the first point
%            on the geodesic curve.
%    latLonEnd The [latitude;longitude] in radians of the second point on
%            the geodesic curve. The curve defined goes from latLonStart to
%            latLonEnd. The shortest way around the world is chosen.
%          a The semi-major axis of the reference ellipsoid (in meters). If
%            this argument is omitted or an empty matrix is passed, the
%            value in Constants.WGS84SemiMajorAxis is used.
%          f The flattening factor of the reference ellipsoid. If this
%            argument is omitted or an empty matrix is passed, the value in
%            Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonHalf The 2X1 [latitude;longitude] point in radians midway
%                    between latLonStart and latLonEnd
%
%Note that the ordering of latLonStart and latLonEnd does not matter.
%
%If f=0, then a reference sphere is used and this function calls
%indirectGreatCircleProb and then calls directGreatCircleProb with half the
%distance returned by indirectGreatCircleProb. If f~=0, then the same thing
%is done using indirectGeodeticProb and directGeodeticProb.
%
%EXAMPLE:
%Two points are shown on the Earth along with a geodesic between them and
%the midpoint is marked with an x.
% latLonStart=deg2rad([19.7216; -155.0849]);
% latLonEnd=deg2rad([34.0736;-118.4004]);
% [azStart,distVal]=indirectGeodeticProb(latLonStart,latLonEnd);
% numPts=300;
% distVals=linspace(0,distVal,numPts);
% %A height above the Earth's surface is used to make the line better appear.
% wayPts=ellips2Cart([directGeodeticProb(latLonStart,azStart,distVals);100e3*ones(1,numPts)]);
% startPt=ellips2Cart(latLonStart);
% halfPt=ellips2Cart([midpointOnGeodesic(latLonStart,latLonEnd);100e3]);
% endPt=ellips2Cart(latLonEnd);
% 
% figure(1)
% clf
% plotMapOnEllipsoid()
% hold on
% plot3(wayPts(1,:),wayPts(2,:),wayPts(3,:),'-r','linewidth',4)
% scatter3(startPt(1),startPt(2),startPt(3),500,'.m')
% scatter3(halfPt(1),halfPt(2),halfPt(3),500,'xc','linewidth',4)
% scatter3(endPt(1),endPt(2),endPt(3),500,'.m')
%
%October 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(f==0)
    %Reference sphere.
     [azStart,distVal]=indirectGreatCircleProb(latLonStart,latLonEnd,a);
     latLonHalf=directGreatCircleProb(latLonStart,azStart,distVal/2,a);
else
    %Reference ellipsoid.
    [azStart,distVal]=indirectGeodeticProb(latLonStart,latLonEnd,a,f);
    latLonHalf=directGeodeticProb(latLonStart,azStart,distVal/2,a,f);
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
