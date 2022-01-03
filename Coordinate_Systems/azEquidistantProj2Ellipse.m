function latLonPts=azEquidistantProj2Ellipse(xyPts,latLonRef,a,f)
%%AZEQUIDISTANTPROJ2ELLIPSE Given points as [x;y] values in a azimuthal
%       equidistant projection, convert the points to latitudes and
%       longitudes on a reference ellipsoid (or sphere).
%
%INPUTS: xyPts A 2XN set of the azimuthal equidistant projection points
%              about latLonRef to convert.
%    latLonRef A 2X1 reference point about which the projection is taken.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonPts The 2XN set of converted [latitude;longitude] points in
%                   radians.
%
%The conversion is described in Chapter 25 of [1], where expressions for
%a spherical Earth are given. However, we do not use those formulae. The
%norm of each xy point corresponds to a distance across the surface of the
%curved Earth and using the inverse tangent function, one can obtain a
%launch azimuth in radians East of North. From there, the latitude and
%logitude of the point can be obtained using the directGreatCircleProb
%function for a spherical Earth or the directGeodeticProb for an
%ellipsoidal Earth.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

numPts=size(xyPts,2);
latLonPts=zeros(2,numPts);
if(f==0)
    %Under a spherical Earth approximation.
    for k=1:numPts
        distVal=norm(xyPts(:,k));
        az=atan2(xyPts(1,k),xyPts(2,k));
        
        latLonPts(:,k)=directGreatCircleProb(latLonRef,az,distVal,a);
    end
else
    %Under an ellipsoidal Earth approximation.
    for k=1:numPts
        distVal=norm(xyPts(:,k));
        az=atan2(xyPts(1,k),xyPts(2,k));
        
        latLonPts(:,k)=directGeodeticProb(latLonRef,az,distVal,a,f);
    end
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
