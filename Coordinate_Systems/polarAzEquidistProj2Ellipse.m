function latLonPts=polarAzEquidistProj2Ellipse(pAzPts,latLonRef,a,f)
%%POLARAZEQUIDISTANTPROJ2ELLIPSE Given points as polar coordinates of an
%       azimuthal equidistant projection, convert the points to latitudes
%       and longitudes on a reference ellipsoid (or sphere). A third height
%       coordinate can also be provided, which doesn't change in the
%       conversion.
%
%INPUTS: pAzPts A 2XN set of the [ground distance; heading] points, with the
%              heading given in radians East of North, to convert.
%              Alternatively, if heights are also given, this can be a 3XN
%              set of points with the height being the third dimension.
%    latLonRef A 2X1 [latitude;longitude] reference point in radians about
%              which the projection is taken.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonPts The 2XN set of converted [latitude;longitude] points in
%                   radians. If heights were given in xyPts, then this is a
%                   3XN set of converted [latitude;longitude;height] points
%                   with the third row the same as in xyPts.
%
%The conversion is described in Chapter 25 of [1], where expressions for
%a spherical Earth are given. However, we do not use those formulae. The
%latitude and longitude of the point can be obtained using the
%directGreatCircleProb function for a spherical Earth or the
%directGeodeticProb for an ellipsoidal Earth, which is just what is done
%here.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

hasHeight=(size(pAzPts,1)==3);
numPts=size(pAzPts,2);
if(hasHeight)
    latLonPts=zeros(3,numPts);
else
    latLonPts=zeros(2,numPts);
end
if(f==0)
    %Under a spherical Earth approximation.
    for k=1:numPts
        distVal=pAzPts(1,k);
        az=pAzPts(2,k);
        
        latLonPts(1:2,k)=directGreatCircleProb(latLonRef,az,distVal,a);
    end
else
    %Under an ellipsoidal Earth approximation.
    for k=1:numPts
        distVal=pAzPts(1,k);
        az=pAzPts(2,k);
        
        latLonPts(1:2,k)=directGeodeticProb(latLonRef,az,distVal,a,f);
    end
end

if(hasHeight)
    latLonPts(3,:)=pAzPts(3,:);
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
