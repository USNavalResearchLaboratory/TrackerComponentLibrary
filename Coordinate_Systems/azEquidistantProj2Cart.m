function xCart=azEquidistantProj2Cart(azEqPts,latLonRef,a,f)
%%AXEQUIDISTANCEPROJ2CART Given points as [x;y] values in a azimuthal
%       equidistant projection, convert the points to 3D Cartesian
%       coordinates on a reference ellipsoid (or sphere). A third height
%       coordinate can also be provided (and is assumed 0 if omitted).
%
%INPUTS: azEqPts A 2XN set of the azimuthal equidistant projection points
%              about latLonRef to convert. Alternatively, if heights are
%              given, this can be a 3XN set of points.
%    latLonRef A 2X1 [latitude;longitude] reference point in radians about
%              which the projection is taken.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%
%OUTPUTS: xCart A 3XN set of N 3D Cartesian points converted from azEqPts.
%
%This function just calls azEquidistantProj2Ellipse and ellips2Cart.
%
%January 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

latLonPts=azEquidistantProj2Ellipse(azEqPts,latLonRef,a,f);
xCart=ellips2Cart(latLonPts,a,f);
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
