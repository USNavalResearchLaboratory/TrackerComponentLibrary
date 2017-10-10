function ellipsLat=spherLat2EllipsLat(spherLat,f)
%%SPHERLAT2ELLIPSLAT  Convert a spherical latitude (geocentric latitude/an
%                     elevation angle) measured up from the x-y plane into
%                     an ellipsoidal (geodetic) latitude, when considering
%                     a point on the surface of the reference ellipsoid.
%
%INPUTS: spherLat A vector or matrix of spherical (geocentric) latitudes
%                 in radians. These range from -pi/2 to pi/2, and are
%                 measured up from the x-y plane.
%               f The flattening factor of the reference ellipsoid. If
%                 this argument is omitted, the value in
%                 Constants.WGS84Flattening is used.
%
%OUTPUTS: ellipsLat The ellipsoidal (geodetic) latitudes corresponding to
%                   the spherical latitudes in spherLat. This ranges from
%                   -pi/2 to pi/2.
%
%The formula for converting between spherical and ellipsoidal latitudes is
%in Section 3.4 of [1].
%
%REFERENCES:
%[1] R. H. Rapp, "Geometric geodesy, part I," Ohio State University
%    Department of Geodetic Science and Surveying, Tech. Rep., Apr. 1991.
%    [Online]. Available: http://hdl.handle.net/1811/24333
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

%The wrapping helps deal with precision limitations that might push the
%input value slightly above or below +/-pi/2. If the wrapping were not
%done, then instead of getting, for example, -pi/2, one will get pi/2 as a
%return value.
spherLat=wrapRange(spherLat,-pi/2,pi/2,true);

ellipsLat=atan((1/(1-e^2))*tan(spherLat));
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
