function ellipsLat=reducedLat2EllipsLat(reducedLat,f)
%%REDUCEDLAT2ELLIPSLAT Convert reduced (parameteric) latitude into an
%                      ellipsoidal (geodetic) latitude. A parameteric
%                      latitude is an angle of elevation drawn from the
%                      center of a sphere whose radius is equal to the
%                      semi-major axis of the ellipsoid. The parametric
%                      latitude is the angle of the point on the sphere
%                      where a vertical line drawn from the major axis of
%                      the ellipse through the point on the ellipse
%                      intersects the sphere. This is a way of mapping an
%                      ellipsoid onto a sphere.
%
%INPUTS: reducedLat A vector or matrix of reduced latitudes in radians.
%                   These typically range from -pi/2 to pi/2.
%                 f The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS: reducedLat The ellipsoidal (geodetic) latitude ranging from
%                   -pi/2 to pi/2, for each reduced latitude.
%
%%The definition of the reduced (parameteric) latitude from Chapter 3 (pg.
%18) of [1] is used. The conversion relation between the latitudes is in
%Section 3.4 of [2].
%
%The reduced latitude has been used for solving geodescic problems on an
%ellipsoid by transforming the problem to an equivalent problem on a
%sphere.
%
%The inverse of this function is ellipsLat2ReducedLat.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] R. H. Rapp, "Geometric geodesy, part I," Ohio State University
%    Department of Geodetic Science and Surveying, Tech. Rep., Apr. 1991.
%    [Online]. Available: http://hdl.handle.net/1811/24333
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

%The wrapping helps deal with precision limitations that might push the
%input value slightly above or below +/-pi/2. If the wrapping were not
%done, then instead of getting, for example, pi/2, one will get pi/2 as a
%return value.
reducedLat=wrapRange(reducedLat,-pi/2,pi/2,true);

ellipsLat=atan(tan(reducedLat)/sqrt(1-e^2));
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
