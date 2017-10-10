function phi=conformLat2EllipsLat(conformalLat,f)
%%CONFORMLAT2ELLIPSELAT Given a conformal latitude, determine an
%           ellipsoidal (geodetic) latitude. A conformal latitude is a
%           latitude in terms of a sphere that is conformal with regard to
%           the ellipsoid. That is, the transformation from the ellipsoid
%           to the sphere preserves angles. This is such that the angle of
%           intersection  between any two lines on the ellipsoid is the
%           same as the corresponding angle on the sphere.
%
%INPUTS: conformalLat The conformal latitudes in radians to convert. This
%               can be a scalar or a matrix.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in Constants.WGS84Flattening
%               is used. This function is only accurate for values from 0
%               to about 0.999. Typically, this will be close to zero, as
%               in Constants.WGS84Flattening.
%
%OUTPUTS: phi The ellipsoidal latitudes in radians. This is the same size as
%             conformallat.
%
%Conformal latitudes are discussed in Chapter 3 of [1]. This function
%implements calls sinCosConformLat2EllipsLat to solve the problem via a
%fixed point iteration. 
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

phi=sinCosConformLat2EllipsLat(sin(conformalLat),cos(conformalLat),f);

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
