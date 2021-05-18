function latLon=Mercator2Ellipse(xy,a,f)
%%MERCATOR2ELLIPSE Convert from Mercator Easting and Northing coordinates
%               to ellipsoidal (geodetic) latitude and longitude. The
%               Mercator projection is useful, because lines of constant
%               heading (i.e. degrees East of North) are straight lines on
%               the map and the projection is conformal. However, the poles
%               are placed at +/-Inf. The Mercator projection is a type of
%               cylindrical projection.
%
%INPUTS: xy A 2XN matrix of [Easting;Northing] points in the Mercator
%           projection. Easting is related to longitude and Northing is
%           related to latitude. The units are those of the semi-major axis
%           and are typically meters.
%         a The semi-major axis of the reference ellipsoid. If this
%           argument is omitted, the value in Constants.WGS84SemiMajorAxis
%           is used.
%         f The flattening factor of the reference ellipsoid. If this
%           argument is omitted, the value in Constants.WGS84Flattening is
%           used.
%
%INPUTS: latLon A 2XN matrix of ellipsoidal [latitude;longitude] points
%               given in radians.
%
%Chapter [7] of [1] discusses the Mercator projection, including its
%relation to conformal latitude, which is defined in Chapter 3. The
%solution here uses the relation to conformal latitude and the
%conformLat2EllipsLat function to obtain the ellipsoidal latitude.
%
%Equations for the conversion from ellipsoidal coordinates into Mercator
%coordinates with a commonly used but very inaccurate spherical Mercator
%conversion (the "Web" Mercator map) are given in [2]. The inverse
%spherical Mercator conversion can be obtained by setting f to zero.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] Office of Geomatics, "National geospatial-intelligence agency
%   standardization document: Implementation practice: Web mercator map
%   projection," National Geospatial-Intelligence Agency, Tech. Rep.
%   NGA.SIG.0011 1.0.0 WEBMERC, 18 Feb. 2014. [Online]. Available:
%   http://earth-info.nga.mil/GandG/wgs84/web_mercator/(U)%20NGA_SIG_0011_1.0.0_WEBMERC.pdf
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

x=xy(1,:);
y=xy(2,:);

lambda=wrapRange(x/a,-pi,pi);
conformalLat=Gudermannian(y/a);
phi=conformLat2EllipsLat(conformalLat,f);

latLon=[phi;lambda];

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
