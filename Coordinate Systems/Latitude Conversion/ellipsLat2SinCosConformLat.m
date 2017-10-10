function [sinChi,cosChi]=ellipsLat2SinCosConformLat(ellipsLat,f)
%%ELLIPSLAT2SINCOSCONFORMLAT Given an ellipsoidal (geodetic) latitude of a
%                    point on the surface of a reference ellipsoid, obtain
%                    the sine and cosine of the corresponding conformal
%                    latitude. A conformal latitude is a latitude in terms
%                    of a sphere that is conformal with regard to the
%                    ellipsoid. That is, the transformation from the
%                    ellipsoid to the sphere preserves angles. This is such
%                    that the angle of intersection  between any two lines
%                    on the ellipsoid is the same as the corresponding
%                    angle on the sphere.
%
%INPUTS: ellipsLat A vector or matrix of ellipsoidal (geodetic) latitudes
%                  in radians. These range from -pi/2 to pi/2.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted, the value in
%                  Constants.WGS84Flattening is used.
%
%OUTPUTS: sinChi The sines of the conformal latitudes corresponding to those 
%                in ellipsLat. This matrix has the same dimensions as
%                ellipsLat.
%         cosChi The cosines of the conformal latitudes corresponding to
%                those in ellipsLat. This matrix has the same dimensions as
%                ellipsLat.
%
%Conformal latitudes are defined in Chapter 3 of [1]. The method for
%converting ellipsoidal latitudes directly to the sine and cosine of the
%conformal latitude is takes from Section 2.8 of [2].
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] Office of Geomatics, "National geospatial-intelligence agency
%    standardization document: Implementation practice: The universal
%    grids and the transverse mercator and polar stereographic
%    map projections," 25 Mar. 2014. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

cosPhi=cos(ellipsLat);
sinPhi=sin(ellipsLat);

P=exp(e*atanh(e*sinPhi));

denom=(1+sinPhi)./P+(1-sinPhi).*P;
cosChi=2*cosPhi./denom;
sinChi=((1+sinPhi)./P-(1-sinPhi).*P)./denom;
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
