function isoLat=ellipsLat2IsoLat(ellipsLat,f)
%ELLIPSLAT2ISOLAT Convert an ellipsoidal (geodetic) latitude to an
%            isometric latitude. The isometric latitude is used in Mercator
%            and transverse Mercator map projections. The isometry comes
%            from the fact that at any point, equal displacements along
%            lines of constant isometric latitude and longitude lead to
%            equal distance displacements along the meridians and
%            parallels. In other words a map in isometric latitude and
%            longitude preserves angles while warping distances. The
%            isometric latitude is sometimes referred to as Psi.
%
%INPUTS: ellipsLat A vector or matrix of ellipsoidal (geodetic) latitudes
%                  in radians. These range from -pi/2 to pi/2.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted, the value in
%                  Constants.WGS84Flattening is used.
%
%OUTPUTS: isoLat The isometric latitudes corresponding to the ellipsoidal
%                latitudes in ellipsLat. This ranges from -Inf to Inf.
%
%The isometric latitude is the integral from 0 to ellipsLat of
%(1-e^2)*sec(L)/(1-e^2*sin(L)^2) dL, where e is the first numerical
%eccentricity of the ellipse. It has an explicit solution, which is
%used below. The explicit solution is taken from Chapter 3 (Equation 3-7)
%of [1] and has been transformed to a slightly simpler form.
%
%The inverse of this function is isoLat2EllipsLat.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

s=sin(ellipsLat);
isoLat=atanh(s)-e*atanh(e*s);
%Note that atanh(sin(ellipsLat))=invGudermannian(ellipsLat).
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
