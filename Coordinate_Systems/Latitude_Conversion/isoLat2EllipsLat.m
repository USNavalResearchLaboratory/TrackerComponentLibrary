function ellipsLat=isoLat2EllipsLat(isoLat,f)
%ISOLAT2ELLIPSLAT Convert an isometric latitude to an ellipsoidal
%                 (geodetic) latitude.The isometric latitude is used in
%                 Mercator and transverse Mercator map projections. The
%                 isometry comes from the fact that at any point, equal
%                 displacements along lines of constant isometric latitude
%                 and longitude lead to equal distance displacements along
%                 the meridians and parallels. In other words a map in
%                 isometric latitude and longitude preserves angles while
%                 warping distances. The isometric latitude is sometimes
%                 referred to as Psi.
%
%INPUTS: isoLat A vector or matrix of real isometric latitudes. These range
%               in value from -Inf to Inf.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in Constants.WGS84Flattening
%               is used.
%
%OUTPUTS: ellipsLat The ellipsoidal latitudes in radians corresponding to
%                   the isometric latitudes in isoLat. These range from
%                   -pi/2 to pi/2.
%
%The isometric latitude is the integral from 0 to ellipsLat of
%(1-e^2)*sec(L)/(1-e^2*sin(L)^2) dL, where e is the first numerical
%eccentricity of the ellipse. It has an explicit solution, given in [1].
%The inverse iteration is taken from the same chapter.
%
%The inverse of this function is ellipsLat2IsoLat.
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

%The initial estimate
ellipsLat=2*atan(exp(isoLat))-pi/2;

%Convergence seems to occur in 7 or fewer iterations when using
%double-precision arithmetic.
maxIter=7;
for curIter=1:maxIter
    s=sin(ellipsLat);
    ellipsLat=2*atan(exp(isoLat)*((1+e*s)/(1-e*s))^(e/2))-pi/2;
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
