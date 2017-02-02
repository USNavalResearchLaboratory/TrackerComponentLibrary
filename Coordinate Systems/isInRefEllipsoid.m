function val=isInRefEllipsoid(z,a,f)
%%ISINREFELLIPSOID Determine whether a point is inside of a reference
%                   ellipsoid centered at the origin and oriented with the
%                   z axis as the axis of rotation. This function is useful
%                   for approximating whether a target estimate is
%                   underground when terrain data are not available. The
%                   WGS-84 reference ellipsoid is a logical choice.
%
%INPUTS:    z       A point in [x;y;z] Cartesian coordinates.
%           a       The semi-major axis of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%           f       The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS: val       True if the point z is within the ellipsoid and false
%                   if the point is not within the reference ellipsoid.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

    b=a*(f-1);
    val=(z(1)^2+z(2)^2)/a^2+z(3)^2/b<1;
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
