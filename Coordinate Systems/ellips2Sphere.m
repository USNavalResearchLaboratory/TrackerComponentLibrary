function points=ellips2Sphere(points,a,f)
%%ELLIPS2SPHERE Convert ellipsoidal coordinates to spherical coordinates.
%               If the ellipsoidal coordinates are not full coordinates,
%               then it is assumed that one is just converting latitudes on
%               the surface of the reference ellipsoid.
%
%INPUTS:    points  One or more points given in geodetic latitude and
%                   longitude, in radians, and height, in meters that are
%                   to be converted to Cartesian coordinates. To convert
%                   N points, points is a 3XN matrix with each column
%                   having the format [latitude;longitude; height].
%                   Alternately, one can just provide points as
%                   [latitude;longitude], where height is implied as zero,
%                   and the return values are 2D geocentric (latitude
%                   longitude pairs). Also, one can just provide points as
%                   latitudes, in which case spherical latitudes for points
%                   on the surface of the reference ellipsoid are returned.
%           a       The semi-major axis of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%           f       The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS:   points  A matrix of the converted points. Each column of the
%                   matrix has the format [r;azimuth;elevation],
%                   with azimuth and elevation given in radians, unless the
%                   input points were given as 2D [latitude;longitude]
%                   points, in which case the output points are of the form
%                   [azimuth;elevation], unless the input points were just
%                   given as latitudes, in which case the output points are
%                   just elevations (geocentric latitudes).
%
%Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
%is the angle above the x-y plane.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

pDim=size(points,1);
if(pDim<3)%If only latitude and longitudes or just latitudes are provided.
    geodetLat=points(1,:);
    
    geocenLat=ellipsLat2SpherLat(geodetLat,f);
        
    if(pDim==2)
        points(1,:)=points(2,:);
        points(2,:)=geocenLat;
    else
        points=geocenLat;
    end
else
    points=Cart2Sphere(ellips2Cart(points,a,f));
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
