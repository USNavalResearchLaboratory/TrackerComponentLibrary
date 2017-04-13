function points=spher2Ellipse(points,a,f)
%%SPHER2CART  Convert points from spherical coordinates to ellipsoidal
%             coordinates. The spherical coordinates can either be full
%             coordinates [range;azimuth; elevation] or just pairs of 
%             [azimuth; elevation] or just [elevation] in which case it is
%             assumed that the points are on the surface of the reference
%             ellipsoid and one just has to convert the elevation
%             (spherical latitude) from geocentric to geodetic.
%
%INPUTS:  points  One or more points given in terms of range, azimuth and
%                 elevation, with the angles in radians. To convert
%                 N points, points is a 3XN matrix with each column
%                 having the format [range;azimuth; elevation].
%                 Alternatively, the points can be given just as [azimuth;
%                 elevation] pairs or just as elevations (geocentric
%                 latitudes). In such an instance, the elevations are
%                 converted to geodetic (ellipsoidal) latitudes for points
%                 located on the surface of the reference ellipsoid.
%           a     The semi-major axis of the reference ellipsoid. If
%                 this argument is omitted, the value in
%                 Constants.WGS84SemiMajorAxis is used.
%           f     The flattening factor of the reference ellipsoid. If
%                 this argument is omitted, the value in
%                 Constants.WGS84Flattening is used.
%
%OUTPUTS:   points  A matrix of the converted points. Each column of the
%                   matrix has the format [latitude;longitude;altitude] or
%                   [latitude;longitude] ot just latitude, respectively
%                   depending on the dimensionality of the input. The
%                   angles are given in radians.
%
%Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
%is the angle above the x-y plane. The formula for directly converting
%between geocentric and geodetic latitudes is described in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An overview of major terrestrial, celestial, and
%    temporal coordinate systems for target tracking", Report, U. S. Naval
%    Research Laboratory, to appear, 2016.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<2||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

pDim=size(points,1);
%If geocentric latitude and longitude (implied on the surface of the
%reference ellipsoid) are given and the latitude must be converted to
%geodetic latitude.
if(pDim<3)
    e2=2*f-f^2;%Get the square of the eccentricity of the reference ellipse.
    if(pDim==2)
        geocenLat=points(2,:);
    else%Assume that the points are 1D
        geocenLat=points;
    end
    
    geodetLat=asin(sign(geocenLat).*sqrt(2*sin(geocenLat).^2./(2-2*e2+e2^2+e2*(e2-2)*cos(2*geocenLat))));
    
    if(pDim==2)
        points(2,:)=points(1,:);
        points(1,:)=geodetLat;
    else
        points=geodetLat;
    end
else
    points=Cart2Ellipse(spher2Cart(points),[],a,f);
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
