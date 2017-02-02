function u=geogHeading2uVec(point,geoEastOfNorth,a,f)
%%GEOHEADING2UVEC Convert geographic headings in radians East of true North
%                 at a point on a reference ellipsoid to unit vectors in
%                 ECEF coordinates.
%
%INPUTS:  point  The location of the point in geodetic latitude and
%                longitude in radians at which the headings are taken.
%                Point can be [latitude;longitude] or
%                [latitude;longitude;height]. The height component is
%                ignored if included since it does not change the result.
%geoEastOfNorth  An NX1 or 1XN array of N geographic headings in radians
%                clockwise from North that should be turned into ECEF
%                unit vectors. A geographic heading is a direction in the
%                local tangent plane of an East-North-Up coordinate system
%                as defined on a specific reference ellipsoid.
%           a    The semi-major axis of the reference ellipsoid. If this
%                argument is omitted, the value in
%                Constants.WGS84SemiMajorAxis is used.
%           f    The flattening factor of the reference ellipsoid. If this
%                argument is omitted, the value in
%                Constants.WGS84Flattening is used.
%
%OUTPUTS:  u     A 3XN matrix of unit vectors pointing in the directions
%                implies by the headings in geoEastOfNorth.
%
%A heading can be turned into components of a local East-North-Up
%coordinate system. sin(geoEastOfNorth) is the component East and
%cos(geoEastOfNorth) is the component North. Thus, a unit vector in the
%appropriate direction in ECEF coordinates can be obtained by multiplying
%the components by the unit vectors in the appropriate directions. The East
%and North unit vectors are obtained using the getENUAxes function.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    f=Constants.WGS84Flattening;
end

if(nargin<3)
    a=Constants.WGS84SemiMajorAxis;
end

%Transpose geoEastOfNorth to make it 1XN in size if it is NX1. This is so
%that bsxfun can be used to multiply the components with the basis vectors
%instead of using a loop.
if(size(geoEastOfNorth,2)==1)
    geoEastOfNorth=geoEastOfNorth';
end

%Get the the unit vectors for East, North, and Up.
uLocal=getENUAxes(point,false,a,f);
uEast=uLocal(:,1);
uNorth=uLocal(:,2);

u=bsxfun(@times,sin(geoEastOfNorth),uEast)+bsxfun(@times,cos(geoEastOfNorth),uNorth);

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
