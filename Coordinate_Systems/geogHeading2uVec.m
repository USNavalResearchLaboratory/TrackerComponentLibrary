function u=geogHeading2uVec(point,geoEastOfNorth,angUpFromLevel)
%%GEOHEADING2UVEC Convert geographic headings in radians East of true North
%                 with elevations above the local tangent plane at a
%                 particular latitude and longitude with respect to a
%                 reference ellipsoid to unit vectors in ECEF coordinates.
%                 This function is useful for determining the Cartesian
%                 direction of a target when one wishes to simulate its
%                 motion in a particular direction over a curved Earth.
%
%INPUTS: point The 2XN or 3XN location matrix of the points in geodetic
%              latitude  and longitude in radians at which the headings are
%              taken. The point can be [latitude;longitude] or
%              [latitude;longitude;height]. The height component is ignored
%              if included because it does not change the result. If all of
%              the points are the same, then a single 2X1 or 3X1 vector can
%              be passed.
% geoEastOfNorth An NX1 or 1XN array of N geographic headings in radians
%              clockwise from North that should be turned into ECEF unit
%              vectors. A geographic heading is a direction in the local
%              tangent plane of an East-North-Up coordinate system as
%              defined on a specific reference ellipsoid. If all headings
%              are the same (but the angles up from level vary), then this
%              can be a scalar value. If this parameter is omitted or an
%              empty matrix is passed, then the default of 0 is used.
% angUpFromLevel An NX1 or 1XN array of N angles of the trajectory above
%              the local tangent plane to the reference ellipsoid. If all
%              elevations are the same (but geographic headings might vary),
%              then this can be a scalar value. If this parameter is
%              omitted or an empty matrix is passed, then the default of 0
%              is used.
%
%OUTPUTS: u A 3XN matrix of unit vectors pointing in the directions
%            implied by the headings in geoEastOfNorth.
%
%A heading can be turned into components of a local East-North-Up
%coordinate system. sin(geoEastOfNorth) is the component East and
%cos(geoEastOfNorth) is the component North. Thus, a unit vector in the
%appropriate direction in ECEF coordinates can be obtained by multiplying
%the components by the unit vectors in the appropriate directions. However,
%if one wishes to also be elevated above the local tangent plane, then it
%can be seen that one is defining  a spherical coordinate system, so the
%East and North components should additionaly be multiplied by
%cos(angUpFromLevel) and the vertical component should be added and
%multiplied by sin(angUpFromLevel). The East, North, and Up unit vectors
%are obtained using the getENUAxes function.
%
%EXAMPLE:
%The unit vector for a target heading can also be obtained by getting the
%rotation matrix to point a radar in a particular direction and rotating
%the z axis (since the convention used in findRFTransParam is that the
%local z axis is the boresight direction of the radar. This example just
%shows that both methods produce the same result.
% point=[0.1;-0.2];%latitude longitude location in radians.
% geoEastOfNorth=-120*(pi/180);
% angUpFromLevel=-4*(pi/180);
% u1=geogHeading2uVec(point,geoEastOfNorth,angUpFromLevel);
% M=findRFTransParam(point,geoEastOfNorth,angUpFromLevel);
% u2=M'*[0;0;1];%The transpose is the inverse rotation.
% norm(u1-u2)
%The result of norm(u1-u2) should be exactly zero.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(angUpFromLevel))
    angUpFromLevel=0; 
end

if(nargin<2||isempty(geoEastOfNorth))
    geoEastOfNorth=0; 
end

Np=size(point,2);
N=max([length(geoEastOfNorth);length(angUpFromLevel);Np]);

if(Np==1)
    point=repmat(point,[1,N]);
end

if(isscalar(geoEastOfNorth))
    geoEastOfNorth=geoEastOfNorth*ones(1,N);
end

if(isscalar(angUpFromLevel))
    angUpFromLevel=angUpFromLevel*ones(1,N);
end

%Transpose geoEastOfNorth to make it 1XN in size if it is NX1. This is so
%that bsxfun can be used to multiply the components with the basis vectors
%instead of using a loop.
geoEastOfNorth=geoEastOfNorth(:).';
%The same thing applies to angUpFromLevel
angUpFromLevel=angUpFromLevel(:).';

u=zeros(3,N);
for k=1:N
    %Get the unit vectors for East, North, and Up.
    uLocal=getENUAxes(point(:,k),false);
    uEast=uLocal(:,1);
    uNorth=uLocal(:,2);
    uUp=uLocal(:,3);
    
    cosEl=cos(angUpFromLevel(k));
    sinEl=sin(angUpFromLevel(k));
    
    u(:,k)=sin(geoEastOfNorth(k)).*cosEl*uEast+cos(geoEastOfNorth(k))*cosEl*uNorth+sinEl*uUp;
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
