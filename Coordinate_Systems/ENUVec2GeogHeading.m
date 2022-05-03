function [heading,angUpFromLevel]=ENUVec2GeogHeading(u)
%%ENUVEC2GEOGHEADING Given a 3D direction vector in local East-North-Up
%           coordinates, get the equivalent heading and elevation in
%           radians East of North and radians above the local  tangent
%           plane.
%
%INPUTS: u A 3XnumVec set of direction vectors in the local ENU coordinate
%          system.
%
%OUTPUTS: heading A 1XnumVec set of headings in radians East of North
%                 corresponding to the values in u.
%  angUpFromLevel A 1XnumVec set of elevations above the local tangent
%                 plane in radians corresponding to the values in u.
%
%The transformation is just a specific definition of the spherical
%coordinate system. There is no need for information on the shape of the
%reference ellipsoid.
%
%EXAMPLE:
%Here, we show that the values produced by this function are consistent
%with those returned by the geogHeading2uVec function. The difference
%between the original heading and elevation and converted back values are
%on the order of finite precision errors.
% headingInit=deg2rad(20);
% angUpFromLevelInit=deg2rad(5);
% 
% latLon=deg2rad([20.775632;-155.975253]);
% uECEF=geogHeading2uVec(latLon,headingInit,angUpFromLevelInit);
% %rotation matrix from local ENU coordinates to global coordinates.
% M=getENUAxes(latLon);
% %Matrix transpose goes from global to local (ENU) coordinates.
% uENU=M'*uECEF;
% [heading,angUpFromLevel]=ENUVec2GeogHeading(uENU);
% abs(heading-headingInit)
% abs(angUpFromLevel-angUpFromLevelInit)
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

x=u(1,:);
y=u(2,:);
z=u(3,:);

heading=atan2(x,y);
r=x.^2+y.^2+z.^2;
angUpFromLevel=asin(z./r);

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
