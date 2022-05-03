function vECEF=ENUVec2ECEFVec(latLonOrigin,vENU,a,f)
%%ENUVEC2ECEFVEC Given a direction vector, such as a velocity vector, in
%           local East-North-Up coordinates with respect to a reference
%           ellipsoid, rotate the vector to be a unit vector in global ECEF
%           coordinates. If a 2D vector is passed, it is assumed to just be
%           in the local East-North tangent plan with the up component
%           zero.
%
%INPUTS: latLonOrigin The 2X1 [latitude; longitude] point in radians
%             defining the local tangent plane with respect to a reference
%             ellipsoid. If a 3X1 vector is passed, the third row is
%             ignored (for example, if one were to pass a
%             [latitude;longitude;height] point.
%        vENU A 3XN set of [East;North;Up] vectors of a 2XN set of
%             [East;North] vectors to rotate into the global coordinate
%             system.
%           a The semi-major axis of the reference ellipsoid. If this
%             argument is omitted, the value in
%             Constants.WGS84SemiMajorAxis is used.
%           f The flattening factor of the reference ellipsoid. If this
%             argument is omitted, the value in Constants.WGS84Flattening
%             is used.
%
%OUTPUTS: vECEF The 3XN set of the vENU vectors rotated into 3D ECEF
%               coordinates.
%
%Such a conversion can be inferred from the rotations in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

uENU=getENUAxes(latLonOrigin,false,a,f);
if(size(vENU,1)==2)
   %Local East-North vectors. 
   vECEF=bsxfun(@times,vENU(1,:),uENU(:,1))+bsxfun(@times,vENU(2,:),uENU(:,2));
else
   %Local East-North-Up vectors.
    vECEF=bsxfun(@times,vENU(1,:),uENU(:,1))+bsxfun(@times,vENU(2,:),uENU(:,2))+bsxfun(@times,vENU(3,:),uENU(:,3));
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
