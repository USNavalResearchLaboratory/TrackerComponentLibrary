function vENU=ECEFVec2ENUVec(latLonOrigin,vECEF,onlyEN,a,f)
%ECEFVEC2ENUVEC Given a direction vector, such as a velocity vector, in a
%           global Earth-Centered-Earth fixed coordinate system, obtain the
%           East-North-Up coordinates of the same vector in the local
%           tangent plane defined at a particular point on a reference
%           ellipsoid.
%
%INPUTS: latLonOrigin The [latitude; longitude] point in radians defining
%             the local tangent plane with respect to a reference
%             ellipsoid.
%       vECEF A 3XN set of [x;y;z] points in the global ECEF coordinate
%             system.
%      onlyEN If the output should be a 2XN set of [East;North] components
%             and no Up, then this is true. Otherwise, the output will be a
%             3XN set of [East;North; Up] components. If omitted or an
%             empty matrix is passed, then the default is false.
%           a The semi-major axis of the reference ellipsoid. If this
%             argument is omitted, the value in
%             Constants.WGS84SemiMajorAxis is used.
%           f The flattening factor of the reference ellipsoid. If this
%             argument is omitted, the value in Constants.WGS84Flattening
%             is used.
%
%OUTPUTS: vENU If onlyEN=false, this is a 3XN matrix of [East;North;Up]
%              vectors. Otherwise, this is a 2XN set of [East;North]
%              vectors.
%
%Such a conversion can be inferred from the rotations in [1].
%
%EXAMPLE:
%Here, we demonstrate that ECEFVec2ENUVec and ENUVec2ECEFVec form a
%consistent pair. A vector is converted from ENU to ECEF and then back. The
%back-converted point agreess without forwar one with a relative error
%indicative of finite precision limits.
% latLonOrigin=deg2rad([24.266997;-152.244802]);
% vEN=[[30;-12],[-3;-1]];
% is2D=true;
% vENBack=ECEFVec2ENUVec(latLonOrigin,ENUVec2ECEFVec(latLonOrigin,vEN),is2D);
% RelErr=max(abs((vEN(:)-vENBack(:))./vEN(:)))
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<4||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<3||isempty(onlyEN))
    onlyEN=false;
end

uENU=getENUAxes(latLonOrigin,false,a,f);
N=size(vECEF,2);
if(onlyEN)
    vENU=zeros(2,N);
    vENU(1,:)=sum(bsxfun(@times,vECEF,uENU(:,1)),1);
    vENU(2,:)=sum(bsxfun(@times,vECEF,uENU(:,2)),1);
else
    vENU=zeros(3,N);
    vENU(1,:)=sum(bsxfun(@times,vECEF,uENU(:,1)),1);
    vENU(2,:)=sum(bsxfun(@times,vECEF,uENU(:,2)),1);
    vENU(3,:)=sum(bsxfun(@times,vECEF,uENU(:,3)),1);
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
