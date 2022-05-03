function CartPoints=TDOASpher2Cart(points,systemType,lRef,lRx,M,c)
%%TDOASPHER2CART Given positions specified as a time-difference-of-arrival
%     (TDOA) and local spherical azimuth and elevation, convert the
%     positions to global Cartesian coordinates. Note that more than one
%     point in TDOA-spherical coordinates maps to the same Cartesian point.
%
%INPUTS: points A 3XN set of points of the form [TDOA;azimuth;elevation]
%           with the angles given in radians.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%      lRef The 3X1 location of the reference receiver in global
%           coordinates.
%       lRx The 3X1 location of the receiver in global coordinates at which
%           the direction of arrival in spherical coordinates is locally
%           defined and whose time of arrival (TOA) has the TOA at lRef
%           subtracted to get the TDOA.
%         M A 3X3 rotation matrix to go from the alignment of the global
%           coordinate system to that at the receiver. If omitted, then it
%           is assumed that the local coordinate system is aligned with the
%           global and M=eye(3) --the identity matrix is used.
%         c The propagation speed in the medium in question. If this
%           parameter is omitted or an empty matrix is passed, the default
%           value of Constants.speedOfLight is used.
%
%OUTPUTS: CartPoints The 3XN set of Cartesian locations corresponding to
%                    the input points.
%
%A derivation of this coordinate system is given in [1].
%
%EXAMPLE:
%We convert a Cartesian position to TDOA-spherical coordinates and then we
%convert it back. The absolute error of the round-trip conversion is within
%finite precision limits.
% lRef=[-3;0;0];
% lRx=[3;0;0];
% xTar=[4;4;10];
% systemType=0;
% z=Cart2TDOASphere(xTar,systemType,lRef,lRx)
% AbsErr=norm(TDOASpher2Cart(z,systemType,lRef,lRx)-xTar)
%
%REFERENCES:
%[1] D. F. Crouse, "Particle flow solutions avoiding stiff integration,"
%    U.S. Naval Research Laboratory, Washington, DC, Tech. Rep.
%    NRL/5340/FR-2021/1, 25 May 2021.
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(c))
    c=Constants.speedOfLight; 
end

if(nargin<5||isempty(M))
    M=eye(3,3); 
end

if(isempty(systemType))
    systemType=0;
end

rDiff=points(1,:)*c;
azimuth=points(2,:);
elevation=points(3,:);

if(systemType==2)
    elevation=pi/2-elevation;
    systemType=0;
elseif(systemType==3)
    azimuth=pi/2-azimuth;
    systemType=0;
end

%Get unit vectors pointing in the direction of the targets from the
%receiver in the LOCAL coordinate system of the receiver.
N=size(points,2);
uVecL=zeros(3,N);
switch(systemType)
    case 0
        uVecL(1,:)=cos(azimuth).*cos(elevation);
        uVecL(2,:)=sin(azimuth).*cos(elevation);
        uVecL(3,:)=sin(elevation);
    case 1
        uVecL(1,:)=sin(azimuth).*cos(elevation);
        uVecL(2,:)=sin(elevation);
        uVecL(3,:)=cos(azimuth).*cos(elevation);
    otherwise
        error('Invalid system type specified.')
end

%Convert the direction vectors into global coordinates.
uVecL=M'*uVecL;

lRefRxDiff=lRef-lRx;

r1=(rDiff.^2-sum(lRefRxDiff.*lRefRxDiff,1))./(2*(rDiff-sum(bsxfun(@times,uVecL,lRefRxDiff),1)));
CartPoints=bsxfun(@plus,lRx,bsxfun(@times,r1,uVecL));

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
