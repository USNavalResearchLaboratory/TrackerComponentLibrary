function CartPoints=TDOAPolar2D2Cart(points,systemType,lRef,lRx,M,c,mirrorNegative)
%%TDOAPOLAR2D2CART Given positions specified as a time-difference-of-
%     arrival (TDOA) and a local azimuth measurement in 2D, convert the
%     positions to global Cartesian coordinates. Note that more than one
%     point in TDOA-azimuth coordinates maps to the same Cartesian point.
%
%INPUTS: points A 2XN set of points of the form [TDOA;azimuth] with the
%           angles given in radians.
%  systemType An optional parameter specifying the axes from which the
%           angles are measured. Possible values are
%           0 (The default if omitted) The azimuth angle is
%             counterclockwise from the x axis.
%           1 The azimuth angle is measured clockwise from the y axis.
%      lRef The 2X1 location of the reference receiver in global
%           coordinates.
%       lRx The 2X1 location of the receiver in global coordinates at which
%           the direction of arrival in spherical coordinates is locally
%           defined and whose time of arrival (TOA) has the TOA at lRef
%           subtracted to get the TDOA.
%         M A 2X2 rotation matrix to go from the alignment of the global
%           coordinate system to that at the receiver. If omitted, then it
%           is assumed that the local coordinate system is aligned with the
%           global and M=eye(2) --the identity matrix is used.
%         c The propagation speed in the medium in question. If this
%           parameter is omitted or an empty matrix is passed, the default
%           value of Constants.speedOfLight is used.
% mirrorNegative If true, then if the computed one-way range to the target
%           is negative (and thus invalid), use a heuristic to change the
%           transformation and make the result usually more accurate. The
%           default if this is omitted or an empty matrix is passed is
%           false.
%
%OUTPUTS: CartPoints The 2XN set of Cartesian locations corresponding to
%                    the input points.
%
%A derivation of a similar coordinate system in 3D is given in [1]. The
%derivation of this follows in a similar manner.
%
%EXAMPLE:
%We convert a Cartesian position to TDOA-polar coordinates and then we
%convert it back. The absolute error of the round-trip conversion is within
%finite precision limits.
% lRef=[-3;0];
% lRx=[3;0];
% xTar=[4;4];
% systemType=0;
% z=Cart2TDOAPolar2D(xTar,systemType,lRef,lRx);
% AbsErr=norm(TDOAPolar2D2Cart(z,systemType,lRef,lRx)-xTar)
%
%REFERENCES:
%[1] D. F. Crouse, "Particle flow solutions avoiding stiff integration,"
%    U.S. Naval Research Laboratory, Washington, DC, Tech. Rep.
%    NRL/5340/FR-2021/1, 25 May 2021.
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(mirrorNegative))
    mirrorNegative=false;
end

if(nargin<6||isempty(c))
    c=Constants.speedOfLight; 
end

if(nargin<5||isempty(M))
    M=eye(2,2); 
end

if(isempty(systemType))
    systemType=0;
end

rDiff=points(1,:)*c;
azimuth=points(2,:);

%Get unit vectors pointing in the direction of the targets from the
%receiver in the LOCAL coordinate system of the receiver.
N=size(points,2);
u=zeros(2,N);
switch(systemType)
    case 0
        u(1,:)=cos(azimuth);
        u(2,:)=sin(azimuth);
    case 1
        u(1,:)=sin(azimuth);
        u(2,:)=cos(azimuth);
    otherwise
        error('Invalid system type specified.')
end

%Convert the direction vectors into global coordinates.
u=M'*u;

lRefRxDiff=lRef-lRx;

r1=(rDiff.^2-sum(lRefRxDiff.*lRefRxDiff,1))./(2*(rDiff-sum(bsxfun(@times,u,lRefRxDiff),1)));

if(mirrorNegative)
    numPts=length(r1);
    CartPoints=zeros(2,numPts);
    for k=1:numPts
        if(r1(k)<0)
            CartPoints(:,k)=-r1(k)*u(:,k)+lRef;
        else
            CartPoints(:,k)=r1(k)*u(:,k)+lRx;
        end
    end
else
    CartPoints=bsxfun(@plus,lRx,bsxfun(@times,r1,u));
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
