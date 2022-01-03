function CartStates=stateRu2D2Cart(x)
%%STATERU2D2CART Convert a target state in 2D space in monostatuc r-u
%               coordinates with first and possibly second time derivatives
%               of the position components into Cartesian coordinates. The
%               state has the format [r;u;rDot;uDot;rDDot;uDDot], where two
%               Ds indicate a second derivative with respect to time. In
%               Cartesian coordinates, the converted state has the
%               form [x;y;xDot;yDot;xDDot;yDDot]. It is assumed that the
%               target is in front of the sensor so the unobserved second
%               component of a u-v unit vector in 2D is positive.
%
%INPUTS: x The 4XN or 6XN set of 2D Ru state vectors consisting of
%          position, velocity and possibly acceleration. The angles are
%          given in radians. The range is a one-way range.
%
%OUTPUTS: CartStates The 4XN or 6XN set of target states given in Cartesian
%                    coordinates.
%
%The derivation is given in [1]. The function aCVRu2D is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateCart2Ru2D.
% x=[100;60;3;12;108;-116];
% xRet=stateRu2D2Cart(stateCart2Ru2D(x));
% max(abs(x(:)-xRet(:)))
%One will see that the error is less than 1e-13, indicating good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic linear Cartesian dynamic models in local
%    coordinates," Naval Research Laboratory, Washington, DC, Tech. Rep.
%    NRL/MR/5344-19-9882, 24 Aug. 2019.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);
N=size(x,2);

CartStates=zeros(numDim,N);

r=x(1,:);
u=x(2,:);
rDot=x(3,:);
uDot=x(4,:);

v2=max(0,1-u.^2);
v=sqrt(v2);
uVecr=[u;v];
uVecu=[v;-u];

%Position
CartStates(1:2,:)=r.*uVecr;

%Velocity
CartStates(3:4,:)=bsxfun(@times,rDot,uVecr)+bsxfun(@times,r.*uDot./v,uVecu);

if(numDim>4)
    rDDot=x(5,:);
    uDDot=x(6,:);

    ar=rDDot-r.*uDot.^2./v2;
    au=(2*rDot.*uDot+r.*uDDot)./v+r.*u.*uDot.^2./v.^3;

    %Acceleration
    CartStates(5:6,:)=bsxfun(@times,ar,uVecr)+bsxfun(@times,au,uVecu);
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
