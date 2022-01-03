function CartStates=stateRuv2Cart(x)
%%STATERUV2CART Convert a target state in 3D space in monostatic
%               r-u-v coordinates with first and possibly second time
%               derivatives of the position components into Cartesian
%               coordinates. The state has the format
%               [r;u;v;rDot;uDot;vDot;rDDot;uDDot;vDDot],
%               where two Ds indicate a second derivative with respect to
%               time. In Cartesian coordinates, the converted state has the
%               form [x;y;z;xDot;yDot;zDot;xDDot;yDDot;zDDot]. The target
%               is assumed to be in front of the radar.
%
%INPUTS: x The 6XN or 9XN set of r-u-v state vectors consisting of
%          position, velocity and possibly acceleration.  The range is a
%          one-way range.
%
%OUTPUTS: CartStates The 6XN or 9XN set of target states given in Cartesian
%                    coordinates.
%
%The derivation is given in [1]. The function aCVRuv is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateCart2Ruv.
% x=[100;-60;200;-3;12;160;108;-116;2];
% xRet=stateRuv2Cart(stateCart2Ruv(x));
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
v=x(3,:);
rDot=x(4,:);
uDot=x(5,:);
vDot=x(6,:);

w2=max(0,1-u.^2-v.^2);
w=sqrt(w2);
diffV2=1-v.^2;
diffV=sqrt(1-v.^2);
denom2=w2.*diffV2;
denom=sqrt(denom2);

u1=[u;v;w];
u2=[w./diffV;zeros(1,N);-u./diffV];
u3=[-u.*v./diffV;diffV;-v.*(w./diffV)];

%Position
CartStates(1:3,:)=bsxfun(@times,r,u1);

c1=(uDot.*diffV2+u.*v.*vDot)./denom;
c2=vDot./diffV;

%Velocity
CartStates(4:6,:)=bsxfun(@times,rDot,u1)+bsxfun(@times,r.*c1,u2)+bsxfun(@times,r.*c2,u3);

if(numDim>6)
    %Acceleration
    rDDot=x(7,:);
    uDDot=x(8,:);
    vDDot=x(9,:);

    c3=-(((w+u.^2./w).*(uDot-uDot.*v.^2+u.*v.*vDot))./diffV.^3);
    c4=v.*(-u.^2.*(1./w)-w).*(-uDot.*diffV2-u.*v.*vDot)./diffV2.^2;
    c5=-vDot./diffV;
    c6=(w./diffV).*(-v.*uDot.*diffV2-u.*vDot)./diffV.^3+u.*(diffV./w).*(-u.*v.*uDot.*diffV2-u.^2.*vDot+vDot.*diffV2.^2)./diffV.^5;
   
    c1Dot=((uDot.*diffV2+u.*v.*vDot).*(v.*vDot.*(2-u.^2-2*v.^2)+u.*uDot.*diffV2))./denom.^3 ...
          +(uDDot.*diffV2-v.*uDot.*vDot+u.*(v.*vDDot+vDot.^2))./denom;
    c2Dot=(vDDot.*diffV2+v.*vDot.^2)./diffV.^3;

    a1=rDDot+r.*(c1.*c3+c2.*c5);
    a2=2*rDot.*c1+r.*(c1Dot+c2.*c6);
    a3=2*rDot.*c2+r.*(c2Dot+c1.*c4);
    
    CartStates(7:9,:)=bsxfun(@times,a1,u1)+bsxfun(@times,a2,u2)+bsxfun(@times,a3,u3);
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
