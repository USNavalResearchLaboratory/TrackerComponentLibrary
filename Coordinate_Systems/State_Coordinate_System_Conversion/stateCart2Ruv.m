function ruvStates=stateCart2Ruv(x)
%%STATECART2RUV Convert a target state in 3D space consiting of
%               position, velocity and possibly acceleration terms into
%               monostatic r-u-v coordinates. With acceleration, this
%               would be going from
%               [x;y;z;xDot;yDot;zDot;xDDot;yDDot;zDDot] to
%               [r;u;v;rDot;uDot;vDot;rDDot;uDDot;vDDot]
%               where two Ds indicates a second derivative with respect to
%               time.
%
%INPUTS: x The 6XN or 9XN set of Cartesian state vectors consisting of
%          position, velocity and possibly acceleration. It is assumed that
%          the z component is positive, which is consistent with
%          measurements being in front of the radar.
%
%OUTPUTS: ruvStates The 6XN or 9XN set of target states given in r-u-v
%                   coordinates. The range is a one-way range.
%
%The derivation is given in [1]. The function aCVRuv is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateRuv2Cart.
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
ruvStates=zeros(numDim,N);

posRUV=Cart2Ruv(x(1:3,:),true);
ruvStates(1:3,:)=posRUV;

r=posRUV(1,:);
u=posRUV(2,:);
v=posRUV(3,:);

w2=max(0,1-u.^2-v.^2);
w=sqrt(w2);

diffV2=1-v.^2;
diffV=sqrt(diffV2);

diffU2=1-u.^2;
diffU=sqrt(diffU2);

u1=[u;v;w];
u2u=[zeros(1,N);w./diffU;-v./diffU];
u2v=[w./diffV;zeros(1,N);-u./diffV];

u3u=[diffU;-u.*v./diffU;-u.*(w./diffU)];
u3v=[-u.*v./diffV;diffV;-v.*(w./diffV)];

rVecVel=x(4:6,:);
rDot=sum(rVecVel.*u1,1);

c2u=sum(rVecVel.*u3u,1)./r;
c2v=sum(rVecVel.*u3v,1)./r;

uDot=c2u.*diffU;
vDot=c2v.*diffV;

ruvStates(4:6,:)=[rDot;uDot;vDot];

if(numDim>6)
    %If acceleration is provided.
    rVecAccel=x(7:9,:);

    a1=sum(rVecAccel.*u1,1);
    a3u=sum(rVecAccel.*u3u,1);
    a3v=sum(rVecAccel.*u3v,1);
    
    c1u=sum(rVecVel.*u2u,1)./r;
    c1v=sum(rVecVel.*u2v,1)./r;
    c2v=vDot./diffV;
    c2u=uDot./diffU;
    c4v=v.*(-u.^2.*(1./w)-w).*(-uDot.*diffV2-u.*v.*vDot)./diffV2.^2;
    c4u=u.*(-v.^2.*(1./w)-w).*(-vDot.*diffU2-u.*v.*uDot)./diffU2.^2;
    
    d2v=2*rDot.*c2v+r.*c1v.*c4v;
    d2u=2*rDot.*c2u+r.*c1u.*c4u;
    
    rDDot=a1+r.*((uDot.^2.*diffV2+2*u.*uDot.*v.*vDot+diffU2.*vDot.^2)./w2);

    uDDot=(a3u-d2u).*diffU./r-u.*uDot.^2./diffU2;
    vDDot=(a3v-d2v).*diffV./r-v.*vDot.^2./diffV2;

    ruvStates(7:9,:)=[rDDot;uDDot;vDDot];
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
