function [aDeriv,aJacob]=aCVRuv(x,t)
%ACVRUV The drift function for a continuous-time motion model where the
%       the target state is given in monostatic range and u-v direction
%       cosine coordinates and the motion is constant-velocity in 3D
%       Cartesian coordinates.
%
%INPUTS: x The 6XN state vector of N targets in the order of
%          [r;u;;v;rDot;uDot;vDot], where r is monostatic range and u and v
%          are direction cosines. The target is assumed to always be in
%          front of the radar and not pass behind the radar.
%        t An unused time component so that aCVSpherical can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
%
%OUTPUTS: aDeriv The 6XN set of time-derivatives of the N state vectors
%                under the Cartesian linear motion model in r-u-v
%                coordinates.
%         AJacob The 6X6XN set of Jacobians of aDeriv, which can be useful
%                in extended Kalman filters. This is the derivative of each
%                component of aDeriv (selected by row) with respect to the
%                elements of the state [r,u,v,rDot,uDot,vDot] selected by
%                column.
%
%The derivation of the model is quite lengthy. A summary is given here. Let
%rVec be a position vector. We shall use the orthonormal basis vectors:
%uVec1=[u;v;sqrt(1-u^2-v^2)]
%uVec2=[sqrt((1-u^2-v^2)/(1-v^2));0;-u/sqrt(1-v^2)]
%uVec3=[-u*v/sqrt(1-v^2);sqrt(1-v^2);-v*sqrt((1-u^2-v^2)/(1-v^2))];
%Note that uVec2 is a normalized version of the derivative of uVec1 with
%respect to u and uVec3 is cross(u1,u2). The time-derivatives of the basis
%vectors can themselves be expressed as basis vectors:
%uVec1Dot=c1*u2+c2*u3
%uVec2Dot=c3*u1+c4*c3
%uVec3Dot=c5*u1+c6*u2
%The coefficients c1 through c6 are omitted for brevity but can be found
%by taking the total derivative of each basis vector with respect to time
%and then taking the dot products of the result with the orthonormal basis
%vectors.
%The vector rVec is written
%rVec=r*u1
%The velocity vector is the derivative of rVec. Using the chain rule and
%substituting in the derivative parameters for the basis vectors,
%rVecDot=rDot*uVec1+r*c1*uVec2+r*c2*uVec3
%The decond derivative is more complicated, but simplifies to
%rVecDDot=a1*uVec1+a2*uVec2+a3*uVec3
%where two D's indicates two time-derivatives.
%a1 is a radial acceleration and a3 and a3 are two other Cartesian
%components of acceleration. The expressions for them are
% a1=rDDot+r*(c1*c3+c2*c5);
% a2=2*rDot*c1+r*(c1Dot+c2*c6);
% a3=2*rDot*c2+r*(c2Dot+c1*c4);
%The expression for the drift function comes from setting a1=a2=a3=0 and
%solving for the second derivative terms. The second derivatives in u and v
%arise through the derivatives of c1 and c2.
%
%EXAMPLE:
%Here, we verify that integrating forward with this model is equivalent to
%linear motion in Cartesian coordinates.
% xInitCart=[400;920;1000;-100;20;64];
% T=50;%The prediction time.
% F=FPolyKal(T,xInitCart,1);
% xInitRuv=stateCart2Ruv(xInitCart);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsRuv=RKAdaptiveOverRange(xInitRuv,[0;T],@(x,t)aCVRuv(x,t),0.01,0,[],[],RelTol,AbsTol);
% xEndRuvRK=xStepsRuv(:,end);
% xEndRuvExact=stateCart2Ruv(xEndCart);
% max(abs(xEndRuvRK-xEndRuvExact)./xEndRuvExact)
%One will observe that the error is less than 1e-11, which is a good
%agreement.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

r=x(1,:);
u=x(2,:);
v=x(3,:);
rDot=x(4,:);
uDot=x(5,:);
vDot=x(6,:);

w2=max(0,1-u.^2-v.^2);
w=sqrt(w2);

diffV2=1-v.^2;
diffV=sqrt(diffV2);
denomV2=w2.*diffV2;
denomV=sqrt(denomV2);

diffU2=1-u.^2;
diffU=sqrt(diffU2);
denomU2=w2.*diffU2;
denomU=sqrt(denomU2);

c1v=(uDot.*diffV2+u.*v.*vDot)./denomV;
c1u=(vDot.*diffU2+u.*v.*uDot)./denomU;
c2v=vDot./diffV;
c2u=uDot./diffU;
c4v=v.*(-u.^2.*(1./w)-w).*(-uDot.*diffV2-u.*v.*vDot)./diffV2.^2;
c4u=u.*(-v.^2.*(1./w)-w).*(-vDot.*diffU2-u.*v.*uDot)./diffU2.^2;

d2v=2*rDot.*c2v+r.*c1v.*c4v;
d2u=2*rDot.*c2u+r.*c1u.*c4u;

rDDot=r.*((uDot.^2.*diffV2+2*u.*uDot.*v.*vDot+diffU2.*vDot.^2)./w2);

uDDot=-d2u.*diffU./r-u.*uDot.^2./diffU2;
vDDot=-d2v.*diffV./r-v.*vDot.^2./diffV2;

aDeriv=[rDot;uDot;vDot;rDDot;uDDot;vDDot];

if(nargout>1)
    ddrDDdr=-(-uDot.^2.*diffV2-2*u.*uDot.*v.*vDot-diffU2.*vDot.^2)./w2;
    ddrDDdu=(2*r.*(u.*uDot+v.*vDot).*(uDot-uDot.*v.^2+u.*v.*vDot))./w2.^2;
    ddrDDdv=-((2*r.*(-u.*uDot.*v-diffU2.*vDot).*(u.*uDot+v.*vDot))./w2.^2);

    ddrDDduDot=((2*r.*(uDot-uDot.*v.^2+u.*v.*vDot))./w2);
    ddrDDdvDot=-(2*r.*(-u.*uDot.*v-diffU2.*vDot))./w2;

    dduDDdr=(2*rDot.*uDot)./r.^2;
    dduDDdu=(-uDot.^2.*(1+u.^2-v.^2).*diffV2-4*u.*uDot.*v.*diffV2.*vDot+(-1+v.^2-u.^2.*(-2+u.^2+3*v.^2)).*vDot.^2)./w2.^2;
    dduDDdv=(2*u.*(-u.*uDot.*v-diffU2.*vDot).*(u.*uDot+v.*vDot))./w2.^2;
    dduDDdrDot=-((2*uDot)./r);
    dduDDduDot=-((2*rDot)./r)-(2*u.*(uDot-uDot.*v.^2+u.*v.*vDot))./w2;
    dduDDdvDot=-(2*u.*(u.*uDot.*v+vDot-u.^2.*vDot))./w2;

    ddvDDdr=(2*rDot.*vDot)./r.^2;
    ddvDDdu=-((2*v.*(u.*uDot+v.*vDot).*(uDot-uDot.*v.^2+u.*v.*vDot))./w2.^2);
    ddvDDdv=(uDot.^2.*(u.^2.*(1-3*v.^2)-diffV2.^2)-4*u.*diffU2.*uDot.*v.*vDot+diffU2.*(-1+u.^2-v.^2).*vDot.^2)./w2.^2;
    ddvDDdrDot=-((2*vDot)./r);
    ddvDDduDot=-(2*v.*(uDot-uDot.*v.^2+u.*v.*vDot))./w2;
    ddvDDdvDot=-((2*rDot)./r)-(2*v.*(u.*uDot.*v+vDot-u.^2.*vDot))./w2;

    
    N=size(x,2);
    aJacob=zeros(6,6,N);
    for k=1:N
        aJacob(:,:,k)=[0,               0,              0,              1,              0,              0;
                       0,               0,              0,              0,              1,              0;
                       0,               0,              0,              0,              0,              1;
                       ddrDDdr(k),      ddrDDdu(k),     ddrDDdv(k),     0,  ddrDDduDot(k),  ddrDDdvDot(k);
                       dduDDdr(k),      dduDDdu(k),     dduDDdv(k),     dduDDdrDot(k),  dduDDduDot(k),  dduDDdvDot(k);
                       ddvDDdr(k),      ddvDDdu(k),     ddvDDdv(k),     ddvDDdrDot(k),  ddvDDduDot(k),  ddvDDdvDot(k)];
    end
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
