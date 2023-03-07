function [aDeriv,aJacob,aHess,papt]=aCVRuv(x)
%ACVRUV The drift function for a continuous-time motion model where the
%       the target state is given in monostatic range and u-v direction
%       cosine coordinates and the motion is constant-velocity in 3D
%       Cartesian coordinates.
%
%INPUTS: x The 6XN state vector of N targets in the order of
%          [r;u;;v;rDot;uDot;vDot], where r is monostatic range and u and v
%          are direction cosines. The target is assumed to always be in
%          front of the radar and not pass behind the radar.
%
%OUTPUTS: aDeriv The 6XN set of time-derivatives of the N state vectors
%                under the Cartesian linear motion model in r-u-v
%                coordinates.
%         aJacob The 6X6XN set of Jacobians of aDeriv, which can be useful
%                in extended Kalman filters. This is the derivative of each
%                component of aDeriv (selected by row) with respect to the
%                elements of the state [r,u,v,rDot,uDot,vDot] selected by
%                column.
%
%A full derivation of the dynamic model is provided in [1]. Summarizing it
%here, let rVec be a position vector. We shall use the orthonormal basis
%vectors:
%uVec1=[u;v;sqrtw]
%uVec2=[sqrt(w/(1-v2));0;-u/sqrt(1-v2)]
%uVec3=[-u*v/sqrt(1-v2);sqrt(1-v2);-v*sqrt(w/(1-v2))];
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
% F=FPolyKal(T,6,1);
% xInitRuv=stateCart2Ruv(xInitCart);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsRuv=RKAdaptiveOverRange(xInitRuv,[0;T],@(x,t)aCVRuv(x),0.01,0,[],[],RelTol,AbsTol);
% xEndRuvRK=xStepsRuv(:,end);
% xEndRuvExact=stateCart2Ruv(xEndCart);
% max(abs(xEndRuvRK-xEndRuvExact)./xEndRuvExact)
%One will observe that the error is less than 1e-11, which is a good
%agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic Linear Dynamic Models in Local Coordinates,"
%    Naval Research Laboratory, Washington, D.C., Tech. Rep.
%    NRL/MR/5344--19-9882, 24 Aug. 2019.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

r=x(1,:);
u=x(2,:);
v=x(3,:);
rDot=x(4,:);
uDot=x(5,:);
vDot=x(6,:);

u2=u.^2;
v2=v.^2;
r2=r.^2;

w2=max(0,1-u2-v2);
w=sqrt(w2);

diffV2=1-v2;
diffV=sqrt(diffV2);

diffU2=1-u2;
diffU=sqrt(diffU2);

uDot2=uDot.*uDot;
vDot2=vDot.*vDot;
w3=w2.*w;

rDDot=-(r*(-uDot2*diffV2-2*u*uDot*v*vDot-diffU2*vDot2))/w2;
uDDot=-((2*rDot*uDot)/r)-u*(vDot2*diffU2+uDot2*diffV2+2*u*v*uDot*vDot)/w2;
vDDot=-((2*rDot*vDot)/r)-v*(uDot2*diffV2+vDot2*diffU2+2*u*v*uDot*vDot)/w2;

aDeriv=[rDot;uDot;vDot;rDDot;uDDot;vDDot];

if(nargout>1)
    N=size(x,2);
    if(N>1)
        error('Derivatives are only available for numPoints=1.')
    end

    w4=w2^2;

    ddrDDdr=-(-uDot2*diffV2-2*u*uDot*v*vDot-diffU2*vDot2)/w2;
    ddrDDdu=(2*r*(u*uDot+v*vDot)*(uDot-uDot*v2+u*v*vDot))/w4;
    ddrDDdv=-((2*r*(-u*uDot*v-diffU2*vDot)*(u*uDot+v*vDot))/w4);

    ddrDDduDot=((2*r*(uDot-uDot*v2+u*v*vDot))/w2);
    ddrDDdvDot=-(2*r*(-u*uDot*v-diffU2*vDot))/w2;

    dduDDdr=(2*rDot*uDot)/r2;
    dduDDdu=(-uDot2*(1+u2-v2)*diffV2-4*u*uDot*v*diffV2*vDot+(-1+v2-u2*(-2+u2+3*v2))*vDot2)/w4;
    dduDDdv=(2*u*(-u*uDot*v-diffU2*vDot)*(u*uDot+v*vDot))/w4;
    dduDDdrDot=-((2*uDot)/r);
    dduDDduDot=-((2*rDot)/r)-(2*u*(uDot-uDot*v2+u*v*vDot))/w2;
    dduDDdvDot=-(2*u*(u*uDot*v+vDot-u2*vDot))/w2;

    ddvDDdr=(2*rDot*vDot)/r2;
    ddvDDdu=-((2*v*(u*uDot+v*vDot)*(uDot-uDot*v2+u*v*vDot))/w4);
    ddvDDdv=(uDot2*(u2*(1-3*v2)-diffV2^2)-4*u*diffU2*uDot*v*vDot+diffU2*(-1+u2-v2)*vDot2)/w4;
    ddvDDdrDot=-((2*vDot)/r);
    ddvDDduDot=-(2*v*(uDot-uDot*v2+u*v*vDot))/w2;
    ddvDDdvDot=-((2*rDot)/r)-(2*v*(u*uDot*v+vDot-u2*vDot))/w2;

    aJacob=[0,               0,              0,              1,              0,              0;
            0,               0,              0,              0,              1,              0;
            0,               0,              0,              0,              0,              1;
            ddrDDdr,         ddrDDdu,        ddrDDdv,        0,              ddrDDduDot,     ddrDDdvDot;
            dduDDdr,         dduDDdu,        dduDDdv,        dduDDdrDot,     dduDDduDot,     dduDDdvDot;
            ddvDDdr,         ddvDDdu,        ddvDDdv,        ddvDDdrDot,     ddvDDduDot,     ddvDDdvDot];

    if(nargout>2)
        aHess=zeros(6,6,6);
        
        u4=u2*u2;
        v4=v2*v2;
        w6=w3*w3;
        diffU4=diffU2*diffU2;
        diffV4=diffV2*diffV2;

        ddrDDdrdr=0;
        ddrDDdudr=(2*(u*uDot+v*vDot)*(uDot-uDot*v2+u*v*vDot))/w4;
        ddrDDdvdr=-((2*(-u*uDot*v-diffU2*vDot)*(u*uDot+v*vDot))/w4);
        ddrDDduDotdr=((2*(uDot-uDot*v2+u*v*vDot))/w2);
        ddrDDdvDotdr=-(-2*u*uDot*v-2*diffU2*vDot)/w2;

        ddrDDdrdu=ddrDDdudr;
        ddrDDdudu=-(2*r*(-uDot2*(1+3*u2-v2)*diffV2-2*u*uDot*v*(3+u2-3*v2)*vDot+v2*(-1-3*u2+v2)*vDot2))/w6;
        ddrDDdvdu=(2*r*(uDot*(1+u2-v2)+2*u*v*vDot)*(2*u*uDot*v+vDot-u2*vDot+v2*vDot))/w6;
        ddrDDduDotdu=(-4*r*u*uDot*(-diffV2)+2*r*v*(1+u2-v2)*vDot)/w4;
        ddrDDdvDotdu=(2*r*v*(uDot*(1+u2-v2)+2*u*v*vDot))/w4;

        ddrDDdrdv=ddrDDdvdr;
        ddrDDdudv=ddrDDdvdu;
        ddrDDdvdv=-(2*r*(u2*uDot2*(-1+u2-3*v2)+2*u*uDot*v*(-3+3*u2-v2)*vDot+diffU2*(-1+u2-3*v2)*vDot2))/w6;
        ddrDDduDotdv=(2*r*u*(2*u*uDot*v+vDot-u2*vDot+v2*vDot))/w4;
        ddrDDdvDotdv=(2*r*u*uDot*(1-u2+v2)+4*r*diffU2*v*vDot)/w4;

        ddrDDdrdRDot=0;
        ddrDDdudRDot=0;
        ddrDDdvdRDot=0;
        ddrDDduDotdRDot=0;
        ddrDDdvDotdRDot=0;

        ddrDDdrduDot=ddrDDduDotdr;
        ddrDDduduDot=ddrDDduDotdu;
        ddrDDdvduDot=ddrDDduDotdv;
        ddrDDduDotduDot=(2*r*diffV2)/w2;
        ddrDDdvDotduDot=(2*r*u*v)/w2;

        ddrDDdrdvDot=ddrDDdvDotdr;
        ddrDDdudvDot=ddrDDdvDotdu;
        ddrDDdvdvDot=ddrDDdvDotdv;
        ddrDDduDotdvDot=ddrDDdvDotduDot;
        ddrDDdvDotdvDot=(2*r*diffU2)/w2;

        %%%
        dduDDdrdr=-((4*rDot*uDot)/r^3);
        dduDDdudr=0;
        dduDDdvdr=0;
        dduDDdrDotdr=(2*uDot)/r2;
        dduDDduDotdr=(2*rDot)/r2;
        dduDDdvDotdr=0;

        dduDDdrdu=dduDDdudr;
        dduDDdudu=-(2*u*uDot2*(3+u2-3*v2)*diffV2-4*uDot*v*diffV2*(-1-3*u2+v2)*vDot+2*u*v2*(3+u2-3*v2)*vDot2)/w6;
        dduDDdvdu=-(2*u2*uDot2*v*(3+u2-3*v2)+4*u*uDot*(1-v4+u2*(-1+3*v2))*vDot+2*v*(1-u4+(-1+3*u2)*v2)*vDot2)/w6;
        dduDDdrDotdu=0;
        dduDDduDotdu=-(2*diffV2*(uDot*(1+u2-v2)+2*u*v*vDot))/w4;
        dduDDdvDotdu=(-4*u*uDot*v*diffV2-2*(diffU4+(-1+3*u2)*v2)*vDot)/w4;

        dduDDdrdv=dduDDdvdr;
        dduDDdudv=dduDDdvdu;
        dduDDdvdv=-(2*u*(u2*uDot2*(1-u2+3*v2)+2*u*uDot*v*(3-3*u2+v2)*vDot-diffU2*(-1+u2-3*v2)*vDot2))/w6;
        dduDDdrDotdv=0;
        dduDDduDotdv=(2*u2*(-2*u*uDot*v+u2*vDot-(1+v2)*vDot))/w4;
        dduDDdvDotdv=(2*u*(u*uDot*(-1+u2-v2)-2*diffU2*v*vDot))/w4;

        dduDDdrdRDot=dduDDdrDotdr;
        dduDDdudRDot=dduDDdrDotdu;
        dduDDdvdRDot=dduDDdrDotdv;
        dduDDdrDotdRDot=0;
        dduDDduDotdRDot=-(2/r);
        dduDDdvDotdRDot=0;

        dduDDdrduDot=dduDDduDotdr;
        dduDDduduDot=dduDDduDotdu;
        dduDDdvduDot=dduDDduDotdv;
        dduDDdrDotduDot=dduDDduDotdRDot;
        dduDDduDotduDot=-((2*u*diffV2)/w2);
        dduDDdvDotduDot=-(2*u2*v)/w2;

        dduDDdrdvDot=dduDDdvDotdr;
        dduDDdudvDot=dduDDdvDotdu;
        dduDDdvdvDot=dduDDdvDotdv;
        dduDDdrDotdvDot=dduDDdvDotdRDot;
        dduDDduDotdvDot=dduDDdvDotduDot;
        dduDDdvDotdvDot=-(2*u*diffU2)/w2;

        %%%
        ddvDDdrdr=-((4*rDot*vDot)/r^3);
        ddvDDdudr=0;
        ddvDDdvdr=0;
        ddvDDdrDotdr=(2*vDot)/r2;
        ddvDDduDotdr=0;
        ddvDDdvDotdr=(2*rDot)/r2;

        ddvDDdrdu=ddvDDdudr;
        ddvDDdudu=-(2*v*(uDot2*(-1+v)*(1+v)*(-1-3*u2+v2)+2*u*uDot*v*(3+u2-3*v2)*vDot+v2*(1+3*u2-v2)*vDot2))/w6;
        ddvDDdvdu=-(2*u*uDot2*(1-v4+u2*(-1+3*v2))+4*uDot*v*(1-u4+(-1+3*u2)*v2)*vDot+2*u*v2*(3-3*u2+v2)*vDot2)/w6;
        ddvDDdrDotdu=0;
        ddvDDduDotdu=(2*v*(2*u*uDot*(-diffV2)+v*(-1-u2+v2)*vDot))/w4;
        ddvDDdvDotdu=-((2*v2*(uDot*(1+u2-v2)+2*u*v*vDot))/w4);

        ddvDDdrdv=ddvDDdvdr;
        ddvDDdudv=ddvDDdvdu;
        ddvDDdvdv=-(2*u2*uDot2*v*(3-3*u2+v2)-4*u*diffU2*uDot*(-1+u2-3*v2)*vDot+2*diffU2*v*(3-3*u2+v2)*vDot2)/w6;
        ddvDDdrDotdv=0;
        ddvDDduDotdv=(-2*uDot*(diffV4+u2*(-1+3*v2))-4*u*diffU2*v*vDot)/w4;
        ddvDDdvDotdv=(2*diffU2*(-2*u*uDot*v+u2*vDot-(1+v2)*vDot))/w4;

        ddvDDdrdRDot=ddvDDdrDotdr;
        ddvDDdudRDot=ddvDDdrDotdu;
        ddvDDdvdRDot=ddvDDdrDotdv;
        ddvDDdrDotdRDot=0;
        ddvDDduDotdRDot=0;
        ddvDDdvDotdRDot=-(2/r);

        ddvDDdrduDot=ddvDDduDotdr;
        ddvDDduduDot=ddvDDduDotdu;
        ddvDDdvduDot=ddvDDduDotdv;
        ddvDDdrDotduDot=ddvDDduDotdRDot;
        ddvDDduDotduDot=-((2*v*diffV2)/w2);
        ddvDDdvDotduDot=-(2*u*v2)/w2;

        ddvDDdrdvDot=ddvDDdvDotdr;
        ddvDDdudvDot=ddvDDdvDotdu;
        ddvDDdvdvDot=ddvDDdvDotdv;
        ddvDDdrDotdvDot=ddvDDdvDotdRDot;
        ddvDDduDotdvDot=ddvDDdvDotduDot;
        ddvDDdvDotdvDot=((2*(-1+u^2)*v)/w2);
        
        aHess(:,:,1)=[0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrdr,       ddrDDdudr,      ddrDDdvdr,      0,              ddrDDduDotdr,   ddrDDdvDotdr;
                      dduDDdrdr,       dduDDdudr,      dduDDdvdr,      dduDDdrDotdr,   dduDDduDotdr,   dduDDdvDotdr;
                      ddvDDdrdr,       ddvDDdudr,      ddvDDdvdr,      ddvDDdrDotdr,   ddvDDduDotdr,   ddvDDdvDotdr];
        aHess(:,:,2)=[0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrdu,       ddrDDdudu,      ddrDDdvdu,      0,              ddrDDduDotdu,   ddrDDdvDotdu;
                      dduDDdrdu,       dduDDdudu,      dduDDdvdu,      dduDDdrDotdu,   dduDDduDotdu,   dduDDdvDotdu;
                      ddvDDdrdu,        ddvDDdudu,     ddvDDdvdu,      ddvDDdrDotdu,   ddvDDduDotdu,   ddvDDdvDotdu];
        aHess(:,:,3)=[0,               0,              0,              0               0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrdv,       ddrDDdudv,      ddrDDdvdv,      0,              ddrDDduDotdv,   ddrDDdvDotdv;
                      dduDDdrdv,       dduDDdudv,      dduDDdvdv,      dduDDdrDotdv,   dduDDduDotdv,   dduDDdvDotdv;
                      ddvDDdrdv,       ddvDDdudv,      ddvDDdvdv,      ddvDDdrDotdv,   ddvDDduDotdv,   ddvDDdvDotdv];
        aHess(:,:,4)=[0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrdRDot,    ddrDDdudRDot,   ddrDDdvdRDot,   0,              ddrDDduDotdRDot,ddrDDdvDotdRDot;
                      dduDDdrdRDot,    dduDDdudRDot,   dduDDdvdRDot,   dduDDdrDotdRDot,dduDDduDotdRDot,dduDDdvDotdRDot;
                      ddvDDdrdRDot,    ddvDDdudRDot,   ddvDDdvdRDot,   ddvDDdrDotdRDot,ddvDDduDotdRDot,ddvDDdvDotdRDot];
        aHess(:,:,5)=[0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrduDot,    ddrDDduduDot,   ddrDDdvduDot,   0,              ddrDDduDotduDot,ddrDDdvDotduDot;
                      dduDDdrduDot,    dduDDduduDot,   dduDDdvduDot,   dduDDdrDotduDot,dduDDduDotduDot,dduDDdvDotduDot;
                      ddvDDdrduDot,    ddvDDduduDot,   ddvDDdvduDot,   ddvDDdrDotduDot,ddvDDduDotduDot,ddvDDdvDotduDot];
        aHess(:,:,6)=[0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      0,               0,              0,              0,              0,              0;
                      ddrDDdrdvDot,    ddrDDdudvDot,   ddrDDdvdvDot,   0,              ddrDDduDotdvDot,ddrDDdvDotdvDot;
                      dduDDdrdvDot,    dduDDdudvDot,   dduDDdvdvDot,   dduDDdrDotdvDot,dduDDduDotdvDot,dduDDdvDotdvDot;
                      ddvDDdrdvDot,    ddvDDdudvDot,   ddvDDdvdvDot,   ddvDDdrDotdvDot,ddvDDduDotdvDot,ddvDDdvDotdvDot];
        
        
        if(nargout>3) 
            papt=zeros(6,1);
        end
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
