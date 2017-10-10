function [aDeriv,aJacob]=aCVSpherical(x,t,systemType)
%%ACVSPHERICAL The drift function for a continuous-time motion model where
%          the target state is given in monostatic spherical coordinates
%          and the motion is constant velocity in 3D Cartesian coordinates.
%
%INPUTS: x The 6XN state vector of N targets in the order of
%          [r;theta;phi;rDot;thetaDot;phiDot], where theta is azimuth and
%          phi elevation. The angles are given in radians.
%        t An unused time component so that aCVSpherical can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
% systemType An optional parameter specifying the axis from which the
%          angles are measured in radians. Possible values are
%          0 (The default if omitted) Azimuth is measured 
%            counterclockwise from the x-axis in the x-y plane. Elevation
%            is measured up from the x-y plane (towards the z-axis). This
%            is consistent with common spherical coordinate systems for
%            specifying longitude (azimuth) and geocentric latitude
%            (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z-axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given
%            elevation, one is given the angle away from the z-axis, which
%            is (pi/2-elevation).
%
%OUTPUTS: aDeriv The 6XN set of time-derivatives of the N state vectors
%                under the Cartesian linear motion model in spherical
%                coordinates.
%         AJacob The 6X6XN set of Jacobians of aDeriv, which can be useful
%                in extended Kalman filters. This is the derivative of each
%                component of aDeriv (selected by row) with respect to the
%                elements of the state [r,theta,phi,rDot,thetaDot,phiDot]
%                selected by column.
%
%A derivation of the dynamic model is provided here for systemType=0. The
%derivations for the other system types are very similar. Let rVec be a
%position vector. We shall use the orthonormal basis vectors:
%u_r=[cos(theta)cos(phi);sin(theta)*cos(phi);sin(phi)]
%u_theta=[-sin(theta);cos(theta);0]
%u_phi=[-cos(theta)*sin(phi);-sin(theta)*sin(phi);cos(phi)];
%Note that
%u_theta=(1/cos(phi))*du_r/dtheta
%and 
%u_phi=du_r/dphi
%We will also note the time derivatives (dot terms):
%uDot_r=thetaDot*cos(phi)*u_theta+phiDot*u_phi
%uDot_theta=-thetaDot*cos(phi)*u_r+thetaDot*sin(phi)*u_phi
%uDot_phi=-phiDot*u_r-thetaDot*sin(phi)*u_theta
%The position vector is just
%rVec=r*u_r
%The velocity vector is the derivative of the psition vector and is
%rVecDot=rDot*u_r+r*uDot_r=rDot*u_r+r*thetaDot*cos(phi)*u_theta+r*phiDot*u_phi
%The acceleration vector is the derivative of the velocity vector and
%simplifies to
%rVecDDot=(rDDot-r*thetaDot*cos(phi)^2-r*phiDot^2)*u_r+((2*rDot*thetaDot+r*thetaDDot)*cos(phi)-2*r*thetaDot*phiDot*sin(phi))*u_theta+(2*rDot*phiDot+r*phiDDot+r*thetaDot^2*cos(phi)*sin(phi))*u_phi
%For a constant velocity model, the coefficients of the u_r, u_theta and
%u_phi vectors in the acceleration equation must be zero. This leads to the
%dynamics:
% rDDot=r*phiDot^2+r*thetaDot^2*cos(phi)^2
% thetaDDot=(1/r)*(-2*rDot*thetaDot+2*r*thetaDot*phiDot*tan(phi))
% phiDDot=(1/r)*(-2*rDot*phiDot-r*thetaDot^2*cos(phi)*sin(phi))
%The above three equations define the dynamic model.
%
%EXAMPLE 1:
%Here, we verify that integrating forward with this model is equivalent to
%linear motion in Cartesian coordinates.
% xInitCart=[1000;40;92;-100;20;64];
% T=50;%The prediction time.
% F=FPolyKal(T,xInitCart,1);
% systemType=0;
% xInitSpher=stateCart2Sphere(xInitCart,systemType);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsSphere=RKAdaptiveOverRange(xInitSpher,[0;T],@(x,t)aCVSpherical(x,t,systemType),0.1,0,[],[],RelTol,AbsTol);
% xEndSphereRK=xStepsSphere(:,end);
% xEndSphereExact=stateCart2Sphere(xEndCart,systemType);
% max(abs(xEndSphereRK-xEndSphereExact))
%One will observe that the error is less than 1e-8, which is a good
%agreement.
%
%EXAMPLE 2:
%Here, we verify that the Jacobian is consistent with the numerical
%derivative of aDeriv.
% x=[1000;0.1;-0.25;10;0.01;-0.02];
% systemType=2;
% [aDeriv,aJacob]=aCVSpherical(x,[],systemType);
% AJacobNumDiff=numDiff(x,@(xState)aCVSpherical(xState,[],systemType),6);
% err=(aJacob-AJacobNumDiff)./AJacobNumDiff;
% max(abs(err(:)))
%One will see that the maximum error is on the order of 6.8e-9, indicating
%good agreement.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(systemType))
    systemType=0; 
end

r=x(1,:);
phi=x(3,:);
rDot=x(4,:);
thetaDot=x(5,:);
phiDot=x(6,:);

cosPhi=cos(phi);
sinPhi=sin(phi);

if(systemType==0||systemType==1)
    tanPhi=sinPhi./cosPhi;
    
    rDDot=r.*phiDot.^2+r.*thetaDot.^2.*cosPhi.^2;
    thetaDDot=(1./r).*(-2*rDot.*thetaDot+2*r.*thetaDot.*phiDot.*tanPhi);
    phiDDot=(1./r).*(-2*rDot.*phiDot-r.*thetaDot.^2.*cosPhi.*sinPhi);
    
    aDeriv=[rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot];

    if(nargout>1)
        drDDdr=phiDot.^2+thetaDot.^2.*cosPhi.^2;
        drDDdphi=-2*r.*thetaDot.^2.*cosPhi.*sinPhi;
        drDDdthetaDot=2*r.*thetaDot.*cosPhi.^2;
        drDDdphiDot=2*phiDot.*r;

        dThetaDDdr=(2*rDot.*thetaDot)./r.^2;
        dThetaDDdphi=2*phiDot.*thetaDot./cosPhi.^2;
        dThetaDDdrDot=-((2*thetaDot)./r);
        dThetaDDdthetaDot=-((2*rDot)./r)+2*phiDot.*tanPhi;
        dThetaDDdphiDot=2*thetaDot.*tanPhi;

        dPhiDDdr=(2*phiDot.*rDot)./r.^2;
        dPhiDDdPhi=(r.*thetaDot.^2.*(sinPhi.^2-cosPhi.^2))./r;
        dPhiDDdrDot=-((2*phiDot)./r);
        dPhiDDdthetaDot=-2*thetaDot.*cosPhi.*sinPhi;
        dPhiDDdphiDot=-((2*rDot)./r);

        N=size(x,2);
        aJacob=zeros(6,6,N);
        for k=1:N
            aJacob(:,:,k)=[ 0,              0,  0,                  1,                  0,                      0;
                            0,              0,  0,                  0,                  1,                      0;
                            0,              0,  0,                  0,                  0,                      1;
                            drDDdr(k),      0,  drDDdphi(k),        0,                  drDDdthetaDot(k),       drDDdphiDot(k);
                            dThetaDDdr(k),  0,  dThetaDDdphi(k),    dThetaDDdrDot(k),   dThetaDDdthetaDot(k),   dThetaDDdphiDot(k);
                            dPhiDDdr(k),    0,  dPhiDDdPhi(k),      dPhiDDdrDot(k),     dPhiDDdthetaDot(k),     dPhiDDdphiDot(k)];
        end
    end
elseif(systemType==2)
    cotPhi=cosPhi./sinPhi;
    
    rDDot=r.*phiDot.^2+r.*thetaDot.^2.*sinPhi.^2;
    thetaDDot=(1./r).*(-2*rDot.*thetaDot-2*r.*thetaDot.*phiDot.*cotPhi);
    phiDDot=(1./r).*(-2*rDot.*phiDot+r.*thetaDot.^2.*cosPhi.*sinPhi);
    
    aDeriv=[rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot];
    
    if(nargout>1)
        drDDdr=phiDot.^2+thetaDot.^2.*sinPhi.^2;
        drDDdphi=2*r.*thetaDot.^2.*cosPhi.*sinPhi;
        drDDdthetaDot=2*r.*thetaDot.*sinPhi.^2;
        drDDdphiDot=2*phiDot.*r;

        dThetaDDdr=(2*rDot.*thetaDot)./r.^2;
        dThetaDDdphi=2*phiDot.*thetaDot./sinPhi.^2;
        dThetaDDdrDot=-((2*thetaDot)./r);
        dThetaDDdthetaDot=-((2*(rDot+phiDot.*r.*cotPhi))./r);
        dThetaDDdphiDot=-2*thetaDot.*cotPhi;

        dPhiDDdr=(2*phiDot.*rDot)./r.^2;
        dPhiDDdPhi=(r.*thetaDot.^2.*cosPhi.^2-r.*thetaDot.^2.*sinPhi.^2)./r;
        dPhiDDdrDot=-((2*phiDot)./r);
        dPhiDDdthetaDot=2*thetaDot.*cosPhi.*sinPhi;
        dPhiDDdphiDot=-((2*rDot)./r);

        N=size(x,2);
        aJacob=zeros(6,6,N);
        for k=1:N
            aJacob(:,:,k)=[ 0,              0,  0,                  1,                  0,                      0;
                            0,              0,  0,                  0,                  1,                      0;
                            0,              0,  0,                  0,                  0,                      1;
                            drDDdr(k),      0,  drDDdphi(k),        0,                  drDDdthetaDot(k),       drDDdphiDot(k);
                            dThetaDDdr(k),  0,  dThetaDDdphi(k),    dThetaDDdrDot(k),   dThetaDDdthetaDot(k),   dThetaDDdphiDot(k);
                            dPhiDDdr(k),    0,  dPhiDDdPhi(k),      dPhiDDdrDot(k),     dPhiDDdthetaDot(k),     dPhiDDdphiDot(k)];
        end
    end
else
    error('Unknown systemType specified.') 
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
