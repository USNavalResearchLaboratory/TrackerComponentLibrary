function [aVal,aJacob,aHess,papt]=aCVSpherical(x,systemType)
%%ACVSPHERICAL The drift function for a continuous-time motion model where
%          the target state is given in monostatic spherical coordinates
%          and the motion is constant velocity in 3D Cartesian coordinates.
%
%INPUTS: x The 6XN state vector of N targets in the order of
%          [r;theta;phi;rDot;thetaDot;phiDot], where theta is azimuth and
%          phi elevation. The angles are given in radians.
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
%OUTPUTS: aVal The 6XN set of time-derivatives of the N state vectors
%              under the Cartesian linear motion model in spherical
%              coordinates.
%       aJacob This and higher partial derivatives can only be requested
%              if N=1. This is the 6X6 matrix of partial derivatives of
%              aVal such that aJacob(:,i) is the partial derivative of aVal
%              with respect to x(i).
%        aHess The 6X6X6  matrix of second derivatives of aVal such
%              that aHess(:,k1,k2) is the second partial derivative of
%              aVal with respect to x(k1) and x(k2).
%         papt The 6X1  partial derivative with resect to time of aVal.
%              This is all zeros, because the model is time invariant.
%
%A full derivation of the dynamic model is provided in [1]. Summarizing it
%here for systemType=0, let rVec be a position vector. We shall use the
%orthonormal basis vectors:
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
%The velocity vector is the derivative of the position vector and is
%rVecDot=rDot*u_r+r*uDot_r=rDot*u_r+r*thetaDot*cos(phi)*u_theta+r*phiDot*u_phi
%The acceleration vector is the derivative of the velocity vector and
%simplifies to
%rVecDDot=(rDDot-r*thetaDot*cos(phi)^2-r*phiDot^2)*u_r+((2*rDot*thetaDot+r*thetaDDot)*cos(phi)-2*r*thetaDot*phiDot*sin(phi))*u_theta+(2*rDot*phiDot+r*phiDDot+r*thetaDot2*cos(phi)*sin(phi))*u_phi
%For a constant velocity model, the coefficients of the u_r, u_theta and
%u_phi vectors in the acceleration equation must be zero. This leads to the
%dynamics:
% rDDot=r*phiDot^2+r*thetaDot2*cos(phi)^2
% thetaDDot=(1/r)*(-2*rDot*thetaDot+2*r*thetaDot*phiDot*tan(phi))
% phiDDot=(1/r)*(-2*rDot*phiDot-r*thetaDot2*cos(phi)*sin(phi))
%The above three equations define the dynamic model.
%
%EXAMPLE 1:
%Here, we verify that integrating forward with this model is equivalent to
%linear motion in Cartesian coordinates.
% xInitCart=[1000;40;92;-100;20;64];
% T=50;%The prediction time.
% F=FPolyKal(T,6,1);
% systemType=0;
% xInitSpher=stateCart2Sphere(xInitCart,systemType);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsSphere=RKAdaptiveOverRange(xInitSpher,[0;T],@(x,t)aCVSpherical(x,systemType),0.1,0,[],[],RelTol,AbsTol);
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
% [aDeriv,aJacob]=aCVSpherical(x,systemType);
% AJacobNumDiff=numDiff(x,@(xState)aCVSpherical(xState,systemType),6);
% err=(aJacob-AJacobNumDiff)./AJacobNumDiff;
% max(abs(err(:)))
%One will see that the maximum error is on the order of 2.4940e-10,
%indicating good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic Linear Dynamic Models in Local Coordinates,"
%    Naval Research Laboratory, Washington, D.C., Tech. Rep.
%    NRL/MR/5344--19-9882, 24 Aug. 2019.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
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
    
    phiDot2=phiDot.^2;
    thetaDot2=thetaDot.^2;
    sinPhi2=sinPhi.^2;
    cosPhi2=cosPhi.^2;
    cosPhiSinPhi=cosPhi.*sinPhi;
    
    rDDot=r.*phiDot2+r.*thetaDot2.*cosPhi2;
    thetaDDot=(-2*rDot.*thetaDot+2*r.*thetaDot.*phiDot.*tanPhi)./r;
    phiDDot=(-2*rDot.*phiDot-r.*thetaDot2.*cosPhiSinPhi)./r;
    
    aVal=[rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot];

    if(nargout>1)
        N=size(x,2);
        if(N>1)
            error('Derivatives are only available for numPoints=1.')
        end
        
        r2=r*r;
        
        drDDdr=phiDot2+thetaDot2*cosPhi2;
        drDDdphi=-2*r*thetaDot2*cosPhiSinPhi;
        drDDdthetaDot=2*r*thetaDot*cosPhi2;
        drDDdphiDot=2*phiDot*r;

        dThetaDDdr=(2*rDot.*thetaDot)/r2;
        dThetaDDdphi=2*phiDot.*thetaDot/cosPhi2;
        dThetaDDdrDot=-((2*thetaDot)/r);
        dThetaDDdthetaDot=-((2*rDot)/r)+2*phiDot*tanPhi;
        dThetaDDdphiDot=2*thetaDot.*tanPhi;

        dPhiDDdr=(2*phiDot*rDot)/r2;
        dPhiDDdPhi=(r*thetaDot2*(sinPhi2-cosPhi2))/r;
        dPhiDDdrDot=-((2*phiDot)/r);
        dPhiDDdthetaDot=-2*thetaDot.*cosPhiSinPhi;
        dPhiDDdphiDot=-((2*rDot)/r);

        aJacob=[ 0,           0,  0,               1,               0,                   0;
                 0,           0,  0,               0,               1,                   0;
                 0,           0,  0,               0,               0,                   1;
                 drDDdr,      0,  drDDdphi,        0,               drDDdthetaDot,       drDDdphiDot;
                 dThetaDDdr,  0,  dThetaDDdphi,    dThetaDDdrDot,   dThetaDDdthetaDot,   dThetaDDdphiDot;
                 dPhiDDdr,    0,  dPhiDDdPhi,      dPhiDDdrDot,     dPhiDDdthetaDot,     dPhiDDdphiDot];

        if(nargout>2)
            aHess=zeros(6,6,6);
            
            r3=r2*r;
            cosPhi2=cosPhi*cosPhi;
            cos2Phi=cosPhi^2-sinPhi^2;

            drDDdrdr=0;
            drDDdphidr=-2*thetaDot2*cosPhiSinPhi;
            drDDdthetaDotdr=2*thetaDot*cosPhi2;
            drDDdphiDotdr=2*phiDot;

            drDDdrdphi=drDDdphidr;
            drDDdphidphi=-2*r*thetaDot2*cos2Phi;
            drDDdthetaDotdphi=-4*r*thetaDot*cosPhiSinPhi;
            drDDdphiDotdphi=0;

            drDDdrdrDot=0;
            drDDdphidrDot=0;
            drDDdthetaDotdrDot=0;
            drDDdphiDotdrDot=0;

            drDDdrdthetaDot=drDDdthetaDotdr;
            drDDdphidthetaDot=drDDdthetaDotdphi;
            drDDdthetaDotdthetaDot=2*r*cosPhi2;
            drDDdphiDotdthetaDot=0;

            drDDdrdphiDot=drDDdphiDotdr;
            drDDdphidphiDot=drDDdphiDotdphi;
            drDDdthetaDotdphiDot=drDDdphiDotdthetaDot;
            drDDdphiDotdphiDot=2*r;

            %%%%
            dThetaDDdrdr=-((4*rDot*thetaDot)/r3);
            dThetaDDdphidr=0;
            dThetaDDdrDotdr=(2*thetaDot)/r2;
            dThetaDDdthetaDotdr=(2*rDot)/r2;
            dThetaDDdphiDotdr=0;

            dThetaDDdrdphi=dThetaDDdphidr;
            dThetaDDdphidphi=4*tanPhi*phiDot*thetaDot/cosPhi2;
            dThetaDDdrDotdphi=0;
            dThetaDDdthetaDotdphi=2*phiDot/cosPhi2;
            dThetaDDdphiDotdphi=2*thetaDot/cosPhi2;

            dThetaDDdrdrDot=dThetaDDdrDotdr;
            dThetaDDdphidrDot=dThetaDDdrDotdphi;
            dThetaDDdrDotdrDot=0;
            dThetaDDdthetaDotdrDot=-(2/r);
            dThetaDDdphiDotdrDot=0;

            dThetaDDdrdthetaDot=dThetaDDdthetaDotdr;
            dThetaDDdphidthetaDot=dThetaDDdthetaDotdphi;
            dThetaDDdrDotdthetaDot=dThetaDDdthetaDotdrDot;
            dThetaDDdthetaDotdthetaDot=0;
            dThetaDDdphiDotdthetaDot=2*tanPhi;

            dThetaDDdrdphiDot=dThetaDDdphiDotdr;
            dThetaDDdphidphiDot=dThetaDDdphiDotdphi;
            dThetaDDdrDotdphiDot=dThetaDDdphiDotdrDot;
            dThetaDDdthetaDotdphiDot=dThetaDDdphiDotdthetaDot;
            dThetaDDdphiDotdphiDot=0;

            %%%
            dPhiDDdrdr=-((4*phiDot*rDot)/r3);
            dPhiDDdPhidr=0;
            dPhiDDdrDotdr=(2*phiDot)/r2;
            dPhiDDdthetaDotdr=0;
            dPhiDDdphiDotdr=(2*rDot)/r2;

            dPhiDDdrdphi=dPhiDDdPhidr;
            dPhiDDdPhidphi=4*thetaDot2*cosPhiSinPhi;
            dPhiDDdrDotdphi=0;
            dPhiDDdthetaDotdphi=-2*thetaDot*cos2Phi;
            dPhiDDdphiDotdphi=0;

            dPhiDDdrdrDot=dPhiDDdrDotdr;
            dPhiDDdPhidrDot=dPhiDDdrDotdphi;
            dPhiDDdrDotdrDot=0;
            dPhiDDdthetaDotdrDot=0;
            dPhiDDdphiDotdrDot=-(2/r);

            dPhiDDdrdthetaDot=dPhiDDdthetaDotdr;
            dPhiDDdPhidthetaDot=dPhiDDdthetaDotdphi;
            dPhiDDdrDotdthetaDot=dPhiDDdthetaDotdrDot;
            dPhiDDdthetaDotdthetaDot=-2*cosPhiSinPhi;
            dPhiDDdphiDotdthetaDot=0;

            dPhiDDdrdphiDot=dPhiDDdphiDotdr;
            dPhiDDdPhidphiDot=dPhiDDdphiDotdphi;
            dPhiDDdrDotdphiDot=dPhiDDdphiDotdrDot;
            dPhiDDdthetaDotdphiDot=dPhiDDdphiDotdthetaDot;
            dPhiDDdphiDotdphiDot=0;
            
            aHess(:,:,1)=[ zeros(3,6);
                           drDDdrdr,      0,  drDDdphidr,        0,                 drDDdthetaDotdr,         drDDdphiDotdr;
                           dThetaDDdrdr,  0,  dThetaDDdphidr,    dThetaDDdrDotdr,   dThetaDDdthetaDotdr,   dThetaDDdphiDotdr;
                           dPhiDDdrdr,    0,  dPhiDDdPhidr,      dPhiDDdrDotdr,     dPhiDDdthetaDotdr,     dPhiDDdphiDotdr];
            aHess(:,:,3)=[ zeros(3,6);
                           drDDdrdphi,      0,  drDDdphidphi,        0,                 drDDdthetaDotdphi,       drDDdphiDotdphi;
                           dThetaDDdrdphi,  0,  dThetaDDdphidphi,    dThetaDDdrDotdphi, dThetaDDdthetaDotdphi,   dThetaDDdphiDotdphi;
                           dPhiDDdrdphi,    0,  dPhiDDdPhidphi,      dPhiDDdrDotdphi,   dPhiDDdthetaDotdphi,     dPhiDDdphiDotdphi];
            aHess(:,:,4)=[ zeros(3,6);
                           drDDdrdrDot,      0,  drDDdphidrDot,        0,                  drDDdthetaDotdrDot,       drDDdphiDotdrDot;
                           dThetaDDdrdrDot,  0,  dThetaDDdphidrDot,    dThetaDDdrDotdrDot, dThetaDDdthetaDotdrDot,   dThetaDDdphiDotdrDot;
                           dPhiDDdrdrDot,    0,  dPhiDDdPhidrDot,      dPhiDDdrDotdrDot,   dPhiDDdthetaDotdrDot,     dPhiDDdphiDotdrDot];
            aHess(:,:,5)=[ zeros(3,6);
                           drDDdrdthetaDot,      0,  drDDdphidthetaDot,        0,                      drDDdthetaDotdthetaDot,       drDDdphiDotdthetaDot;
                           dThetaDDdrdthetaDot,  0,  dThetaDDdphidthetaDot,    dThetaDDdrDotdthetaDot, dThetaDDdthetaDotdthetaDot,   dThetaDDdphiDotdthetaDot;
                           dPhiDDdrdthetaDot,    0,  dPhiDDdPhidthetaDot,      dPhiDDdrDotdthetaDot,   dPhiDDdthetaDotdthetaDot,     dPhiDDdphiDotdthetaDot];
            aHess(:,:,6)=[ zeros(3,6);
                           drDDdrdphiDot,      0,  drDDdphidphiDot,        0,                    drDDdthetaDotdphiDot,       drDDdphiDotdphiDot;
                           dThetaDDdrdphiDot,  0,  dThetaDDdphidphiDot,    dThetaDDdrDotdphiDot, dThetaDDdthetaDotdphiDot,   dThetaDDdphiDotdphiDot;
                           dPhiDDdrdphiDot,    0,  dPhiDDdPhidphiDot,      dPhiDDdrDotdphiDot,   dPhiDDdthetaDotdphiDot,     dPhiDDdphiDotdphiDot];
            
            
            if(nargout>3)
                papt=zeros(6,1);
            end
        end
    end
elseif(systemType==2)
    cotPhi=cosPhi./sinPhi;
    phiDot2=phiDot.^2;
    thetaDot2=thetaDot.^2;
    sinPhi2=sinPhi.^2;
    cosPhiSinPhi=cosPhi.*sinPhi;

    rDDot=r.*phiDot2+r.*thetaDot2.*sinPhi2;
    thetaDDot=(1./r).*(-2*rDot.*thetaDot-2*r.*thetaDot.*phiDot.*cotPhi);
    phiDDot=(1./r).*(-2*rDot.*phiDot+r.*thetaDot2.*cosPhiSinPhi);
    
    aVal=[rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot];
    
    if(nargout>1)
        N=size(x,2);
        if(N>1)
            error('Derivatives are only available for numPoints=1.')
        end

        r2=r*r;
        cosPhi2=cosPhi*cosPhi;
        
        drDDdr=phiDot2+thetaDot2*sinPhi2;
        drDDdphi=2*r*thetaDot2*cosPhiSinPhi;
        drDDdthetaDot=2*r*thetaDot*sinPhi2;
        drDDdphiDot=2*phiDot*r;

        dThetaDDdr=(2*rDot*thetaDot)/r2;
        dThetaDDdphi=2*phiDot*thetaDot/sinPhi2;
        dThetaDDdrDot=-((2*thetaDot)/r);
        dThetaDDdthetaDot=-((2*(rDot+phiDot*r*cotPhi))/r);
        dThetaDDdphiDot=-2*thetaDot.*cotPhi;

        dPhiDDdr=(2*phiDot.*rDot)/r2;
        dPhiDDdPhi=(r*thetaDot2.*cosPhi2-r.*thetaDot2.*sinPhi2)/r;
        dPhiDDdrDot=-((2*phiDot)/r);
        dPhiDDdthetaDot=2*thetaDot*cosPhiSinPhi;
        dPhiDDdphiDot=-((2*rDot)/r);

        aJacob=[ 0,           0,  0,               1,               0,                   0;
                 0,           0,  0,               0,               1,                   0;
                 0,           0,  0,               0,               0,                   1;
                 drDDdr,      0,  drDDdphi,        0,               drDDdthetaDot,       drDDdphiDot;
                 dThetaDDdr,  0,  dThetaDDdphi,    dThetaDDdrDot,   dThetaDDdthetaDot,   dThetaDDdphiDot;
                 dPhiDDdr,    0,  dPhiDDdPhi,      dPhiDDdrDot,     dPhiDDdthetaDot,     dPhiDDdphiDot];

        if(nargout>2)
            aHess=zeros(6,6,6);
 
            sin2Phi=2*sinPhi*cosPhi;
            cos2Phi=cosPhi^2-sinPhi^2;
            r3=r2*r;

            drDDdrdr=0;
            drDDdphidr=thetaDot2*sin2Phi;
            drDDdthetaDotdr=2*thetaDot*sinPhi2;
            drDDdphiDotdr=2*phiDot;

            drDDdrdphi=drDDdphidr;
            drDDdphidphi=2*r*thetaDot2*cos2Phi;
            drDDdthetaDotdphi=4*r*thetaDot*cosPhiSinPhi;
            drDDdphiDotdphi=0;

            drDDdrdrDot=0;
            drDDdphidrDot=0;
            drDDdthetaDotdrDot=0;
            drDDdphiDotdrDot=0;

            drDDdrdthetaDot=drDDdthetaDotdr;
            drDDdphidthetaDot=drDDdthetaDotdphi;
            drDDdthetaDotdthetaDot=2*r*sinPhi2;
            drDDdphiDotdthetaDot=0;

            drDDdrdphiDot=drDDdphiDotdr;
            drDDdphidphiDot=drDDdphiDotdphi;
            drDDdthetaDotdphiDot=drDDdphiDotdthetaDot;
            drDDdphiDotdphiDot=2*r;

            %%%

            dThetaDDdrdr=-((4*rDot*thetaDot)/r3);
            dThetaDDdphidr=0;
            dThetaDDdrDotdr=(2*thetaDot)/r2;
            dThetaDDdthetaDotdr=(2*rDot)/r2;
            dThetaDDdphiDotdr=0;

            dThetaDDdrdphi=dThetaDDdphidr;
            dThetaDDdphidphi=-4*phiDot*thetaDot*cotPhi/sinPhi2;
            dThetaDDdrDotdphi=0;
            dThetaDDdthetaDotdphi=2*phiDot/sinPhi2;
            dThetaDDdphiDotdphi=2*thetaDot/sinPhi2;

            dThetaDDdrdrDot=dThetaDDdrDotdr;
            dThetaDDdphidrDot=dThetaDDdrDotdphi;
            dThetaDDdrDotdrDot=0;
            dThetaDDdthetaDotdrDot=-(2/r);
            dThetaDDdphiDotdrDot=0;

            dThetaDDdrdthetaDot=dThetaDDdthetaDotdr;
            dThetaDDdphidthetaDot=dThetaDDdthetaDotdphi;
            dThetaDDdrDotdthetaDot=dThetaDDdthetaDotdrDot;
            dThetaDDdthetaDotdthetaDot=0;
            dThetaDDdphiDotdthetaDot=-2*cotPhi;

            dThetaDDdrdphiDot=dThetaDDdphiDotdr;
            dThetaDDdphidphiDot=dThetaDDdphiDotdphi;
            dThetaDDdrDotdphiDot=dThetaDDdphiDotdrDot;
            dThetaDDdthetaDotdphiDot=dThetaDDdphiDotdthetaDot;
            dThetaDDdphiDotdphiDot=0;

            %%%
            dPhiDDdrdr=-((4*phiDot*rDot)/r3);
            dPhiDDdPhidr=0;
            dPhiDDdrDotdr=(2*phiDot)/r2;
            dPhiDDdthetaDotdr=0;
            dPhiDDdphiDotdr=(2*rDot)/r2;

            dPhiDDdrdphi=dPhiDDdPhidr;
            dPhiDDdPhidphi=-4*thetaDot2*cosPhiSinPhi;
            dPhiDDdrDotdphi=0;
            dPhiDDdthetaDotdphi=2*thetaDot*cos2Phi;
            dPhiDDdphiDotdphi=0;

            dPhiDDdrdrDot=dPhiDDdrDotdr;
            dPhiDDdPhidrDot=dPhiDDdrDotdphi;
            dPhiDDdrDotdrDot=0;
            dPhiDDdthetaDotdrDot=0;
            dPhiDDdphiDotdrDot=0;

            dPhiDDdrdthetaDot=dPhiDDdthetaDotdr;
            dPhiDDdPhidthetaDot=dPhiDDdthetaDotdphi;
            dPhiDDdrDotdthetaDot=dPhiDDdthetaDotdrDot;
            dPhiDDdthetaDotdthetaDot=sin2Phi;
            dPhiDDdphiDotdthetaDot=0;

            dPhiDDdrdphiDot=dPhiDDdphiDotdr;
            dPhiDDdPhidphiDot=dPhiDDdphiDotdphi;
            dPhiDDdrDotdphiDot=dPhiDDdphiDotdrDot;
            dPhiDDdthetaDotdphiDot=dPhiDDdphiDotdthetaDot;
            dPhiDDdphiDotdphiDot=0;
            
            aHess(:,:,1)=[ zeros(3,6);
                           drDDdrdr,      0,  drDDdphidr,        0,               drDDdthetaDotdr,       drDDdphiDotdr;
                           dThetaDDdrdr,  0,  dThetaDDdphidr,    dThetaDDdrDotdr, dThetaDDdthetaDotdr,   dThetaDDdphiDotdr;
                           dPhiDDdrdr,    0,  dPhiDDdPhidr,      dPhiDDdrDotdr,   dPhiDDdthetaDotdr,     dPhiDDdphiDotdr];
            aHess(:,:,3)=[ zeros(3,6);
                           drDDdrdphi,      0,  drDDdphidphi,        0,                 drDDdthetaDotdphi,       drDDdphiDotdphi;
                           dThetaDDdrdphi,  0,  dThetaDDdphidphi,    dThetaDDdrDotdphi, dThetaDDdthetaDotdphi,   dThetaDDdphiDotdphi;
                           dPhiDDdrdphi,    0,  dPhiDDdPhidphi,      dPhiDDdrDotdphi,   dPhiDDdthetaDotdphi,     dPhiDDdphiDotdphi];
            aHess(:,:,4)=[ zeros(3,6);
                           drDDdrdrDot,      0,  drDDdphidrDot,        0,                  drDDdthetaDotdrDot,       drDDdphiDotdrDot;
                           dThetaDDdrdrDot,  0,  dThetaDDdphidrDot,    dThetaDDdrDotdrDot, dThetaDDdthetaDotdrDot,   dThetaDDdphiDotdrDot;
                           dPhiDDdrdrDot,    0,  dPhiDDdPhidrDot,      dPhiDDdrDotdrDot,   dPhiDDdthetaDotdrDot,     dPhiDDdphiDotdrDot];
            aHess(:,:,5)=[ zeros(3,6);
                           drDDdrdthetaDot,      0,  drDDdphidthetaDot,        0,                      drDDdthetaDotdthetaDot,       drDDdphiDotdthetaDot;
                           dThetaDDdrdthetaDot,  0,  dThetaDDdphidthetaDot,    dThetaDDdrDotdthetaDot, dThetaDDdthetaDotdthetaDot,   dThetaDDdphiDotdthetaDot;
                           dPhiDDdrdthetaDot,    0,  dPhiDDdPhidthetaDot,      dPhiDDdrDotdthetaDot,    dPhiDDdthetaDotdthetaDot,    dPhiDDdphiDotdthetaDot];
            aHess(:,:,6)=[ zeros(3,6);
                           drDDdrdphiDot,      0,  drDDdphidphiDot,        0,                    drDDdthetaDotdphiDot,       drDDdphiDotdphiDot;
                           dThetaDDdrdphiDot,  0,  dThetaDDdphidphiDot,    dThetaDDdrDotdphiDot, dThetaDDdthetaDotdphiDot,   dThetaDDdphiDotdphiDot;
                           dPhiDDdrdphiDot,    0,  dPhiDDdPhidphiDot,      dPhiDDdrDotdphiDot,   dPhiDDdthetaDotdphiDot,     dPhiDDdphiDotdphiDot];
            
            
            if(nargout>3)
                papt=zeros(6,1);
            end
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
