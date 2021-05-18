function [aVals,aJacob,aHess,papt]=aSpiral(xPoints)
%%ASPIRAL The drift function for a spiraling target motion model in 3
%         dimensions. This formulation of the dynamic model is more
%         difficult to use to make a target go in a desired direction than
%         the function aSpiralSimp.
%
%INPUTS: xPoints The 10XN target state at time t. It consists of position,
%                velocity, acceleration, and a scalar spiral rate. If x is
%                an xDim X numStates matrix, then the spiraling model is
%                evaluated for all of the state vectors.
%
%OUTPUTS: aVals The 10XN flat-Earth time-derivative(s) of the state.
%        aJacob This and higher partial derivatives can only be requested
%               if N=1. This is the 10X10 matrix of partial derivatives of
%               aVals such that aJacob(:,i) is the partial derivative of
%               aVals with respect to xPoints(i).
%         aHess The 10X10X10 matrix of second partial derivatives of
%               aVals such that aHess(:,k1,k2) is the second partial
%               derivative of aVals with respect to xPoints(k1) and
%               xPoints(k2).
%          papt The 10X1 partial derivative with resect to time of aVals.
%               This is all zeros, because the model is time invariant.
%
%A derivation of the flat-Earth spiraling dynamic model is given in [1],
%where various orthogonality assumptions are present to keep the speed of
%the target from increasing or decreasing.
%
%As noted in [1], for the model to be a spiral, it is required that
%dot(a,v)=0, where a=xPoints(7:9) is the acceleration vector and
%v=xPoints(3:6) is the velocity vector.
%
%EXAMPLE:
%Here, we demonstrate that the Jacobian is consistent with numerically
%differentiating the function value.
% x=[0;%x
%    10e3;%y
%    20e3;%z
%    0;%vx
%    0;%vy
%    -2e3;%vz
%    3e3;%ax
%    0;%ay
%    0;%az
%    2.25];%omega
% [~,aJacob]=aSpiral(x);
% f=@(x)aSpiral(x);
% epsVal=1e-5*[20e3;20e3;20e3;2e2;2e2;2e2;3e3;3e3;3e3;2.25];
% aJacobNumDiff=numDiff(x,f,10,1,epsVal);
% max(max(abs((aJacob-aJacobNumDiff)./aJacob)))
%The maximum error should be about 6.5753e-11.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, Part II, pp. 4-41, Feb. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(xPoints,2);
    aVals=zeros(10,numPoints);
    
    aVals(1:6,:)=xPoints(4:9,:);
    
    rdot=xPoints(4:6,:);
    rddot=xPoints(7:9,:);
    omega=xPoints(10,:);
    
    normRDot2=xPoints(4,:).^2+xPoints(5,:).^2+xPoints(6,:).^2;
    normRDot=sqrt(normRDot2);
    normRDDotSqrd=xPoints(7,:).^2+xPoints(8,:).^2+xPoints(9,:).^2;

    aVals(7:9,:)=bsxfun(@rdivide,bsxfun(@times,bsxfun(@cross,rddot,rdot),omega),normRDot)-bsxfun(@times,rdot,normRDDotSqrd./(normRDot.^2));

    %Deal with errors that arise if the velocity is near zero, so 0/0
    %occurs. This might be more numericallty stable if when close to 0, one
    %were to use a Taylor series expansion.
    aVals(~isfinite(aVals))=0;
    
    if(nargout>1)
        if(numPoints>1)
            error('Derivatives are only available for N=1.')
        end
        
        vx=xPoints(4);
        vy=xPoints(5);
        vz=xPoints(6);
        ax=xPoints(7);
        ay=xPoints(8);
        az=xPoints(9);
        
        vx2=vx*vx;
        vy2=vy*vy;
        vz2=vz*vz;

        normRDot3=normRDot2*normRDot;
        normRDot4=normRDot2*normRDot2;
        temp0=(az*vy-ay*vz);
        temp=omega*temp0*normRDot;
                
        dXvx=(2*normRDDotSqrd*vx2+vx*temp-normRDDotSqrd*normRDot2)/normRDot4;
        dXvy=(2*normRDDotSqrd*vx*vy+vy*temp-az*omega*normRDot3)/normRDot4;
        dXvz=(2*normRDDotSqrd*vx*vz+vz*temp+ay*omega*normRDot3)/normRDot4;

        dXax=-((2*ax*vx)/normRDot2);
        dXay=-((2*ay*vx)/normRDot2)+(omega*vz)/normRDot;
        dXaz=-((2*az*vx)/normRDot2)-(omega*vy)/normRDot;
        
        dXOmega=-temp0/normRDot;
        
        temp0=(-az*vx+ax*vz);
        temp=omega*temp0*normRDot;

        dYvx=(2*normRDDotSqrd*vx*vy+vx*temp+az*omega*normRDot3)/normRDot4;
        dYvy=-((-2*normRDDotSqrd*vy2-vy*temp+normRDDotSqrd*normRDot2)/normRDot4);
        dYvz=(2*normRDDotSqrd*vy*vz+vz*temp-ax*omega*normRDot3)/normRDot4;

        dYax=-((2*ax*vy)/normRDot2)-(omega*vz)/normRDot;
        dYay=-((2*ay*vy)/normRDot2);
        dYaz=-((2*az*vy)/normRDot2)+(omega*vx)/normRDot;
        
        dYOmega=-temp0/normRDot;
        
        temp0=(ay*vx-ax*vy);
        temp=omega*temp0*normRDot;
        
        dZvx=(2*normRDDotSqrd*vx*vz+vx*temp-ay*omega*normRDot3)/normRDot4;
        dZvy=(2*normRDDotSqrd*vy*vz+vy*temp+ax*omega*normRDot3)/normRDot4;
        dZvz=(2*normRDDotSqrd*vz2+vz*temp-normRDDotSqrd*normRDot2)/normRDot4;
        
        dZax=-((2*ax*vz)/normRDot2)+(omega*vy)/normRDot;
        dZay=-((2*ay*vz)/normRDot2)-(omega*vx)/normRDot;
        dZaz=-((2*az*vz)/normRDot2);
        
        dZOmega=-temp0/normRDot;

        aJacob=[0,0,0,1,    0,      0,      0,      0,      0,      0;
                0,0,0,0,    1,      0,      0,      0,      0,      0;
                0,0,0,0,    0,      1,      0,      0,      0,      0;
                0,0,0,0,    0,      0,      1,      0,      0,      0;
                0,0,0,0,    0,      0,      0,      1,      0,      0;
                0,0,0,0,    0,      0,      0,      0,      1,      0;
                0,0,0,dXvx, dXvy,   dXvz,   dXax,   dXay,   dXaz,   dXOmega;
                0,0,0,dYvx, dYvy,   dYvz,   dYax,   dYay,   dYaz,   dYOmega;
                0,0,0,dZvx, dZvy,   dZvz,   dZax,   dZay,   dZaz,   dZOmega;
                0,0,0,0,    0,      0,      0,      0,      0,      0];
        if(nargout>2)
            aHess=zeros(10,10,10);

            normRDot6=normRDot3*normRDot3;

            vx3=vx2*vx;
            vy3=vy2*vy;
            vz3=vz2*vz;

            dXvxvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx3+3*omega*vx2*(-az*vy+ay*vz)*normRDot+6*normRDDotSqrd*vx*normRDot2+omega*(az*vy-ay*vz)*normRDot3);
            dXvyvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx2*vy+3*omega*vx*vy*(-az*vy+ay*vz)*normRDot+2*normRDDotSqrd*vy*normRDot2+az*omega*vx*normRDot3);
            dXvzvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx2*vz+3*omega*vx*vz*(-az*vy+ay*vz)*normRDot+2*normRDDotSqrd*vz*normRDot2-ay*omega*vx*normRDot3);
            dXaxvx=(2*ax*(vx2-vy2-vz2))/normRDot4;
            dXayvx=(2*ay*(vx2-vy2-vz2)-omega*vx*vz*normRDot)/normRDot4;
            dXazvx=(2*az*(vx2-vy2-vz2)+omega*vx*vy*normRDot)/normRDot4;
            dXOmegavx=(vx*(az*vy-ay*vz))/normRDot3;

            dXvxvy=dXvyvx;
            dXvyvy=(1/(normRDot6))*(-8*normRDDotSqrd*vx*vy2+3*omega*vy2*(-az*vy+ay*vz)*normRDot+2*normRDDotSqrd*vx*normRDot2+2*az*omega*vy*normRDot3+omega*(az*vy-ay*vz)*normRDot3);
            dXvzvy=(-8*normRDDotSqrd*vx*vy*vz+3*omega*vy*vz*(-az*vy+ay*vz)*normRDot-ay*omega*vy*normRDot3+az*omega*vz*normRDot3)/normRDot6;
            dXaxvy=(4*ax*vx*vy)/normRDot4;
            dXayvy=(vy*(4*ay*vx-omega*vz*normRDot))/normRDot4;
            dXazvy=(4*az*vx*vy)/normRDot4-(omega*(vx2+vz2))/normRDot3;
            dXOmegavy=(-ay*vy*vz-az*(vx2+vz2))/normRDot3;

            dXvxvz=dXvzvx;
            dXvyvz=dXvzvy;
            dXvzvz=(1/(normRDot6))*(-8*normRDDotSqrd*vx*vz2+3*omega*vz2*(-az*vy+ay*vz)*normRDot+2*normRDDotSqrd*vx*normRDot2-2*ay*omega*vz*normRDot3+omega*(az*vy-ay*vz)*normRDot3);
            dXaxvz=(4*ax*vx*vz)/normRDot4;
            dXayvz=(4*ay*vx*vz)/normRDot4+(omega*(vx2+vy2))/normRDot3;
            dXazvz=(vz*(4*az*vx+omega*vy*normRDot))/normRDot4;
            dXOmegavz=(ay*(vx2+vy2)+az*vy*vz)/normRDot3;

            dXvxax=dXaxvx;
            dXvyax=dXaxvy;
            dXvzax=dXaxvz;
            dXaxax=-((2*vx)/normRDot2);
            dXayax=0;
            dXazax=0;
            dXOmegaax=0;

            dXvxay=dXayvx;
            dXvyay=dXayvy;
            dXvzay=dXayvz;
            dXaxay=dXayax;
            dXayay=-((2*vx)/normRDot2);
            dXazay=0;
            dXOmegaay=vz/normRDot;

            dXvxaz=dXazvx;
            dXvyaz=dXazvy;
            dXvzaz=dXazvz;
            dXaxaz=dXazax;
            dXayaz=dXazay;
            dXazaz=-((2*vx)/normRDot2);
            dXOmegaaz=-(vy/normRDot);

            dXvxOmega=dXOmegavx;
            dXvyOmega=dXOmegavy;
            dXvzOmega=dXOmegavz;
            dXaxOmega=dXOmegaax;
            dXayOmega=dXOmegaay;
            dXazOmega=dXOmegaaz;
            dXOmegaOmega=0;

            %%%
            dYvxvx=(1/(normRDot6))*(2*ax^2*vy*(-3*vx2+vy2+vz2)+2*ay^2*vy*(-3*vx2+vy2+vz2)+ax*omega*vz*(-2*vx2+vy2+vz2)*normRDot+az*(2*az*vy*(-3*vx2+vy2+vz2)-3*omega*vx*(vy2+vz2)*normRDot));
            dYvyvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx*vy2+3*omega*vx*vy*(az*vx-ax*vz)*normRDot+2*normRDDotSqrd*vx*normRDot2-az*omega*vy*normRDot3);
            dYvzvx=(-8*normRDDotSqrd*vx*vy*vz+3*omega*vx*vz*(az*vx-ax*vz)*normRDot+ax*omega*vx*normRDot3-az*omega*vz*normRDot3)/normRDot6;
            dYaxvx=(vx*(4*ax*vy+omega*vz*normRDot))/normRDot4;
            dYayvx=(4*ay*vx*vy)/normRDot4;
            dYazvx=(4*az*vx*vy)/normRDot4+(omega*(vy2+vz2))/normRDot3;
            dYOmegavx=(ax*vx*vz+az*(vy2+vz2))/normRDot3;

            dYvxvy=dYvyvx;
            dYvyvy=(1/(normRDot6))*(-8*normRDDotSqrd*vy3+3*omega*vy2*(az*vx-ax*vz)*normRDot+6*normRDDotSqrd*vy*normRDot2-omega*(az*vx-ax*vz)*normRDot3);
            dYvzvy=(1/(normRDot6))*(-8*normRDDotSqrd*vy2*vz+3*omega*vy*vz*(az*vx-ax*vz)*normRDot+2*normRDDotSqrd*vz*normRDot2+ax*omega*vy*normRDot3);
            dYaxvy=(-2*ax*(vx2-vy2+vz2)+omega*vy*vz*normRDot)/normRDot4;
            dYayvy=-((2*ay*(vx2-vy2+vz2))/normRDot4);
            dYazvy=-((2*az*(vx2-vy2+vz2)+omega*vx*vy*normRDot)/normRDot4);
            dYOmegavy=(-az*vx*vy+ax*vy*vz)/normRDot3;

            dYvxvz=dYvzvx;
            dYvyvz=dYvzvy;
            dYvzvz=(1/(normRDot6))*(-8*normRDDotSqrd*vy*vz2+3*omega*vz2*(az*vx-ax*vz)*normRDot+2*normRDDotSqrd*vy*normRDot2+2*ax*omega*vz*normRDot3-omega*(az*vx-ax*vz)*normRDot3);
            dYaxvz=(4*ax*vy*vz)/normRDot4-(omega*(vx2+vy2))/normRDot3;
            dYayvz=(4*ay*vy*vz)/normRDot4;
            dYazvz=(vz*(4*az*vy-omega*vx*normRDot))/normRDot4;
            dYOmegavz=(-ax*(vx2+vy2)-az*vx*vz)/normRDot3;

            dYvxax=dYaxvx;
            dYvyax=dYaxvy;
            dYvzax=dYaxvz;
            dYaxax=-((2*vy)/normRDot2);
            dYayax=0;
            dYazax=0;
            dYOmegaax=-(vz/normRDot);

            dYvxay=dYayvx;
            dYvyay=dYayvy;
            dYvzay=dYayvz;
            dYaxay=dYayax;
            dYayay=-((2*vy)/normRDot2);
            dYazay=0;
            dYOmegaay=0;

            dYvxaz=dYazvx;
            dYvyaz=dYazvy;
            dYvzaz=dYazvz;
            dYaxaz=dYazax;
            dYayaz=dYazay;
            dYazaz=-((2*vy)/normRDot2);
            dYOmegaaz=vx/normRDot;

            dYvxOmega=dYOmegavx;
            dYvyOmega=dYOmegavy;
            dYvzOmega=dYOmegavz;
            dYaxOmega=dYOmegaax;
            dYayOmega=dYOmegaay;
            dYazOmega=dYOmegaaz;
            dYOmegaOmega=0;

            %%%
            dZvxvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx2*vz+3*omega*vx2*(-ay*vx+ax*vy)*normRDot+2*normRDDotSqrd*vz*normRDot2+2*ay*omega*vx*normRDot3+omega*(ay*vx-ax*vy)*normRDot3);
            dZvyvx=(-8*normRDDotSqrd*vx*vy*vz+3*omega*vx*vy*(-ay*vx+ax*vy)*normRDot-ax*omega*vx*normRDot3+ay*omega*vy*normRDot3)/normRDot6;
            dZvzvx=(1/(normRDot6))*(-8*normRDDotSqrd*vx*vz2+3*omega*vx*(-ay*vx+ax*vy)*vz*normRDot+2*normRDDotSqrd*vx*normRDot2+ay*omega*vz*normRDot3);
            dZaxvx=(vx*(4*ax*vz-omega*vy*normRDot))/normRDot4;
            dZayvx=(4*ay*vx*vz)/normRDot4-(omega*(vy2+vz2))/normRDot3;
            dZazvx=(4*az*vx*vz)/normRDot4;
            dZOmegavx=(-ax*vx*vy-ay*(vy2+vz2))/normRDot3;

            dZvxvy=dZvyvx;
            dZvyvy=(1/(normRDot6))*(-8*normRDDotSqrd*vy2*vz+3*omega*vy2*(-ay*vx+ax*vy)*normRDot+2*normRDDotSqrd*vz*normRDot2-2*ax*omega*vy*normRDot3+omega*(ay*vx-ax*vy)*normRDot3);
            dZvzvy=(1/(normRDot6))*(-8*normRDDotSqrd*vy*vz2+3*omega*vy*(-ay*vx+ax*vy)*vz*normRDot+2*normRDDotSqrd*vy*normRDot2-ax*omega*vz*normRDot3);
            dZaxvy=(4*ax*vy*vz)/normRDot4+(omega*(vx2+vz2))/normRDot3;
            dZayvy=(vy*(4*ay*vz+omega*vx*normRDot))/normRDot4;
            dZazvy=(4*az*vy*vz)/normRDot4;
            dZOmegavy=(ay*vx*vy+ax*(vx2+vz2))/normRDot3;

            dZvxvz=dZvzvx;
            dZvyvz=dZvzvy;
            dZvzvz=(1/(normRDot6))*(-8*normRDDotSqrd*vz3+3*omega*(-ay*vx+ax*vy)*vz2*normRDot+6*normRDDotSqrd*vz*normRDot2+omega*(ay*vx-ax*vy)*normRDot3);
            dZaxvz=-((2*ax*(vx2+vy2-vz2)+omega*vy*vz*normRDot)/normRDot4);
            dZayvz=(-2*ay*(vx2+vy2-vz2)+omega*vx*vz*normRDot)/normRDot4;
            dZazvz=-((2*az*(vx2+vy2-vz2))/normRDot4);
            dZOmegavz=((ay*vx-ax*vy)*vz)/normRDot3;

            dZvxax=dZaxvx;
            dZvyax=dZaxvy;
            dZvzax=dZaxvz;
            dZaxax=-((2*vz)/normRDot2);
            dZayax=0;
            dZazax=0;
            dZOmegaax=vy/normRDot;

            dZvxay=dZayvx;
            dZvyay=dZayvy;
            dZvzay=dZayvz;
            dZaxay=dZayax;
            dZayay=-((2*vz)/normRDot2);
            dZazay=0;
            dZOmegaay=-(vx/normRDot);

            dZvxaz=dZazvx;
            dZvyaz=dZazvy;
            dZvzaz=dZazvz;
            dZaxaz=dZazax;
            dZayaz=dZazay;
            dZazaz=-((2*vz)/normRDot2);
            dZOmegaaz=0;

            dZvxOmega=dZOmegavx;
            dZvyOmega=dZOmegavy;
            dZvzOmega=dZOmegavz;
            dZaxOmega=dZOmegaax;
            dZayOmega=dZOmegaay;
            dZazOmega=dZOmegaaz;
            dZOmegaOmega=0;
            
            aHess(:,:,4)=[zeros(6,10);
                          0,0,0,dXvxvx, dXvyvx,   dXvzvx,   dXaxvx,   dXayvx,   dXazvx,   dXOmegavx;
                          0,0,0,dYvxvx, dYvyvx,   dYvzvx,   dYaxvx,   dYayvx,   dYazvx,   dYOmegavx;
                          0,0,0,dZvxvx, dZvyvx,   dZvzvx,   dZaxvx,   dZayvx,   dZazvx,   dZOmegavx;
                          zeros(1,10)];
            aHess(:,:,5)=[zeros(6,10);
                          0,0,0,dXvxvy, dXvyvy,   dXvzvy,   dXaxvy,   dXayvy,   dXazvy,   dXOmegavy;
                          0,0,0,dYvxvy, dYvyvy,   dYvzvy,   dYaxvy,   dYayvy,   dYazvy,   dYOmegavy;
                          0,0,0,dZvxvy, dZvyvy,   dZvzvy,   dZaxvy,   dZayvy,   dZazvy,   dZOmegavy;
                          zeros(1,10)];
            aHess(:,:,6)=[zeros(6,10);
                          0,0,0,dXvxvz, dXvyvz,   dXvzvz,   dXaxvz,   dXayvz,   dXazvz,   dXOmegavz;
                          0,0,0,dYvxvz, dYvyvz,   dYvzvz,   dYaxvz,   dYayvz,   dYazvz,   dYOmegavz;
                          0,0,0,dZvxvz, dZvyvz,   dZvzvz,   dZaxvz,   dZayvz,   dZazvz,   dZOmegavz;
                          zeros(1,10)];
            aHess(:,:,7)=[zeros(6,10);
                          0,0,0,dXvxax, dXvyax,   dXvzax,   dXaxax,   dXayax,   dXazax,   dXOmegaax;
                          0,0,0,dYvxax, dYvyax,   dYvzax,   dYaxax,   dYayax,   dYazax,   dYOmegaax;
                          0,0,0,dZvxax, dZvyax,   dZvzax,   dZaxax,   dZayax,   dZazax,   dZOmegaax;
                          zeros(1,10)];
            aHess(:,:,8)=[zeros(6,10);
                          0,0,0,dXvxay, dXvyay,   dXvzay,   dXaxay,   dXayay,   dXazay,   dXOmegaay;
                          0,0,0,dYvxay, dYvyay,   dYvzay,   dYaxay,   dYayay,   dYazay,   dYOmegaay;
                          0,0,0,dZvxay, dZvyay,   dZvzay,   dZaxay,   dZayay,   dZazay,   dZOmegaay;
                          zeros(1,10)];
            aHess(:,:,9)=[zeros(6,10);
                          0,0,0,dXvxaz, dXvyaz,   dXvzaz,   dXaxaz,   dXayaz,   dXazaz,   dXOmegaaz;
                          0,0,0,dYvxaz, dYvyaz,   dYvzaz,   dYaxaz,   dYayaz,   dYazaz,   dYOmegaaz;
                          0,0,0,dZvxaz, dZvyaz,   dZvzaz,   dZaxaz,   dZayaz,   dZazaz,   dZOmegaaz;
                          zeros(1,10)];
            aHess(:,:,10)=[zeros(6,10);
                          0,0,0,dXvxOmega, dXvyOmega,   dXvzOmega,   dXaxOmega,   dXayOmega,   dXazOmega,   dXOmegaOmega;
                          0,0,0,dYvxOmega, dYvyOmega,   dYvzOmega,   dYaxOmega,   dYayOmega,   dYazOmega,   dYOmegaOmega;
                          0,0,0,dZvxOmega, dZvyOmega,   dZvzOmega,   dZaxOmega,   dZayOmega,   dZazOmega,   dZOmegaOmega;
                          zeros(1,10)];
            
            if(nargout>3)
                papt=zeros(10,1);
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
