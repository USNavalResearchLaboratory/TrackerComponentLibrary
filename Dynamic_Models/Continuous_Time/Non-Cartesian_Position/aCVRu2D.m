function [aVal,aJacob,aHess,papt]=aCVRu2D(x)
%%ACVRU2D The drift function for a continuous-time motion model where the
%         target state is given in 2D range and a single directionc cosine 
%         coordinates and the motion is constant velocity in 2D Cartesian
%         coordinates.
%
%INPUTS: x The 4XN state vector of N targets in the order of
%          [r;u;rDot;uDot]. it is required that abs(u)<1. The target is
%          assumed to be in front of the sensor (hence the lack of
%          explicitly providing v).
%
%OUTPUTS: aVal The 4XN set of time-derivatives of the N state vectors
%              under the Cartesian linear motion model in r-u coordinates.
%       aJacob This and higher partial derivatives can only be requested
%              if N=1. This is the 4X4 matrix of partial derivatives of
%              aVal such that aJacob(:,i) is the partial derivative of aVal
%              with respect to x(i).
%        aHess The 4X4X$ matrix of second partial derivatives of aVal such
%              that aHess(:,k1,k2) is the second partial derivative of aVal
%              with respect to xPoints(k1) and x(k2).
%         papt The 4X1 partial derivative with resect to time of aVal. This
%              is all zeros, because the model is time invariant.
%
%A full derivation of the dynamic model is provided in [1]. Summarizing it
%here, let rVec be a position vector. We shall use the orthonormal basis
%vectors:
%uVec_r=[u;sqrt(1-u2)]; and uVec_u=[sqrt(1-u2);-u];
%Note that d/du uVec_r=(1/sqrt(1-u2))*uvec_u,
%the time derivative of uVec_r is
%uVecDot_r=(uDot/sqrt(1-u2))*u_u
%and the time derivative of uVec_u is
%uVecDot_u=-uDot/sqrt(1-u2)*uVec_r
%A position can be expressed as
%rVec=r*uVec_r
%Taking the derivative, one gets the velocity vector
%rVecDot=rDot*uVec_r+r*uVecDot_r=rDot*uVec_r+r*uDot/9sqrt(1-u2)*uVec_u
%Taking another derivative, one gets expressions for the acceleration
%rVecDDot=(rDDot-r*uDot2/(1-u2))*uVec_r+((2*rDot*uDot+r*uDDot)/sqrt(1-u2)+r*u*uDot2/(1-u2)^(3/2))*uVec_u
%For a zero acceleration model, we note that the coefficient of uVec_r and
%the coefficient of uVec_u must be zero. Thus, we get the differential
%equations. 
%rDDot=r*uDot2/(1-u2)
%uDDot=(-2/r)*rDot*uDot-u*uDot2/(1-u2)
%These form the basis of the dynamic model.
%
%EXAMPLE 1:
%Here, we verify that integrating forward with this model is equivalent to
%linear motion in Cartesian coordinates.
% xInitCart=[1000;40;-5;20];
% T=10;%The prediction time.
% F=FPolyKal(T,4,1);
% xInitRu2D=stateCart2Ru2D(xInitCart);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsRu2D=RKAdaptiveOverRange(xInitRu2D,[0;T],@(x,t)aCVRu2D(x),0.1,0,[],[],RelTol,AbsTol);
% xEndRu2DRK=xStepsRu2D(:,end);
% xEndRu2DExact=stateCart2Ru2D(xEndCart);
% max(abs(xEndRu2DRK-xEndRu2DExact))
%One will observe that the error is less than 1e-9, which is a good
%agreement.
%
%EXAMPLE 2:
%Here, we verify that the Jacobian is consistent with the numerical
%derivative of aDeriv.
% x=[1000;-0.1;10;0.01];
% [aDeriv,aJacob]=aCVRu2D(x);
% AJacobNumDiff=numDiff(x,@(xState)aCVRu2D(xState),4);
% err=(aJacob-AJacobNumDiff)./AJacobNumDiff;
% max(abs(err(:)))
%One will see that the maximum error is on the order of 1e-9, indicating
%good agreement.
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
rDot=x(3,:);
uDot=x(4,:);

r2=r.*r;
u2=u.^2;
v2=max(0,1-u2);
uDot2=uDot.*uDot;

rDDot=r.*uDot2./v2;
uDDot=-2*rDot.*uDot./r-u.*uDot2./v2;

aVal=[rDot;uDot;rDDot;uDDot];

if(nargout>1)
    N=size(x,2);
    if(N>1)
        error('Derivatives are only available for numPoints=1.')
    end
    v4=v2.*v2;
    
    drDDdr=uDot2./v2;
    drDDdu=(2*r.*u.*uDot2)./v4;
    drDDduDot=(2*r.*uDot)./v2;
    
    duDDdr=(2*rDot.*uDot)./r2;
    duDDdu=-(((1+u2).*uDot2)./v4);
    duDDdrDot=-((2.*uDot)./r);
    duDDduDot=-((2.*rDot)./r)-(2*u.*uDot)./v2;
    
    aJacob=[0,           0,         1,           0;
            0,           0,         0,           1;
            drDDdr,   drDDdu,       0,           drDDduDot;
            duDDdr,   duDDdu,       duDDdrDot,   duDDduDot];
        
    if(nargout>2)
        aHess=zeros(4,4,4);
        v6=v4*v2;
        
        drDDdrdr=0;
        drDDdudr=(2*u*uDot2)/v4;
        drDDduDotdr=-((2*uDot)/v2);
        
        drDDdrdu=drDDdudr;
        drDDdudu=-((2*r*(1+3*u2)*uDot2)/v6);
        drDDduDotdu=(4*r*u*uDot)/v4;
        
        drDDdrdrDot=0;
        drDDdudrDot=0;
        drDDduDotdrDot=0;
        
        drDDdrduDot=drDDduDotdr;
        drDDduduDot=drDDduDotdu;
        drDDduDotduDot=-((2*r)/v2);

        duDDdrdr=-((4*rDot*uDot)/r^3);
        duDDdudr=0;
        duDDdrDotdr=(2*uDot)/r2;
        duDDduDotdr=(2*rDot)/r2;
        
        duDDdrdu=duDDdudr;
        duDDdudu=(2*u*(3+u2)*uDot2)/v6;
        duDDdrDotdu=0;
        duDDduDotdu=-((2*v2*uDot)/v4);
        
        duDDdrdrDot=duDDdrDotdr;
        duDDdudrDot=duDDdrDotdu;
        duDDdrDotdrDot=0;
        duDDduDotdrDot=-(2/r);
        
        duDDdrduDot=duDDduDotdr;
        duDDduduDot=duDDduDotdu;
        duDDdrDotduDot=duDDduDotdrDot;
        duDDduDotduDot=(2*u)/v2;
        
        aHess(:,:,1)=[zeros(2,4);
                        drDDdrdr,       drDDdudr,       0,                  drDDduDotdr;
                        duDDdrdr,       duDDdudr,       duDDdrDotdr,        duDDduDotdr];
        aHess(:,:,2)=[zeros(2,4);
                        drDDdrdu,       drDDdudu,       0,                  drDDduDotdu;
                        duDDdrdu,       duDDdudu,       duDDdrDotdu,        duDDduDotdu];
        aHess(:,:,3)=[zeros(2,4);
                        drDDdrdrDot,    drDDdudrDot,    0,                  drDDduDotdrDot;
                        duDDdrdrDot,    duDDdudrDot,    duDDdrDotdrDot,     duDDduDotdrDot];
        aHess(:,:,4)=[zeros(2,4);
                        drDDdrduDot,    drDDduduDot,    0,                  drDDduDotduDot;
                        duDDdrduDot,    duDDduduDot,    duDDdrDotduDot,     duDDduDotduDot];
        
        if(nargout>3)
            papt=zeros(4,1);
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
