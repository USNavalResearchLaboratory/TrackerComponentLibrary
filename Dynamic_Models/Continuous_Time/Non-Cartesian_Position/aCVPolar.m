function [aVal,aJacob,aHess,papt]=aCVPolar(x)
%%ACVPOLAR The drift function for a continuous-time motion model where the
%          target state is given in polar coordinates and the motion is
%          constant velocity in 2D Cartesian coordinates. This function
%          works regardless of whether the angle is measured
%          counterclockwise from the x axis or clockwise form the y axis.
%
%INPUTS: x The 4XN state vector of N targets in the order of
%          [r;theta;rDot;thetaDot]. The angle theta is given in radians.
%
%OUTPUTS: aVal The 4XN set of time-derivatives of the N state vectors under
%              the Cartesian linear motion model in polar coordinates.
%       aJacob This and higher partial derivatives can only be requested
%              if N=1. This is the 4X4  matrix of partial derivatives
%              of aVal such that aJacob(:,i) is the partial derivative of
%              aVal with respect to x(i).
%        aHess The 4X4X4  matrix of second derivatives of aVal such
%              that aHess(:,k1,k2) is the second partial derivative of
%              aVal with respect to x(k1) and x(k2).
%         papt The 4X1  partial derivative with resect to time of aVal.
%              This is all zeros, because the model is time invariant.
%
%A full derivation of the dynamic model is provided in [1]. Summarizing it
%here, let rVec be a position vector. We shall use the orthonormal basis
%vectors
%u_r=[cos(theta);sin(theta)] and u_theta=[-sin(theta);cos(theta)]
%if one is measuring theta from the x axis counter clockwise and 
%u_r=[sin(theta);cos(theta)] and u_theta=[cos(theta);-sin(theta)]
%if one is measuring theta from the y axis clockwise.
%Note that u_theta is the derivative of u_r with respect to theta and the
%two vectors are orthonormal.
%A position can be expressed
%rVec=r*u_r
%Note the identities for the derivatives (dot terms):
%uDot_r=thetaDot*du_r/dtheta=dTheta*u_theta
%and uDot_theta=-thetaDor*u_r
%The velocity vector is
%rVecDot=rDot*u_r+r*uDot_r=rDot*u_r+r*thetaDot*u_theta
%Taking one more derivative to get an acceleration vectors, one has, after
%simplification:
%rDDot=(rDDot-r*thetaDot^2)*u_r+(r*thetaDDot+2*rDot*thetaDot)*u_theta
%where the coefficient of u_r is the radial acceleration term and the
%coefficient of u_theta is the angular acceleration term. Noting that for
%linear motion both accelerations must be zero, one gets the following two
%equations for the second derivative terms:
%rDDot=r*thetaDot^2
%thetaDDot=-(2/r)*rDot*thetaDot
%Thus, one has a linear dynamic model in 2D polar coordinates.
%
%EXAMPLE 1:
%Here, we verify that integrating forward with this model is equivalent to
%linear motion in Cartesian coordinates.
% xInitCart=[1000;40;-100;20];
% T=50;%The prediction time.
% F=FPolyKal(T,4,1);
% systemType=0;
% xInitSpher=stateCart2Pol(xInitCart,systemType);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsPolar=RKAdaptiveOverRange(xInitSpher,[0;T],@(x,t)aCVPolar(x),0.1,0,[],[],RelTol,AbsTol);
% xEndPolarRK=xStepsPolar(:,end);
% xEndPolarExact=stateCart2Pol(xEndCart,systemType);
% max(abs(xEndPolarRK-xEndPolarExact))
%One will observe that the error is less than 1e-7, which is a good
%agreement.
%
%EXAMPLE 2:
%Here, we verify that the Jacobian is consistent with the numerical
%derivative of aDeriv.
% x=[1000;-0.1;10;0.01];
% [aDeriv,aJacob]=aCVPolar(x);
% AJacobNumDiff=numDiff(x,@(xState)aCVPolar(xState),4);
% err=(aJacob-AJacobNumDiff)./AJacobNumDiff;
% max(abs(err(:)))
%One will see that the maximum error is on the order of 1.6e-11, indicating
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
rDot=x(3,:);
thetaDot=x(4,:);

rDDot=r.*thetaDot.^2;
thetaDDot=-(2./r).*rDot.*thetaDot;

aVal=[rDot;thetaDot;rDDot;thetaDDot];

if(nargout>1)
    N=size(x,2);
    
    if(N>1)
        error('Derivatives are only available for numPoints=1.')
    end
    
    aJacob=[0                    0,  1,                 0;
            0,                   0,  0,                 1;
            thetaDot^2,          0,  0,                 2*r*thetaDot;
            2*rDot*thetaDot/r^2, 0,  -((2*thetaDot)/r), -((2*rDot)/r)];

    if(nargout>2)
        aHess=zeros(4,4,4);

        aHess(:,:,1)=[zeros(2,4);
                      0,                        0,  0,                2*thetaDot;
                      -((4*rDot*thetaDot)/r^3), 0,  (2*thetaDot)/r^2, (2*rDot)/r^2];
        aHess(:,:,3)=[zeros(3,4);
                      2*thetaDot/r^2,0,0,-(2/r)];
        aHess(:,:,4)=[zeros(2,4);
                      2*thetaDot, 0,  0,        2*r;
                      2*rDot/r^2,   0,  -((2)/r), -((2*rDot)/r)];

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
