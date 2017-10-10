function [aDeriv,aJacob]=aCVPolar(x,t)
%%ACVPOLAR The drift function for a continuous-time motion model where the
%          target state is given in polar coordinates and the motion is
%          constant velocity in 2D Cartesian coordinates. This function
%          works regardless of whether the angle is measured
%          counterclockwise from the x axis or clockwise form the y axis.
%
%INPUTS: x The 4XN state vector of N targets in the order of
%          [r;theta;rDot;thetaDot]. The angle theta is given in radians.
%        t An unused time component so that aCVPolar can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
%
%OUTPUTS: aDeriv The 4XN set of time-derivatives of the N state vectors
%                under the Cartesian linear motion model in polar
%                coordinates.
%         AJacob The 4X4XN set of Jacobians of aDeriv, which can be useful
%                in extended Kalman filters. This is the derivative of each
%                component of aDeriv (selected by row) with respect to the
%                elements of the state [r,theta,rDot,thetaDot] selected by
%                column.
%
%A derivation of the dynamic model is provided here. Let rVec be a position
%vector. We shall use the orthonormal basis vectors
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
% F=FPolyKal(T,xInitCart,1);
% systemType=0;
% xInitSpher=stateCart2Pol(xInitCart,systemType);
% xEndCart=F*xInitCart;
% RelTol=1e-10;
% AbsTol=1e-13;
% xStepsPolar=RKAdaptiveOverRange(xInitSpher,[0;T],@(x,t)aCVPolar(x,t),0.1,0,[],[],RelTol,AbsTol);
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
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

r=x(1,:);
rDot=x(3,:);
thetaDot=x(4,:);

rDDot=r.*thetaDot.^2;
thetaDDot=-(2./r).*rDot.*thetaDot;

aDeriv=[rDot;thetaDot;rDDot;thetaDDot];

if(nargout>1)
    N=size(x,2);
    aJacob=zeros(4,4,N);
    for k=1:N
        aJacob(:,:,k)=[0                                0,  1,                          0;
                       0,                               0,  0,                          1;
                       thetaDot(k)^2,                   0,  0,                          2*r(k)*thetaDot(k);
                       2*rDot(k)*thetaDot(k)/r(k)^2,    0,  -((2*thetaDot(k))/r(k)),    -((2*rDot(k))/r(k))];
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
