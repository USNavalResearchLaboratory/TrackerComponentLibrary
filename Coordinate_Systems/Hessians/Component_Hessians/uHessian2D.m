function H=uHessian2D(xG,lRx,M,includeV)
%%UHESSIAN2D Determine the Hessian (a matrix of second derivatives) of a
%          direction cosine in 2D with respect to 2D position. Relativity
%          and atmospheric effects are not taken into account.
%
%INPUTS: xG The 2X1 target position vector in the global coordinate
%          system with [x;y] components for which a Hessian matrix is
%          desired.
%      lRx The 2X1 position vector of the receiver. If omitted, the
%          receiver is placed at the origin.
%        M A 2X2 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. If
%          omitted, it is assumed to be the identity matrix.
% includeV An optional boolean value indicating whether a second direction
%          cosine component should be included. The u direction cosine is
%          one parts of a 2D unit vector. The default if this parameter is
%          omitted or an empty matrix is passed is false.
%
%OUTPUTS: H A 2X2X1 (or 2X2X2 is includeV is true) set of Hessian matrices.
%           The partial derivatives for each i of H(:,:,i) are ordered: 
%           [d2/(dxdx), d2/(dxdy);
%           d2/(dydx), d2/(dyd)];
%           where H(:,:,1) hold partial derivatives of u and H(:,:,2)
%           partial derivatives of v. Each matrix is symmetric, because
%           d2/(dydx)=d2/(dxdy).
%
%A derivation of the components of the Jacobian of a 2D u measurement is
%given in [1] (remove a component form the 3D measurement). The results
%here are similar and are simply one derivative higher.
%
%EXAMPLE:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% x=[100;-1000];
% lRx=[500;20;];
% epsVal=1e-5;
% M=randRotMat(2);
% includeV=true;
% 
% H=uHessian2D(x,lRx,M,includeV);
% J=uGradient2D(x,lRx,M,includeV);
% JdX=uGradient2D(x+[epsVal;0],lRx,M,includeV);
% JdY=uGradient2D(x+[0;epsVal],lRx,M,includeV);
% HNumDiff=zeros(2,2,1+includeV);
% for k=1:(1+includeV)
%     HNumDiff(1,1,k)=(JdX(k,1)-J(k,1))/epsVal;
%     HNumDiff(2,2,k)=(JdY(k,2)-J(k,2))/epsVal;
%     HNumDiff(1,2,k)=(JdX(k,2)-J(k,2))/epsVal;
%     HNumDiff(2,1,k)=HNumDiff(1,2,k);
% end
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will tend to be on the order of 1e-8, indicating good
%agreement between the numerical Hessian matrix and the actual Hessian
%matrix.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(includeV))
    includeV=false;
end

if(nargin<3||isempty(M))
    M=eye(2,2); 
end

if(nargin<2||isempty(lRx))
    lRx=zeros(2,1); 
end

%Convert the state into the local coordinate system.
xLocal=M*(xG(1:2)-lRx(1:2));

r5=norm(xLocal)^5;
x=xLocal(1);
y=xLocal(2);
    
H=zeros(2,2,1+includeV);

%Fill in the Hessian for the u component.
H(1,1,1)=-((3*x*y^2)/r5);
H(1,2,1)=-((y*(-2*x^2+y^2))/r5);
H(2,1,1)=H(1,2,1);
H(2,2,1)=-((x*(x^2-2*y^2))/r5);
%Convert back to global coordinates.
H(:,:,1)=M'*H(:,:,1)*M;

if(includeV)
    H(1,1,2)=-y*(-2*x^2+y^2)/r5;
    H(1,2,2)=-x*(x^2-2*y^2)/r5;
    H(2,1,2)=H(1,2,2);
    H(2,2,2)=-3*x^2*y/r5;
    %Convert back to global coordinates.
    H(:,:,2)=M'*H(:,:,2)*M;
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
