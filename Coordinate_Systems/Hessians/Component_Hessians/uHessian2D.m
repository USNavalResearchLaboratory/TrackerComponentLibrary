function HTotal=uHessian2D(xG,lRx,M)
%%UHESSIAN2D Determine the Hessian (a matrix of second derivatives) of a
%          direction cosine in 2D with respect to 2D position. Relativity
%          and atmospheric effects are not taken into account.
%
%INPUTS: xG The 2XN set of target position vectors in the global coordinate
%          system with [x;y] components for which a Hessian matrix is
%          desired.
%      lRx The 2X1 position vector of the receiver. If omitted, the
%          receiver is placed at the origin.
%        M A 2X2 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. If
%          omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: HTotal A 2X2XN set of Hessian matrices, one for each value in xG.
%                The partial derivatives for each i of HTotal(:,:,i) are
%                ordered:  [d2/(dxdx), d2/(dxdy);
%                           d2/(dydx), d2/(dyd)];
%                The matrix is symmetric, because d2/(dydx)=d2/(dxdy).
%
%A derivation of the components of the Jacobian of a 3D u-v measurement is
%given in [1]. The results here are similar and are simply one derivative
%higher.
%
%EXAMPLE:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% x=[100;-1000];
% lRx=[500;20;];
% epsVal=1e-6;
% M=randRotMat(2);
% 
% H=uHessian2D(x,lRx,M);
% J=uGradient2D(x,lRx,M);
% JdX=uGradient2D(x+[epsVal;0],lRx,M);
% JdY=uGradient2D(x+[0;epsVal],lRx,M);
% HNumDiff=zeros(2,2);
% HNumDiff(1,1)=(JdX(1)-J(1))/epsVal;
% HNumDiff(2,2)=(JdY(2)-J(2))/epsVal;
% HNumDiff(1,2)=(JdX(2)-J(2))/epsVal;
% HNumDiff(2,1)=HNumDiff(1,2);
% 
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will be on the order of 5e-7, indicating good
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

N=size(xG,2);

if(nargin<3||isempty(M))
   M=eye(2,2); 
end

if(nargin<2||isempty(lRx))
   lRx=zeros(2,1); 
end

HTotal=zeros(2,2,N);
for curPoint=1:N
    %Convert the state into the local coordinate system.
    xLocal=M*(xG(1:2,curPoint)-lRx(1:2));

    r5=norm(xLocal)^5;
    x=xLocal(1);
    y=xLocal(2);
    
    H=zeros(2,2);
    H(1,1)=-((3*x*y^2)/r5);
    H(2,2)=-((x*(x^2-2*y^2))/r5);
    H(1,2)=-((y*(-2*x^2+y^2))/r5);
    H(2,1)=H(1,2);
    
    %Convert back to global coordinates.
    HTotal(:,:,curPoint)=M'*H*M;
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
