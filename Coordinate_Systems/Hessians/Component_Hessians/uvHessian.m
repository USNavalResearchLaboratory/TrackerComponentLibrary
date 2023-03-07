function HTotal=uvHessian(xG,lRx,M,includeW)
%%UVHESSIAN Determine the Hessian matrix (second derivative matrix) of a
%           direction cosine measurement with respect to 3D position.
%           Relativity and atmospheric effects are not taken into account.
%
%INPUTS: xG The 3XN target position vectors in the global coordinate system
%           with [x;y;z] components for which gradients are desired.
%       lRx The 3X1  position vector of the receiver. If omitted, the
%           receiver is placed at the origin.
%         M A 3X3 rotation matrix from the global coordinate system to the
%           orientation of the coordinate system at the receiver. If
%           omitted, it is assumed to be the identity matrix.
%  includeW An optional boolean value indicating whether a third direction
%           cosine component should be included. The u and v direction
%           cosines are two parts of a 3D unit vector. Generally, one might
%           assume that the target is in front of the sensor, so the third
%           component would be positive and is not needed. However, the
%           third component can be included if ambiguity exists. The
%           default if this parameter is omitted or an empty matrix is
%           passed is false.
%
%OUTPUTS: HTotal A 3X3X(2+includeW)XN matrix such that for the ith measurement,
%          HTotal(:,:,1,i) is the Hessian matrix with respect to the u
%          component and HTotal(:,:,2,i) is the Hessian matrix with respect
%          to the v component and, if included, HTotal(:,:,3,i) is the
%          Hessian matrix with respect to the w component. The elements in
%          the matrices for each component/ point are
%          ordered [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
%                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
%                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
%          note that each matrix is symmetric (i.e.
%           d^2/(dydx)=d^2/(dxdy) ).
%
%A derivation of the components of the Jacobian is given in [1]. The
%Hessian is just one derivative higher.
%
%EXAMPLE:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% xG=[100;-1000;500];
% lRx=[500;20;-400];
% epsVal=1e-6;
% M=randRotMat(3);
% includeW=true;
% 
% H=uvHessian(xG,lRx,M,includeW);
% J=uvGradient(xG,lRx,M,includeW);
% JdX=uvGradient(xG+[epsVal;0;0],lRx,M,includeW);
% JdY=uvGradient(xG+[0;epsVal;0],lRx,M,includeW);
% JdZ=uvGradient(xG+[0;0;epsVal],lRx,M,includeW);
% HNumDiff=zeros(3,3,2);
% for k=1:(2+includeW)
%     HNumDiff(1,1,k)=(JdX(k,1)-J(k,1))/epsVal;
%     HNumDiff(1,2,k)=(JdX(k,2)-J(k,2))/epsVal;
%     HNumDiff(2,1,k)=HNumDiff(1,2,k);
%     HNumDiff(1,3,k)=(JdX(k,3)-J(k,3))/epsVal;
%     HNumDiff(3,1,k)=HNumDiff(1,3,k);
%     HNumDiff(2,2,k)=(JdY(k,2)-J(k,2))/epsVal;
%     HNumDiff(2,3,k)=(JdY(k,3)-J(k,3))/epsVal;
%     HNumDiff(3,2,k)=HNumDiff(2,3,k);
%     HNumDiff(3,3,k)=(JdZ(k,3)-J(k,3))/epsVal;
% end
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will be on the order of 1e-6 or better, indicating good
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

if(nargin<4||isempty(includeW))
    includeW=false;
end

if(nargin<3||isempty(M))
    M=eye(3,3); 
end

if(nargin<2||isempty(lRx))
    lRx=zeros(3,1); 
end

N=size(xG,2);
HTotal=zeros(3,3,2+includeW,N);
for curPoint=1:N
    %Convert the state into the local coordinate system.
    xLocal=M*(xG(1:3)-lRx(1:3));
    
    x=xLocal(1);
    y=xLocal(2);
    z=xLocal(3);

    x2=x*x;
    y2=y*y;
    z2=z*z;
    
    r5=norm(xLocal)^5;

    H=zeros(3,3,2);

    %u Hessian values
    %du/(dxdx)
    H(1,1,1)=-((3*x*(y2+z2))/r5);
    %du/(dydy)
    H(2,2,1)=-((x*(x2-2*y2+z2))/r5);
    %du/(dzdz)
    H(3,3,1)=-((x*(x2+y2-2*z2))/r5);
    %du/(dxdy)
    H(1,2,1)=-((y*(-2*x2+y2+z2))/r5);
    %du/(dydx)
    H(2,1,1)=H(1,2,1);
    %du/(dxdz)
    H(1,3,1)=-((z*(-2*x2+y2+z2))/r5);
    %du/(dzdx)
    H(3,1,1)=H(1,3,1);
    %du/(dydz)
    H(2,3,1)=(3*x*y*z)/r5;
    %du/(dzdy)
    H(3,2,1)=H(2,3,1);

    %v Hessian values
    %dv/(dxdx)
    H(1,1,2)=-((y*(-2*x2+y2+z2))/r5);
    %dv/(dydy)
    H(2,2,2)=-((3*y*(x2+z2))/r5);
    %dv/(dzdz)
    H(3,3,2)=-((y*(x2+y2-2*z2))/r5);
    %dv/(dxdy)
    H(1,2,2)=-((x*(x2-2*y2+z2))/r5);
    %dv/(dydx)
    H(2,1,2)=H(1,2,2);
    %dv/(dxdz)
    H(1,3,2)=(3*x*y*z)/r5;
    %dv/(dzdx)
    H(3,1,2)=H(1,3,2);
    %dv/(dydz)
    H(2,3,2)=-((z*(x2-2*y2+z2))/r5);
    %dv/(dzdy)
    H(3,2,2)=H(2,3,2);

    %Rotate the values back into global coordinates
    HTotal(:,:,1,curPoint)=M'*H(:,:,1)*M;
    HTotal(:,:,2,curPoint)=M'*H(:,:,2)*M;

    %w Hessian values
    if(includeW)
        %dw/(dxdx)
        H(1,1,3)=-z*(-2*x2+y2+z2)/r5;
        %dw/(dydy)
        H(2,2,3)=-z*(x2-2*y2+z2)/r5;
        %dw/(dzdz)
        H(3,3,3)=-3*(x2+y2)*z/r5;
        %dw/(dxdy)
        H(1,2,3)=3*x*y*z/r5;
        %dw/(dydx)
        H(2,1,3)=H(1,2,3);
        %dw/(dxdz)
        H(1,3,3)=-x*(x2+y2-2*z2)/r5;
        %dw/(dzdx)
        H(3,1,3)=H(1,3,3);
        %dw/(dydz)
        H(2,3,3)=-y*(x2+y2-2*z2)/r5;
        %dw/(dzdy)
        H(3,2,3)=H(2,3,3);

        HTotal(:,:,3,curPoint)=M'*H(:,:,3)*M;
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
