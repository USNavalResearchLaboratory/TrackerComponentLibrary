function J=uGradient3D(xG,lRx,M)
%%UGRADIENT3D Determine the gradient of just the u component of a direction
%          cosine measurement in 3D with respect to 3D position. Relativity
%          and atmospheric effects are not taken into account.
%
%INPUTS: xG The 3XN target position vectors in the global coordinate system
%           with [x;y;z] components for which gradients are desired.
%       lRx The 3X1 position vector of the receiver. If omitted, the
%           receiver is placed at the origin.
%         M A 3X3 rotation matrix from the global Coordinate system to the
%           orientation of the coordinate system at the receiver. If
%           omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: J A 1X3XN set of N Jacobian matrices where the rows is u and the
%           columns take the partial derivative of u with respect to
%           [x,y,z] in that order.
%
%A derivation of the components of the Jacobian is given in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2021 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(M))
    M=eye(3,3); 
end

if(nargin<2||isempty(lRx))
    lRx=zeros(3,1); 
end

N=size(xG,2);
J=zeros(1,3,N);

du=zeros(1,3);
for curPoint=1:N
    %Convert the state into the local coordinate system.
    xLocal=M*(xG(1:3,curPoint)-lRx(1:3));

    r=norm(xLocal);
    x=xLocal(1);
    y=xLocal(2);
    z=xLocal(3);

    %The gradient of u in local coordinates.
    du(1)=(y^2+z^2)/r^3;
    du(2)=-x*y/r^3;
    du(3)=-x*z/r^3;
    
    %Rotate the gradient back into the global coordinate system.
    du=du*M;

    J(:,:,curPoint)=du;
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
