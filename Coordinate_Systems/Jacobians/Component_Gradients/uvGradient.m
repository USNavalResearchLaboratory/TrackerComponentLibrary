function J=uvGradient(xG,lRx,M,includeW)
%%UVGRADIENT Determine the gradient of a direction cosine measurement with
%          respect to 3D position. Relativity and atmospheric effects
%          are not taken into account.
%
%INPUTS: xG The 3XN target position vectors in the global coordinate system
%           with [x;y;z] components for which gradients are desired.
%       lRx The 3X1 position vector of the receiver. If omitted, the
%           receiver is placed at the origin.
%         M A 3X3 rotation matrix from the global Coordinate system to the
%           orientation of the coordinate system at the receiver. If
%           omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: J A 2X3XN set of N Jacobian matrices where the rows are [u;v] in
%           that order and the columns take the partial derivative of the
%           rows with respect to [x,y,z] in that order. If includeW is
%           true, then this is a 3X3XN set of Jacobian matrices, where the
%           final row has partial derivatives of w.
%
%A derivation of the components of the Jacobian is given in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
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

if(includeW)
    J=zeros(3,3,N);
else
    J=zeros(2,3,N);
end
for curPoint=1:N
    %Convert the state into the local coordinate system.
    xLocal=M*(xG(1:3,curPoint)-lRx(1:3));

    r=norm(xLocal);
    x=xLocal(1);
    y=xLocal(2);
    z=xLocal(3);

    %The gradient of u in local coordinates.
    du=zeros(1,3);
    du(1)=(y^2+z^2)/r^3;
    du(2)=-x*y/r^3;
    du(3)=-x*z/r^3;

    %The gradient of v in local coordinates.
    dv=zeros(1,3);
    dv(1)=-x*y/r^3;
    dv(2)=(x^2+z^2)/r^3;
    dv(3)=-y*z/r^3;
  
    %Now, the gradient vectors for the angular components must be
    %rotated back into the global coordinate system.
    du=du*M;
    dv=dv*M;

    if(includeW)
        dw=zeros(1,3);
        dw(1)=-x*z/r^3;
        dw(2)=-y*z/r^3;
        dw(3)=(x^2+y^2)/r^3;
        dw=dw*M;
        J(:,:,curPoint)=[du;dv;dw];
    else
        J(:,:,curPoint)=[du;dv];
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
