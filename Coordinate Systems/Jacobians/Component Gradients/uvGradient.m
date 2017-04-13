function J=uvGradient(x,lRx,M)
%%UVGRADIENT Determine the gradient of a 3D direction cosine measurement
%          with respect to 3D position. Higher order gradient terms are
%          not provided and are zero. Relativity and atmospheric effects
%          are not taken into account.
%
%INPUTS: x The 3X1 target position vector in the global coordinate system
%          with [x;y;z] components.
%      lRx The 3X1  position vector of the receiver. If omitted, the
%          receiver is placed at the origin.
%        M A 3X3 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. If
%          omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: J A 2X3 Jacobian matrix where the rows are [u;v] in that order
%           and the columns take the derivative of the rows component with
%           respect to [x,y,z] in that order.
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

if(nargin<3||isempty(M))
   M=eye(3,3); 
end

if(nargin<2||isempty(lRx))
   lRx=zeros(3,1); 
end

%Convert the state into the local coordinate system.
xLocal=M*(x(1:3)-lRx(1:3));

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
    
J=[du;dv];

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
