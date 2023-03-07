function J=uGradient2D(xTar,lRx,M,includeV)
%%UGRADIENT2D Determine the gradient of a direction cosine measurement in
%          2D with respect to 2D position. Relativity and atmospheric
%          effects are not taken into account.
%
%INPUTS: xTar The 2XN set of target position vectors in the global
%          coordinate system with [x;y] components for which a gradient is
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
%OUTPUTS: J A 1X2XN set of N Jacobian matrices (or a 2X2XN set if includeV
%           is true) where the row represents u (and v) and the columns
%           take the partial derivatives of the rows with respect to [x,y]
%           in that order.
%
%A derivation of the components of the Jacobian of a 3D u-v measurement is
%given in [1]. The results here are similar.
%
%EXAMPLE:
%Here, we validate the results using numerical differentiation.
% x=[1000;4];
% lRx=[-40;2];
% M=eye(2,2);
% J=uGradient2D(x,lRx,M);
% f=@(x)getUDirection2D(x,lRx,M);
% N=3;
% epsVal=[1;0.01];
% numJAug=numDiff(x,f,2,N,epsVal);
% JNumDiff=numJAug(2,:);
% RelErr=max(abs((J-JNumDiff)./JNumDiff))
%One will see that J and JNumDiff agree to more than 7 digits, indicating
%good agreeement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%April 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(xTar,2);

if(nargin<4||isempty(includeV))
    includeV=false;
end

if(nargin<3||isempty(M))
   M=eye(2,2); 
end

if(nargin<2||isempty(lRx))
   lRx=zeros(2,1); 
end

J=zeros(1+includeV,2,N);
for curPoint=1:N
    %Convert the state into the local coordinate system.
    xLocal=M*(xTar(1:2,curPoint)-lRx(1:2));

    r3=norm(xLocal)^3;
    x=xLocal(1);
    y=xLocal(2);

    %The gradient of u in local coordinates.
    du=zeros(1,2);
    du(1)=y^2/r3;
    du(2)=-x*y/r3;

    J(1,:,curPoint)=du*M;

    if(includeV)
        dv=zeros(1,2);
        dv(1)=-x*y/r3;
        dv(2)=x^2/r3;
        %Now, the gradient vectors for the angular components must be
        %rotated back into the global coordinate system.
        J(2,:,curPoint)=dv*M;
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
