function [zProj,AProj]=projEllipsOntoSpace(z,A,gammaVal,d,T)
%%PROJELLIPSONTOSPACE Determine the orthogonal projection of an ellipsoid
%                    onto an affine space. The projection itself will be an
%                    ellipsoid. For example, in 3D, the projection onto a
%                    plane will be an ellipse.
%
%INPUTS: z  The numDimX1 center of the ellipsoid.
%        A  A numDimXnumDim symmetric, positive definite matrix that
%           specifies the size and shape of the ellipse or ellipsoid, where
%           a point zp is on the ellipse/ellipsoid if
%           (zp-z)'*A*(zp-z)=gammaVal.
%  gammaVal The threshold for declaring a point to be in the ellipsoid. If
%           an empty matrix is passed, the default value of 1 is used.
%      d,T  A numDimX1 point d and a numDImXn matrix such that points in
%           the affine subspace are defined by x=d+T*t where t is an nX1
%           real vector, n<numDim.
%
%OUTPUTS: zProj, AProj Parameters defining the projected ellipsoid. zProj
%           is numDImX1 and AProj is numDimXnumDim. However, AProj will be
%           a singular matrix as everything resides on a subspace. A point
%           is on the projected ellipsoid if
%           (zp-zProj)'*AProj*(zp-zProj)=1
%           The example below shows gives an example of using such a
%           singular matrix to display the projection of a 3D ellipsoid
%           onto a 2D plane.
%
%The algorithm is taken from Section 13 of [1].
%
%EXAMPLE:
%We will project an ellipsoid in 3D onto a plane.
%Define the ellipsoid:
% A=[10,0,0;
%     0,20,0;
%     0,0,1];
% %A rotation matrix
% M=[-0.039622269933159, -0.378261925022885, -0.924850253718582;
%     0.765018603761816, -0.606901366847761,  0.215446668149313;
%    -0.642788154545232, -0.698991163746939,  0.313424219517309];
% A=M*A*M';
% c=[0;0;0];%The center of the ellipsoid
% %Define a plane.
% d=[1;0;3];
% T=[1,1;
%    0,1;
%    0,1];
% T=orth(T);%Must be othonormal.
% %Determine the projection of the ellipsoid onto the plane.
% [cProj,AProj]=projEllipsOntoSpace(c,A,1,d,T);
%
% figure(1)
% clf
% hold on
% axis([-4,4,-4,4,-4,4])
% %Draw the ellipsoid
% drawEllipse(c,A,1,'EdgeColor','none')
% view(104,-29)
% %Compute points on the plane to display.
% numPoints=50;
% tSpan=linspace(-10,10,numPoints);
% [tx,ty]=meshgrid(tSpan,tSpan);
% t=[tx(:)';ty(:)'];
% placeSurfVals=bsxfun(@plus,d,T*t);
% X=reshape(placeSurfVals(1,:),numPoints,numPoints);
% Y=reshape(placeSurfVals(2,:),numPoints,numPoints);
% Z=reshape(placeSurfVals(3,:),numPoints,numPoints);
% %Draw the plane
% surf(X,Y,Z)
% %Draw the projection of the ellipsoid onto the plane
% drawEllipse(cProj,AProj,1)
% %This makes the ellipsoid shiny
% light('Position',[0;1;-1])
%
%REFERENCES:
%[1] S. B. Pope, "Algorithms for ellipsoids," Cornell University, Tech.
%    Rep. FDA-08-01, Feb. 2008. [Online].
%    Available: https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(gammaVal))
   gammaVal=1; 
end

A=A/gammaVal;

L=chol(A,'lower');

%Make T into an orthonormal matrix. This does not rotate or shift the
%affine space at all. It just changes the parameterization as in x=d+T*t
%still holds, but if T is orthonormalized, then a different value of t
%is needed to get the same x.
T=orth(T);%Must be othonormal.

%Extracted from Equation 126
zProj=d+T*T'*(z-d);
[U,Sigma,~]=svd(T'/(L'),'econ');

L=T*U/Sigma;
AProj=L*L';

end