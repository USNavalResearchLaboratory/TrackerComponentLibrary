function AInterior=ellipsInsideConcEllips(A1,A2,gammaVal)
%%ELLIPSINSIDECONCELLIPS Given two concentric ellipsoids (they have the
%               same center), find the largest ellipsoid that fits inside
%               both of them. It will also have the same center.
%
%INPUTS: A1,A2 Two numDimXnumDim symmetric, positive definite matrices that
%           specify the size and shape of the ellipse or ellipsoid, where
%           a point zp is on the ellipse/ellipsoid if
%           (zp-z)'*A*(zp-z)=gammaVal.
%           %Both ellipsoids share the same center point, z.
%  gammaVal The threshold for declaring a point to be in the ellipsoid. If
%           this parameter is omitted or an empty matrix is passed, the
%           default value of 1 is used.
%
%OUTPUTS: AInterior A numDimXnumDIm symmetric, positive definite matrix
%            that specifies the size and shape of the largest ellipse that
%            is within both of the ellipses.  
%
%This function implements the algorithm that is given in Section 17 of [1].
%
%EXAMPLE:
% A1=[30,0;
%     0,1];
% A2=[1,0;
%     0,21];
% M=[0.413074198133900,  0.910697373904216;
%   -0.910697373904216,   0.413074198133900];%A rotation matrix
% A2=M*A2*M';
% AInterior=ellipsInsideConcEllips(A1,A2);
% figure(1)
% clf
% hold on
% drawEllipse([0;0],A1,1,'b','linewidth',4)
% drawEllipse([0;0],A2,1,'b','linewidth',4)
% drawEllipse([0;0],AInterior,1,'r','linewidth',2)
%
%REFERENCES:
%[1] S. B. Pope, "Algorithms for ellipsoids," Cornell University, Tech.
%    Rep. FDA-08-01, Feb. 2008. [Online].
%    Available: https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(gammaVal))
   gammaVal=1; 
end

A1=A1/gammaVal;
A2=A2/gammaVal;

L1=chol(A1,'lower');
L2=chol(A2,'lower');

L2Prime=L1\L2;

[U,Sigma,~]=svd(L2Prime);
SigmaTilde=diag(max(diag(Sigma).^2,1));
APrime=U*SigmaTilde*U';
LPrime=chol(APrime,'lower');

L=L1*LPrime;
AInterior=L*L';

end