function ANew=ellipsCoveringEllipsAndPoint(z,A,p,gammaVal,interiorOpt)
%%ELLIPSCOVERINGELLIPSEANDPOINT Given an ellipse (or ellipsoid) and a
%           point, find the minimum volume ellipse/ellipsoid that is
%           centered at the center of the original ellipse and that covers
%           the original ellipse and the point.
%
%INPUTS: z  The numDimX1 center of the ellipsoid.
%        A  A numDimXnumDim symmetric, positive definite matrix that
%           specifies the size and shape of the ellipse or ellipsoid, where
%           a point zp is on the ellipse/ellipsoid if
%           (zp-z)'*A*(zp-z)=gammaVal.
%        p  A numDimX1 point.
%  gammaVal The threshold for declaring a point to be in the ellipsoid. If
%           this parameter is omitted or an empty matrix is passed, the
%           default value of 1 is used.
% interiorOpt An optional parameter specifying what should be done if the
%           point p is inside of the ellipsoid. Possible values are
%           0 (The default if omitted or an empty matrix is passed) The
%             original ellipsoid is returned if the point is inside it.
%           1 The ellipse that is returned is the maximum volume ellipsoid
%             such that p is on the boundary and the ellipse is entirely
%             contained in the original ellipsoid.
%
%OUTPUTS: ANew The numDimXnumDim new positive definite matrix specifying
%              the shape of the minimum volume ellipsoid covering the
%              original ellipsoid and the point.
%
%The algorithm is taken from Section 10 of [1]. The basic idea is that
%first, we divide A by gammaVal so that the equation for the ellipse is
%(zp-z)'*A*(zp-z)=1. Next, let L be the lower-triangular Cholesky
%decomposition of A. We then do the transformation y=L'*(zp-z) and in the
%transformed domain, the ellipse is the unit sphere norm(y)=1. We can then
%apply the transformation to p to get pTilde=L'*(p-z). We basically now
%just want to transform p so that it is on the x-axis, and then expand the
%sphere in that direction minimally to intersect p and then transform back.
%
%The sphere in y has an equvalent matrix ATilde=eye(numDiom,numDim). We
%want to expand it just enough in a single direction to touch pTilde. This
%is a rank-one modification problem, which is described in a little detail
%in Sections 4 and 10 of [1].
%
%EXAMPLE 1:
%Here we have an ellipse and a point outside of the ellipse.
% A=[1,0;
%    0,10];
% M=[0.413074198133900,  0.910697373904216;
%    -0.910697373904216,   0.413074198133900];%A rotation matrix
% A=M*A*M';
% z=[5;6];
% p=[6;7];
% ANew=ellipsCoveringEllipsAndPoint(z,A,p);
% figure(1)
% clf
% hold on 
% drawEllipse(z,A,1,'b','linewidth',2)
% drawEllipse(z,ANew,1,'m','linewidth',2)
% scatter(p(1),p(2),'xr','linewidth',2)
%The new ellipse engulf the original ellipse and the point.
%
%EXAMPLE 2:
%In this instance, we have an ellipse and a point inside the ellipse.
% A=[1,0;
%    0,10];
% M=[0.413074198133900,  0.910697373904216;
%    -0.910697373904216,   0.413074198133900];%A rotation matrix
% A=M*A*M';
% z=[5;6];
% p=[5.2;6];
% ANew=ellipsCoveringEllipsAndPoint(z,A,p);
% figure(1)
% clf
% hold on 
% drawEllipse(z,A,1,'b','linewidth',2)
% drawEllipse(z,ANew,1,'m','linewidth',2)
% scatter(p(1),p(2),'xr','linewidth',2)
%The new ellipse is the same as the orignal ellipse, because the point is
%inside the original ellipse.
%
%REFERENCES:
%[1] S. B. Pope, "Algorithms for ellipsoids," Cornell University, Tech.
%    Rep. FDA-08-01, Feb. 2008. [Online].
%    Available: https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(interiorOpt))
    interiorOpt=0;
end

if(nargin<4||isempty(gammaVal))
   gammaVal=1; 
end

n=length(z);

A=A/gammaVal;

%Check for the point being inside the ellipse.
diff=p-z;
if(interiorOpt==0&&diff'*A*diff<=1)
    %The point is already inside the ellipse.
    ANew=A;
    return;
end

L=chol(A,'lower');

pTilde=L'*(p-z);%Eq 62.

pN=norm(pTilde);
g=(1/pN-1)*(1/pN^2);
G=eye(n,n)+g*(pTilde*pTilde');
ANew=L*G*(L*G)'*gammaVal;

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
