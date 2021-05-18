function [z,exitCode]=nearestPointOnEllipsoid(zIn,A,pIn,gammaVal,varargin)
%%NEARESTPOINTONELLIPSOID Given an ellipsoid in an arbitrary number of
%        dimensions such that a point zp on the surface of the ellipsoid
%        satisfies the equation (zIn-zp)'*A*(zIn-zp)=gammaVal find the
%        point on the ellipsoid that is closest to another point p. All
%        inputs are real.
%
%INPUTS: zIn The numDimX1 center of the ellipsoid.
%          A A numDimXnumDim symmetric, positive definite matrix that
%            specifies the size and shape of the ellipse or ellipsoid,
%            where a point zp is on the ellipse/ellipsoid if
%            (zIn-zp)'*A*(zIn-zp)=gammaVal.
%        pIn A numDimX1 point.
%   gammaVal The threshold for declaring a point to be on the ellipsoid. If
%            this parameter is omitted or an empty matrix is passed, the
%            default value of 1 is used.
%   varargin Any parameters that one might wish to pass to the fzero
%            function. These are generally comma-separated things, such as
%           'TolX',1e-12.
%
%OUTPUTS: z The numDimX1 point on the ellipse that is closest to p. When
%           multiple solutions exist, only one is chosen.
%  exitCode This function calls exitCode. This is the same as the exit flag
%           returned by the fzero function.
%
%The optimization problem being solved is
%minimize norm(z-Pin)^2
%such that (z-zIn)'*A*(z-zIn)=gammaVal
%To solve this problem, we will perofrm a few changes of variables. First,
%we perform an upper triangular Cholesky decomposition of A, so A=S'*S.
%Next, we make the substitutions
%zTilde=S*(z-zIn) (in reverse, z=inv(S)*zTilde+zIn).
%PInTilde=pIn-zIn
%This means that the minimization problem becomes
%norm(inv(S)*zTilde-PInTilde)^2
%such that norm(zTilde)^2=gammaVal.
%That formulation of the optimization problem can be solved using the
%eqConstLSSpher function. The solution to the original problem is then
%obtained by undoing the transformations.
%
%EXAMPLE 1:
%Here, we find the nearest point on an ellipsoid in 3D.
% A=[27,  4,  10;
%     4, 21, 16;
%    10, 16, 15];
% z=[12;24;36];
% p=[1000;-1000;2000];
% gammaVal=2.5;
% zp=nearestPointOnEllipsoid(z,A,p,gammaVal)
%One gets zp=[11.606225228425036;22.829459247176374;37.594179506992610]
%One can verify that (z-zp)'*A*(z-zp)=gammaVal within finite precision
%limits.
%
%EXAMPLE 2:
%This example is used in the function nearestPointInEllipsoid. There, since
%the point chosen in inside of the 2D ellipse, the function returns that
%point. However, this function returns the point on the outside of the
%ellipse that is closest.
% A=[1,0;
%    0,5];
% z=[0;10];
% p=[0.5;10];
% gammaVal=1;
% zp=nearestPointOnEllipsoid(z,A,p,gammaVal);
% figure(1)
% clf
% hold on
% axis([-6,6,-6+10,6+10])
% axis square
% drawEllipse(z,A,gammaVal,'g','linewidth',2)
% plot([p(1),zp(1)],[p(2),zp(2)],'--c')
% scatter(p(1),p(2),'ok','linewidth',2)
% scatter(zp(1),zp(2),'xr','linewidth',2)
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(gammaVal))
    gammaVal=1;
end

S=chol(A,'upper');
PInTilde=pIn-zIn;

[zTilde,exitCode]=eqConstLSSpher(inv(S),PInTilde,sqrt(gammaVal),varargin{:});
z=S\zTilde+zIn;

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
