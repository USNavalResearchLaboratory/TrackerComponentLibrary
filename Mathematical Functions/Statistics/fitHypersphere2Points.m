function [r,xCen,didConverge]=fitHypersphere2Points(x,relTol,absTol,maxIter)
%%FITHYPERSPHERE2POINTS Given a set of n-dimensional points, fit an
%           n-dimensional hypersphere to the points. In 2D, this is a
%           circle and in 3D, a sphere, but n can be >3. the fit is in an
%           unweighted least-squared sense. The algorithm assumes that the
%           points have sufficient diversity for the estimation problem. In
%           2D, this would mean that they are not all collinear, in 3D not
%           all coplanar, etc.
%
%INPUTS: x A nXnumPoints matrix of points in nD.
%   relTol Convergence is determined when the norm of the difference in r
%          and xCen across iterations is less than
%          relTol*norm([r;xCen]). The default if omitted or an empty matrix
%          is passed is 1e-10.
%   absTol Convergence is determined when the norm of the difference in r
%          and xCen across iterations is less than absTol. The default if
%          omitted or an empty matrix is passed is 1e-13.
%  maxIter The maximum number of iterations for the algorithm to perform.
%          if this parameter is omitted or an empty matrix is passed, then
%          maxIter=5000 is used.
%
%OUTPUTS: r The scalar radius of the estimated hypersphere.
%      xCen The nX1 center coordinate of the fitted hypersphere.
% didConverge This is true (1) if convergence occurred. Otherwise, this is
%           false (0). 
%
%This function minimizes the average squared difference between r and the
%range obtained from all of the points and the specified center of the
%circle. The algorithm is described in 2D in Appendix A.7.5  of [1] and in
%3D in Appendix A.7.6. the extension to an arbitrary nujmber of dimensions
%is straightforward.
%
%EXAMPLE 1:
%Here, we randomly sample a circle centered at (5,4) with radius 3. We add
%noise to the points but still get a decent fit. If noise were not added,
%then the result would be very close.
% numPoints=200;
% theta=2*pi*rand(1,numPoints);
% rTrue=3;
% xCenTrue=[5;4];
% x=0.1*randn(2,numPoints)+bsxfun(@plus,xCenTrue,rTrue*[cos(theta);sin(theta)]);
% [r,xCen,didConverge]=fitHypersphere2Points(x)
%
%EXAMPLE 2:
%Here, we fit to 4D data without noise added. the result is very close.
% xCenTrue=[5;4;-18;2];
% numDim=length(xCenTrue);
% numPoints=200;
% u=randDirVec(numDim,numPoints);
% rTrue=3;
% x=bsxfun(@plus,xCenTrue,rTrue*u);
% [r,xCen,didConverge]=fitHypersphere2Points(x)
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(maxIter))
    maxIter=5000;
end

if(nargin<3||isempty(absTol))
    absTol=1e-13;
end

if(nargin<2||isempty(relTol))
    relTol=1e-10;
end

xBar=mean(x,2);
xCen=xBar;

rPrev=Inf;
numDim=size(x,1);
xCenPrev=zeros(numDim,1);

didConverge=false;
for curIter=1:maxIter
    L=sqrt(sum(bsxfun(@minus,x,xCen).^2,1));
    r=mean(L);

    LBar=mean(bsxfun(@rdivide,xCen-x,L),2);
    
    xCen=xBar+r*LBar;
    
    penalty=norm([rPrev-r;xCenPrev-xCen]);
    if(penalty<=relTol*norm([r;xCen])||penalty<=absTol)
        didConverge=true;
        break; 
    end
    rPrev=r;
    xCenPrev=xCen;
end
end
