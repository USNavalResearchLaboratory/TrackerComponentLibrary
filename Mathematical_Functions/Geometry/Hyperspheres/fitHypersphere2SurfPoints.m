function [r,center]=fitHypersphere2SurfPoints(points)
%%FITHYPERSPHERE2SURFPOINTS Given at least n+1 points on the surface of a
%  hypersphere in n-dimensional space, obtain the center and radius of the
%  sphere. The solution is exact for n+1 points and is a linear
%  approximation when given more than n+1 points. Note that if all 4 points
%  are coplanar, then an infinite number of solutions exist, though this
%  function only returns a single solution.
%
%INPUTS: points An nXnumPts set of n-dimensional points to which a sphere
%               should be fit. numPts>=n+1.
%
%OUTPUTS: r The scalar radius of the fitted sphere.
%    center The nX1 location of the center of the sphere.
%
%The equation for a sphere is r=norm(x-xc), where x is a point on the
%sphere, xc is the center, and r is the radius. Squaring both sides, one
%gets
%r^2=x'*x+xc*xc-2*x'*xc
%Let xi be the ith point on the sphere and x1 the first point. Subtracting
%the equation for the first point from the ith point, one gets
%0=xi'*xi-x1'*x1-2*(xi-x1)'*xc
%This is a linear equation in xc and with enough points forms a system that
%can be solved for xc with a linear least squares solution (use a
%pseudoinverse).
%
%To get r given xc, we minimize the sum of (r-norm(xi-xc)) for i=1 to
%numPts  with respect to r. That just ends up being
%r=(1/numPts)*sum_{i=1}^numPts norm(xi-xc)
%
%EXAMPLE:
%Here, we take 4 random points on a sphere and recover the sphere.
% center=[10;-20;12];
% r=6;
% numPts=4;
% points=zeros(3,numPts);
% for k=1:numPts
%     %Random points on the sphere.
%     points(:,k)=center+r*randRotMat(3)*[0;0;1];
% end
% %One can see that rBack=r and centerBack=center (within finite precision
% %limits).
% [rBack,centerBack]=fitHypersphere2SurfPoints(points)
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts=size(points,2);

A=bsxfun(@minus,points(:,2:numPts),points(:,1)).';

x11=points(:,1)'*points(:,1);
b=(sum(points(:,2:numPts).*points(:,2:numPts),1)-x11).'/2;

center=pinv(A)*b;
r=sum(sqrt(sum(bsxfun(@minus,points,center).^2,1)),2)/numPts;

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
