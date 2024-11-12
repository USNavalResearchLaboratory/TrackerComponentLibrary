function extrema=findEllipsoidCartExtrema(R,zBar,gammaVal)
%%FINDELLIPSOIDCARTEXTREMA Given an ellipsoid defined by the region
%       specified by (z-zBar)'*inv(R)*(z-zBar)<=gamma, determine the points
%       on the ellipsoid that are at extrema in the Cartesian dimensions.
%       These are not the farthest points from the center of the ellipsoid,
%       just the farthest in each Cartesian direction.
%
%INPUTS: R A real NXN positive definite symmetric matrix.
%     zBar The nX1 location fo the center of the ellipsoid. if this
%          parameter is omitted or an empty matrix is passed, then the
%          default of zeros(N,1) is used.
% gammaVal The threshold determining the size of the ellipsoid. If this
%          parameter is omitted or an empty matrix is passed, the default
%          of 1 is used.
%
%OUTPUTS: An NX(2*N) set of all of the extreme points on the ellipsoid.
%
%To find the farthest point in the positive direction in the ith dimension,
%then, initially taking zar=0, we are solving
%arg max_z e_i'*z such that z'*inv(R)*z=gamma.
%where e_i is a vector of zeros with a one in the ith element.
%The solution is zMax=R(:,i)*sqrt(gammaVal/R(i,i)), which is used here for
%the solution. Given this optimal z, -z is also a solution and to all the
%solutions, one must add zBar to solve the original problem.
%
%EXAMPLE:
%Here, we draw an ellipse in 2D and then scatter plot the extrema to show
%that they agree with the locations of the 
% R=[1,1.3;
%    1.3,3];
% zBar=[1;2]*0;
% gammaVal=ChiSquareD.invCDF(0.999,2);
% figure(1)
% clf
% hold on
% drawEllipse(zBar,inv(R),gammaVal,'linewidth',2)
% vertices=findEllipsoidCartExtrema(R,zBar,gammaVal);
% scatter(vertices(1,:),vertices(2,:),600,'.')
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(gammaVal))
    gammaVal=1;
end

if(nargin<2||isempty(zBar))
    numDim=length(R);
    zBar=zeros(numDim,1);
end

Rii=diag(R);
extrema=bsxfun(@plus,zBar,bsxfun(@times,R,sqrt(gammaVal./Rii(:).')));
extrema=[extrema,-extrema];

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
