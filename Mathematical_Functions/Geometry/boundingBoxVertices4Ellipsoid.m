function vertices=boundingBoxVertices4Ellipsoid(R,zBar,gamma)
%%BOUNDINGBOXVERTICES4ELLIPSOID  Given an ellipsoid defined by the region
%       specified by (z-zBar)'*inv(R)*(z-zBar)<=gamma, find the vertices of
%       an axis-aligned bounding box that encloses the ellipsoid.
% 
%INPUTS: R A real NXN positive definite symmetric matrix.
%     zBar The nX1 location fo the center of the ellipsoid. if this
%          parameter is omitted or an empty matrix is passed, then the
%          default of zeros(N,1) is used.
% gammaVal The threshold determining the size of the ellipsoid. If this
%          parameter is omitted or an empty matrix is passed, the default
%          of 1 is used.
%
%OUTPUTS: vertices A 2X(2^N) matrix of the vertice sof a hyperrectangle
%                  that encloses the ellipsoid.
%
%To find the farthest point in the positive direction in the ith dimension,
%then, initially taking zBar=0, we are solving
%arg max_z e_i'*z such that z'*inv(R)*z=gamma.
%where e_i is a vector of zeros with a one in the ith element.
%The solution is zMax=R(:,i)*sqrt(gammaVal/R(i,i)), which is used here for
%the solution. The distance of this extremal point from the center of the
%ellipsoid is e_i*z=sqrt(gammaVal*R(i,i). Thus, knowing the distance of
%each dimension from the center, all vertices can be formed by putting all
%the distance sin a vector and calling PMCombos, which generates a vector
%of all signes. Then, one must add zBar back to shift the mean to the true
%mean location.
%
%EXAMPLE:
%Here, a 2D ellipse and its bounding box are plotted. Also plotted at the
%extrema of the ellipse in each dimensions. One can see that the bounding
%box intersects the extrema.
% R=[1,1.3;
%    1.3,3];
% zBar=[1;2]*0;
% gammaVal=ChiSquareD.invCDF(0.999,2);
% 
% figure(1)
% clf
% hold on
% drawEllipse(zBar,inv(R),gammaVal,'linewidth',2)
% vertices=boundingBoxVertices4Ellipsoid(R,zBar,gammaVal);
% vertices=vertices(:,[1,2,4,3]);
% vertices=[vertices,vertices(:,1)];
% plot(vertices(1,:),vertices(2,:),'-k','linewidth',2)
% extrema=findEllipsoidCartExtrema(R,zBar,gammaVal);
% scatter(extrema(1,:),extrema(2,:),600,'.')
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(gamma))
    gamma=1;
end

if(nargin<2||isempty(zBar))
    numDim=length(R);
    zBar=zeros(numDim,1);
end

Rii=diag(R);

%The maximum offset from the mean in each dimension.
deltaMax=sqrt(gamma*Rii);
vertices=bsxfun(@plus,zBar,PMCombos(deltaMax));

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
