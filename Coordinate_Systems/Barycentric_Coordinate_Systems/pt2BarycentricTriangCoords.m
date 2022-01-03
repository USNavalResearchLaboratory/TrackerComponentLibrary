function [phi,dphi]=pt2BarycentricTriangCoords(x,v1,v2,v3)
%%PT2BARYCENTRICTIANGCOORDS Convert a point into barycentric points with
%           respect to the vertices of a triangle in 2D. This solves for a
%           specific set of phi such that
%           x=phi(1)*v1+phi(2)*v2+phi(3)*v3
%           and sum(phi)=1. If x is inside the triangle, then phi>=0.
%
%INPUTS: x The 2X1 point whose barycentric coordiantes are desired.
% v1,v2,v3 The three 2X1 vertices of the triangle. It doens't matter if
%          they are in clockwise or counterclockwise order.
%
%OUTPUTS: phi The 3X1 set of barycentric coordinates of x.
%        dPhi A 3X2 set of partial derivatives of phi (rows) with respect
%             to the components of x (columns).
%
%The barycentric coordinate system on a triangle is described in Section 2
%of [1].
%
%EXAMPLE:
%In this example, coordinates inside a triangle are obtained. It is
%verified that the returned derivatives are (within finite precision
%limits) within the bounds of what one gets with finite diferencing. Also,
%the normalization of the weights is verified and the relative error of the
%interpolated values is computed and is also negligible within finite
%precision limits.
% v1=[-2;3];
% v2=[11;0];
% v3=[1;1];
% x=[4;1.5];
% [phi,dphi]=pt2BarycentricTriangCoords(x,v1,v2,v3);
% 
% f=@(y)pt2BarycentricTriangCoords(y,v1,v2,v3);
% dPhiNumDiff=numDiff(x,f,3);
% 
% RelDerivError=max(max(abs((dPhiNumDiff-dphi)./dPhiNumDiff)))
% normalizationError=sum(phi)-1
% 
% interpPoint=phi(1)*v1+phi(2)*v2+phi(3)*v3;
% RelErrVal=max(abs((interpPoint-x)./x))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

Av=triangleArea(v1,v2,v3);

phi=zeros(3,1);
phi(1)=triangleArea(x,v2,v3)/Av;
phi(2)=triangleArea(x,v3,v1)/Av;
phi(3)=triangleArea(x,v1,v2)/Av;

if(nargout>1) 
    dphi=zeros(3,2);
    
    dphi(1,:)=[(v2(2)-v3(2)),(v3(1)-v2(1))]/(2*Av);
    dphi(2,:)=[(v3(2)-v1(2)),(v1(1)-v3(1))]/(2*Av);
    dphi(3,:)=[(v1(2)-v2(2)),(v2(1)-v1(1))]/(2*Av);
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
