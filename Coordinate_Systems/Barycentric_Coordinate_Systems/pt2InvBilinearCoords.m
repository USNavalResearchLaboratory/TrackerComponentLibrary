function [phi,dphi]=pt2InvBilinearCoords(x,v1,v2,v3,v4)
%%PT2INVBILINEARCOORDS Inverse bilinear coordinates are a type of 2D
%           barycentric coordinates on quadrilateral. This converts x,
%           which must be inside the specified quadrilateral, into inverse
%           bilinear coordinates. This solves for a specific set of phi
%           such that
%           x=phi(1)*v1+phi(2)*v2+phi(3)*v3+phi(4)*v4.
%           and such that sum(phi)=1 with all(phi>=0). This can be useful
%           for the interpolation of irregularly sampled multidimensional
%           functions of two parameters. One can interpolate to a point
%           inside of a specified region by weighting the function outputs
%           at the sample points with phi.
%
%INPUTS: x The 2X1 point whose inverse bilinear coordinates are desired.
% v1,v2,v3,v4 The four 2X1 vertices of the convex quadrilateral given in
%          COUNTERCLOCKWISE order. The ordering is important. Convexity can
%          be checked using the polygonIsConvex function and whether the
%          points are in counterclockwise order can be deduced from the
%          sign of signedPolygonArea. 
%
%OUTPUTS: phi The 4X1 set of bilinear coordinates of x.
%        dPhi A 4X2 set of partial derivatives of phi (rows) with respect
%             to the components of x (columns).
%
%Inverse bilinear coordinates are defined in Section 3 of [1].
%
%EXAMPLE:
%In this example, coordinates inside a simple quadrilateral are obtained.
%It is verified that the returned derivatives are (within finite precision
%limits) within the bounds of what one gets with finite diferencing. Also,
%the normalization of the weights is verified and the relative error of the
%interpolated values is computed and is also negligible within finite
%precision limits.
% v1=[0;0];
% v2=[11;0];
% v3=[1;1];
% v4=[0;1];
% 
% x=[0.2;0.75];
% [phi,dphi]=pt2InvBilinearCoords(x,v1,v2,v3,v4);
% f=@(y)pt2InvBilinearCoords(y,v1,v2,v3,v4);
% dPhiNumDiff=numDiff(x,f,4);
% 
% RelDerivError=max(max(abs((dPhiNumDiff-dphi)./dPhiNumDiff)))
% normalizationError=sum(phi)-1
% 
% interpPoint=phi(1)*v1+phi(2)*v2+phi(3)*v3+phi(4)*v4;
% RelErrVal=max(abs((interpPoint-x)./x))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A1=triangleArea(x,v1,v2);
A2=triangleArea(x,v2,v3);
A3=triangleArea(x,v3,v4);
A4=triangleArea(x,v4,v1);

B2=triangleArea(x,v1,v3);
B3=triangleArea(x,v2,v4);
B1=-B3;
B4=-B2;

sD=sqrt(abs(B1^2+B2^2+2*(A1*A3+A2*A4)));

E1=2*A1-B1-B2+sD;
E2=2*A2-B2-B3+sD;
E3=2*A3-B3-B4+sD;
E4=2*A4-B4-B1+sD;

phi=zeros(4,1);
phi(1)=4*A2*A3/(E2*E3);
phi(2)=4*A3*A4/(E3*E4);
phi(3)=4*A4*A1/(E4*E1);
phi(4)=4*A1*A2/(E1*E2);

if(nargout>1)
    %Derivatives with respect to [x(1),x(2)].
    dA1=[(v1(2)-v2(2))/2,(v2(1)-v1(1))/2];
    dA2=[(v2(2)-v3(2))/2,(v3(1)-v2(1))/2];
    dA3=[(v3(2)-v4(2))/2,(v4(1)-v3(1))/2];
    dA4=[(v4(2)-v1(2))/2,(v1(1)-v4(1))/2];
    
    dB2=[(v1(2)-v3(2))/2,(v3(1)-v1(1))/2];
    dB3=[(v2(2)-v4(2))/2,(v4(1)-v2(1))/2];
    dB1=-dB3;
    dB4=-dB2;
    
    dsD=(A3*dA1+A4*dA2+A1*dA3+A2*dA4+B1*dB1+B2*dB2)/sD;
    
    dE1=2*dA1-dB1-dB2+dsD;
    dE2=2*dA2-dB2-dB3+dsD;
    dE3=2*dA3-dB3-dB4+dsD;
    dE4=2*dA4-dB4-dB1+dsD;
    
    dphi=zeros(4,2);
    dphi(1,:)=phiDeriv(A2,dA2,A3,dA3,E2,dE2,E3,dE3);
    dphi(2,:)=phiDeriv(A3,dA3,A4,dA4,E3,dE3,E4,dE4);
    dphi(3,:)=phiDeriv(A4,dA4,A1,dA1,E4,dE4,E1,dE1);
    dphi(4,:)=phiDeriv(A1,dA1,A2,dA2,E1,dE1,E2,dE2);
end
end

function val=phiDeriv(Ai,dAi,Ai1,dAi1,Ei,dEi,Ei1,dEi1)
%%PHIDERIV The derivative of phi from the chain rule.
    val=(4.*Ai.*Ei.*Ei1.*dAi1-4.*Ai1.*(Ai.*Ei1.*dEi+Ei.*(Ai.*dEi1-Ei1.*dAi)))./(Ei.^2.*Ei1.^2);
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
