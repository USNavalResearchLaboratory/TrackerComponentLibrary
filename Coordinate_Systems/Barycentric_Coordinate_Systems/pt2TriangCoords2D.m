function [phi,dphi]=pt2TriangCoords2D(x,v)
%%PT2TRIANGCOORDS2D Convert a 2D point to barycentric coordinates in a
%      triangle using the relative areas of subtriangles, as discussed at
%      the beginning of [1]. A signed triangular area is used so that
%      points outside of the triangle can also be represented (albeit with
%      some negative weights).
%
%INPUTS: x The 2XnumPts set of points to convert into triangular
%          coordinates.
%        v The 2X3 set of vertices of the triangle.
%
%OUTPUTS: phi A 3XnumPts set of the coordinates in the triangle. The
%             coordinates sum to 1, but the magnitudes of the individual
%             elements can be over 1.
%        dphi A 3X2XnumPts set of derivatives of the elements of phi
%             (rows) taken with respect to the elements of x (columns) for
%             each input measurement (3rd dimension).
%
%EXAMPLE 1:
%In this example, the coordinates of three points are found and the points
%are recreated from the phi values. The residual error of the reverse
%conversion in this example is zero.
% v=[0,1,3;
%    0,2,3];
% x=[[1;1],[1;1.5],[2;0.5]];
% phi=pt2TriangCoords2D(x,v);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-x)))
%
%EXAMPLE 2:
%This example demonstrates that the gradient returned by pt2TriangCoords2D
%is consistent with that computed via numerical differentiation. The
%relative error is on the order of what one might expect due to finite
%precision limitations.
% v=[0,1,3;
%    0,2,3];
% x=[[1;1],[1;1.5],[2;0.5]];
% numX=size(x,2);
% [~,dphi]=pt2TriangCoords2D(x,v);
% f=@(y)pt2TriangCoords2D(y,v);
% dPhiNumDiff=zeros(3,2,numX);
% for k=1:numX
%     dPhiNumDiff(:,:,k)=numDiff(x(:,k),f,3);
% end
% RelErr=max(max(max(abs((dphi-dPhiNumDiff)./dPhiNumDiff))))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

phiAll=triangleArea(v(:,1),v(:,2),v(:,3));

numX=size(x,2);
phi=zeros(3,numX);
dphi=zeros(3,2,numX);
for curX=1:numX
    %The phi values will be negative or positive depending on the ordering
    %of the vertices.
    [phi1,dphi1]=triangleArea(x(:,curX),v(:,2),v(:,3));
    [phi2,dphi2]=triangleArea(x(:,curX),v(:,3),v(:,1));
    [phi3,dphi3]=triangleArea(x(:,curX),v(:,1),v(:,2));

    phi(:,curX)=[phi1;phi2;phi3]/phiAll;
    dphi(:,:,curX)=[dphi1;
                    dphi2;
                    dphi3]/phiAll;
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
