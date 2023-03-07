function [phi,dphi]=pt2TetrahedronCoords3D(x,v)
%%PT2TETRAHEDRONCOORDS3D Convert a 3D point to barycentric coordinates in
%      a tetrahedron using the relative areas of subtetrahedra, as
%      discussed at the beginning of [1]. A signed tetrahedral volume area
%      is used so that points outside of the tetrahedron can also be
%      represented (albeit with some negative weights).
%
%INPUTS: x The 3XnumPts set of points to convert into tetrahedral
%          coordinates.
%        v The 3X4 set of vertices of the tetrahedron.
%
%OUTPUTS: phi A 4XnumPts set of the coordinates in the triangle.
%        dphi A 4X3XnumPts set of derivatives of the elements of phi
%             (rows) taken with respect to the elements of x (columns) for
%             each input measurement (3rd dimension).
%
%In [1], the coordinates are defined as the ratio of the volume of a number
%of sub-tetrahedra (with x as a vertex) and the tetrahedron given in v.
%Here, things differ a little bit in that signed tetrahedral volumes are
%used with a very specific ordering. This allows this function to work
%with points inside the tetrahedron (with all positive weights) as well as
%with points outside the tetrahdron (with some negative weights).
%
%EXAMPLE 1:
%In this example, the coordinates of three points are found and the points
%are recreated from the phi values. The residual error of the reverse
%conversion is seen to be on the order of finite precision errors.
% x=sqrt(3)/3;
% h=sqrt(6)/3;
% d=sqrt(3)/6;
% v=[x,  -d,  -d, 0;
%    0, 1/2,-1/2, 0;
%    0,   0,   0, h];
% xPts=[[0;0;0.5],[x;0;0],[1;1;1]];
% phi=pt2TetrahedronCoords3D(xPts,v);
% xBack=barycentricCoords2Pt(phi,v);
% ResidualErr=max(max(abs(xBack-xPts)))
%
%EXAMPLE 2:
%In this example, we verify that the derivative is consistent with
%numerical differentiation. The absolute error is consistent with what one
%might expect due to finite precision limitations.
% v=randn(3,4);
% x=randn(3,1);
% [~,dphi]=pt2TetrahedronCoords3D(x,v);
% f=@(y)pt2TetrahedronCoords3D(y,v);
% dPhiNumDiff=numDiff(x,f,4);
% RelErr=max(max(abs((dphi-dPhiNumDiff)./dPhiNumDiff)))
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

a=v(:,1);
b=v(:,2);
c=v(:,3);
d=v(:,4);

numX=size(x,2);
phi=zeros(4,numX);
dphi=zeros(4,3,numX);
for curX=1:numX
    [Va,dVa]=tetrahedronVolume([b,x(:,curX),d,c],false,2);
    [Vb,dVb]=tetrahedronVolume([a,x(:,curX),c,d],false,2);
    [Vc,dVc]=tetrahedronVolume([a,x(:,curX),d,b],false,2);
    [Vd,dVd]=tetrahedronVolume([a,x(:,curX),b,c],false,2);

    %Nominally, V=tetrahedronVolume(v,false); should equal Va+Vb+Vc+Vd.
    %However, to ensure that phi is normalized to more significant digits,
    %we compute V directly from the sum rather than computing it once
    %outside this loop and reusing it.
    V=Va+Vb+Vc+Vd;

    phi(:,curX)=[Va;Vb;Vc;Vd]/V;

    dphi(:,:,curX)=[dVa;
                    dVb;
                    dVc;
                    dVd]/V;
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
