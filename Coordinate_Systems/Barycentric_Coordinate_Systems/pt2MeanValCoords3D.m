function phi=pt2MeanValCoords3D(x,v,vertInFaces,epsValEq,epsValPi,epsValS)
%%PT2MANVALCOORDS3D Convert a point in or near a (not necessarily convex)
%       polyhedron to mean value coordinates in 3D  This is a type of
%       barycentric coordinate system that is described in [1] and [2].
%       These coordinates can be useful for Barycentric mapping and for
%       interpolation. This function only supports polyhedra with
%       triangular faces.
%
%INPUTS: x The 3XN set of points to convert to 3D Wachspress coordinates.
%        v The 3XnumVert set of vertices in the polyhedron.
% vertInFaces A 3XnumFaces array of indices on which vertices are in which
%          faces of the polyhedron. The ordering of the vertices within
%          each triangl should be counterclockwise as viewed from outside
%          the polyhedron.
% epsValEq An optional threshold parameter to avoid finite precision
%          issues. If dMax is the maximum distance between x and a vertex,
%          and di is the distance between x and the ith vertex, if di/dMax
%          is less than or equal to this value, then the point is just
%          taken to be the vertex. The default if omitted or an empty
%          matrix is passed is 2*eps().
% epsValPi An optional threshold parameter to avoid finite precision
%          issues. This is a threshold for how close a value is to pi to
%          determine whether or not a point lies on a face. The default if
%          omitted or an empty matrix is passed is 1024*eps(pi).
%  epsValS An optional threshold parameter to avoid finite precision
%          issues. This value is a threshold for skipping triangles that
%          would otherwise lead to infinite weights, It is a test related
%          to whether or not a point is coplanar with a particular face,
%          but does not lie in that face. The default if omitted or an mpty
%          matrix is passed is 2*eps().
%          
%OUTPUTS: phi A numVertXN set of points converted into 3D Wachspress
%             coordinates. These sould be valid inside and outside th
%             polyhedron.
%
%This function implements the algorithm of Fig. 4 of [1], except the
%section handling the 2D barycentic case has been replaced by
%Equation 2, because the formulae in Fig. 4 are incorrect.
%
%EXAMPLE:
%In this example, we convert points to 3D mean-value coordinates and then
%convert them back. A non-convex polyhedron is created, plotted and used.
%Sample points inside the polyhedron as well as on a vertex, outside the
%polyhedron and on a side of the polyhedron are used. The absolute error is
%computed and displayed and is on the order of what one would expect due to
%finite precision limitations. Then, a point just above a side of the
%polyhedron (but not on the side) is used and the absolute error is
%computed. It is far enough from the side that the algorithm doesn't just
%assume it is on the side. One can see a loss of precision (essentially
%from double precision to single precision) is present for points very near
%but not on the sides.
% %List of vertices.
% v=[-3,1,1,-3,1,1,-3,-3,1/2;
%     0,0,0, 0,1,1, 1, 1,1/2;
%     0,0,4, 1,0,1, 1, 0,2];
% x=[[-1;0.5;0.75],[0.5;0.5;1.5],[0.1;0.5;0],v(:,3),[-1;1;2]];
% %How the vertices combine to form triangular faces. All vertices are
% %ordered in counterclockwise order when viewed from outside the polyhedron.
% vertInFaces=[8,1,1,1,2,6,5,8,8,4,4,3,6,7;
%              5,5,2,3,5,3,8,7,1,7,3,6,7,4;
%              1,2,3,4,6,2,6,6,4,8,9,9,9,9];
% numFaces=size(vertInFaces,2);
% 
% figure(1)
% clf
% hold on
% %Draw the edges in the faces.
% for k=1:numFaces
%     idx=[vertInFaces(:,k);vertInFaces(1,k)];
%     plot3(v(1,idx),v(2,idx),v(3,idx),'-b','linewidth',2)
% end
% %Scatter the vertices.
% scatter3(v(1,:),v(2,:),v(3,:),400,'.k')
% %Scatter the points being tested
% scatter3(x(1,:),x(2,:),x(3,:),400,'.r')
% view(25,9)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% 
% phi=pt2MeanValCoords3D(x,v,vertInFaces);
% xBack=barycentricCoords2Pt(phi,v);
% AbsErrBack=max(max(abs((xBack-x))))
% 
% x=[0.1;0.5;2e-6];%Point near but not on a side.
% scatter3(x(1),x(2),x(3),200,'.g')
% phi=pt2MeanValCoords3D(x,v,vertInFaces);
% xBack=barycentricCoords2Pt(phi,v);
% AbsErrBack=max(max(abs((xBack-x))))
%
%REFERENCES:
%[1] T. Ju, S. Schaefer, and J. Warren, "Mean value coordinates for closed
%    triangular meshes," in Proceedings of the Special Interest Group
%    on Computer Graphics and Interactive Techniques Conference, Los
%    Angeles, CA, 31 Jul. - 4 Aug. 2005, pp. 561-566.
%[2] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(epsValS))
    epsValS=2*eps();
end

if(nargin<5||isempty(epsValPi))
    epsValPi=1024*eps(pi);
end

if(nargin<5||isempty(epsValEq))
    epsValEq=2*eps(max(abs(v(:))));
end

numPts=size(x,2);
numVertices=size(v,2);

phi=zeros(numVertices,numPts);
for k=1:numPts
    phi(:,k)=meanVal3DOnePt(x(:,k),v,vertInFaces,epsValEq,epsValPi,epsValS);
end
end

function phi=meanVal3DOnePt(x,v,vertInFaces,epsValEq,epsValPi,epsValS)

numVertices=size(v,2);
numFaces=size(vertInFaces,2);

u=bsxfun(@minus,v,x);
d=sqrt(sum(u.*u,1));

eqIdx=find(d/max(d)<=epsValEq,1);
%If the input is essentially a vertex.
if(~isempty(eqIdx))
    phi=zeros(numVertices,1);
    phi(eqIdx)=1;
    return;
end
u=bsxfun(@rdivide,u,d);

l=zeros(3,1);
c=zeros(3,1);
phi=zeros(numVertices,1);
for curFace=1:numFaces
    vertIdx=vertInFaces(:,curFace);
    i1=vertIdx(1);
    i2=vertIdx(2);
    i3=vertIdx(3);
    
    l(1)=norm(u(:,i2)-u(:,i3));
    l(2)=norm(u(:,i3)-u(:,i1));
    l(3)=norm(u(:,i1)-u(:,i2));
    theta=2*asin(l/2);
    h=sum(theta)/2;
    sinTheta=sin(theta);
    
    if(pi-h<=epsValPi)
        %Use 2D barycentric coordinates for this face. The solution in Fig.
        %4 of [1] for this step is incorrect. Instead, we use the solution
        %from Eq. 2.
        phi=zeros(numVertices,1);
        alpha(1)=angBetweenVecs(u(:,i1),u(:,i2));
        alpha(2)=angBetweenVecs(u(:,i2),u(:,i3));
        alpha(3)=angBetweenVecs(u(:,i3),u(:,i1));
        
        tanAlpha2=tan(alpha/2);
        phi(i1)=(tanAlpha2(3)+tanAlpha2(1))/d(i1);
        phi(i2)=(tanAlpha2(1)+tanAlpha2(2))/d(i2);
        phi(i3)=(tanAlpha2(2)+tanAlpha2(3))/d(i3);

        phi=phi/sum(phi);
        return;
    end
    c(1)=sin(h-theta(1))/prod(sinTheta([2,3]));
    c(2)=sin(h-theta(2))/prod(sinTheta([3,1]));
    c(3)=sin(h-theta(3))/prod(sinTheta([1,2]));
    c=c*2*sin(h)-1;
    s=sign(det([u(:,i1),u(:,i2),u(:,i3)])).*sqrt(1-c.^2);
    
    if(all(abs(s))<=epsValS)
         %Skip this triangle.
         continue;
     end
    
    phi(i1)=phi(i1)+(theta(1)-c(2)*theta(3)-c(3)*theta(2))/(d(i1)*sinTheta(2)*s(3));
    phi(i2)=phi(i2)+(theta(2)-c(3)*theta(1)-c(1)*theta(3))/(d(i2)*sinTheta(3)*s(1));
    phi(i3)=phi(i3)+(theta(3)-c(1)*theta(2)-c(2)*theta(1))/(d(i3)*sinTheta(1)*s(2));
end
phi=phi/sum(phi);

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
