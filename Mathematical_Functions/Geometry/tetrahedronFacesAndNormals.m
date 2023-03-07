function [un,vertInFaces,fAdj,numAdj]=tetrahedronFacesAndNormals(v)
%%TETRAHEDRONFACESANDNORMALS Given the 4 vertices of a tetrahedron,
%     determine normal vectors to each face. Also, arrange the vertices
%     into the four faces in a manner such that the vertices are in
%     counterclockwise order when viewed from outside of the tetrahedron.
%     Additionally, sort faces about each vertex in a counterclockwise
%     order (as viewed into the vertex from outside the tetrahedron).
%
%INPUTS: v A 3X4 set of points in 3D.
%
%OUTPUTS: un A 3X4 set of unit normal vectors for each triangular face in
%            the tetrahedron. The unit normals are pointing out of the
%            tetrahedron's faces.
% vertInFaces A 3X4 collection of indices speciying the vertices in each
%           triangle in the tetrahedron. Also, the vertices are sorted in a
%           counterclockwise manner when looking at the face from outwise
%           the tetrahedron.
% fAdj, numAdj The matrix fAdj is a 3X4 matrix, where for the kth vertex,
%           fAdj(:,k) is the set of indices of the faces that are adjacent
%           to the kth vertex in v. It effectively selects a column of
%           vertInFaces. numAdj is a 3X1 vector of 3s, which is just the
%           number of faces adjacent to each vertex.
%
%Normal vectors are obtained by a cross product relation on triangles. A
%dot product relation is then used involving the normal, one point on the
%triangle and the other point in the tetrahedron to determine whether the
%normal is pointing into or out of the tetrahedron. The
%sortTriangFacesAroundVertices function sorts the faces.
%
%EXAMPLE:
%Generate a tetrehedron from 4 random points, compute the nornmals
%and then plot the tetrahedron and the normals in 3D. One might have to
%rotate the figure in order to fully see that the normals are pointing away
%from the faces and not into the tetrahedron.
% v=randn(3,4);
% [un,vertInFaces]=tetrahedronFacesAndNormals(v);
% 
% %Draw each face and then draw a normal arrow coming out of the face.
% figure(1)
% clf
% hold on
% for curFace=1:4
%     v1=v(:,vertInFaces(1,curFace));
%     v2=v(:,vertInFaces(2,curFace));
%     v3=v(:,vertInFaces(3,curFace));
% 
%     plot3([v1(1),v2(1),v3(1),v1(1)],[v1(2),v2(2),v3(2),v1(2)],[v1(3),v2(3),v3(3),v1(3)],'-k')
% 
%     %Get the center of the triangle.
%     vCent=(v1+v2+v3)/3;
% 
%     quiver3(vCent(1),vCent(2),vCent(3),un(1,curFace),un(2,curFace),un(3,curFace),'-b','filled','maxHeadSize',1)
% end
% axis equal
%
%April 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%All 4 faces contain 3 points (are triangles).
un=zeros(3,4);
vertInFaces=zeros(3,4);

%Normal to the 4-3-2 triangular face.
un(:,1)=cross(v(:,2)-v(:,4),v(:,3)-v(:,4));
%Determine whether the normal has the correct direction (it points away
%from the other point). If not, then flip the sign of the un(:,1) and swap
%the final two vertices of the face so that the triangle obeys a
%counterclockwise vertex ordering when looking at the vertices from outside
%the tetrahedron.
if(accurateDotProdK(un(:,1),v(:,1)-v(:,4))<0)
    %The normal vector is in the correct direction.
    vertInFaces(:,1)=[4;3;2];
else
    %The normal vector is in the wrong direction.
    un(:,1)=-un(:,1);
    vertInFaces(:,1)=[4;2;3];
end

%Normal to the 3-4-1 triangular face.
un(:,2)=cross(v(:,1)-v(:,3),v(:,4)-v(:,3));
if(accurateDotProdK(un(:,2),v(:,2)-v(:,3))<0)
    %The normal vector is in the correct direction.
    vertInFaces(:,2)=[3;4;1];
else
    %The normal vector is in the wrong direction.
    un(:,2)=-un(:,2);
    vertInFaces(:,2)=[3;1;4];
end

%Normal to the 2-1-4 triangular face.
un(:,3)=cross(v(:,4)-v(:,2),v(:,1)-v(:,2));
if(accurateDotProdK(un(:,3),v(:,3)-v(:,2))<0)
    %The normal vector is in the correct direction.
    vertInFaces(:,3)=[2;1;4];
else
    %The normal vector is in the wrong direction.
    un(:,3)=-un(:,3);
    vertInFaces(:,3)=[2;4;1];
end

%Normal to the 1-2-3 triangular face.
un(:,4)=cross(v(:,3)-v(:,1),v(:,2)-v(:,1));
if(accurateDotProdK(un(:,4),v(:,4)-v(:,1))<0)
    %The normal vector is in the correct direction.
    vertInFaces(:,4)=[1;2;3];
else
    %The normal vector is in the wrong direction.
    un(:,4)=-un(:,4);
    vertInFaces(:,4)=[1;3;2];
end

for k=1:4
    un(:,k)=un(:,k)/norm(un(:,k));%Normalize
end

%The vertices in the faces have been ordered in counterclockwise order when
%looking at the face from outside the tetrahedron. Now, record which
%faces are adjacent to which vertices, but store it in a sorted manner.
%First, we record it without sorting, then we sort it.
if(nargout>2)
    numAdj=[3;3;3;3];
    fAdj=[2,1,1,1;
          3,3,2,2;
          4,4,4,3];
    [fAdj,numAdj]=sortTriangFacesAroundVertices(vertInFaces,fAdj,numAdj);
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
