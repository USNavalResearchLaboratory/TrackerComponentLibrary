function u=getPolygonNormals(v)
%GETPOLYGONNORMALS Given the vertices in a polygon in 2D, obtain the
%     normals to the edges of the polygon. If the vertices are provided in
%     counterclockwise order, then the normals point outward. Otherwise,
%     the normals point inward. One can determine the order of the vertices
%     using the signedPolygonArea function.
%
%INPUTS: v A 2XnumVertices set of vertices of the polygon. numVertices>=3.
%
%OUTPUTS: u The normals to the edges of the polygon.
%
%EXAMPLE:
%With a set of vertices in counter-clockwise order, the normals are
%obtained and are seen all going out of the polygon. With the vertices in
%clockwise order, the normals are observed all going into the polygon. This
%plots a polygon and its normals with both vertex orderings.
% v=zeros(2,5);
% %Vertices in counter-clockwise order.
% v(:,1)=[-2;3];
% v(:,2)=[11;0];
% v(:,3)=[9;10];
% v(:,4)=[7;11];
% v(:,5)=[0;10];
% u=getPolygonNormals(v);
% figure(1)
% clf
% hold on
% plot([v(1,:),v(1,1)],[v(2,:),v(2,1)],'linewidth',2)
% midPoints=[(v(:,1:4)+v(:,2:5))/2,(v(:,5)+v(:,1))/2];
% quiver(midPoints(1,:),midPoints(2,:),u(1,:),u(2,:),0.5,'linewidth',2)
% axis equal
% %Vertices in clockwise order.
% v=v(:,5:-1:1);
% u=getPolygonNormals(v);
% figure(2)
% clf
% hold on
% plot([v(1,:),v(1,1)],[v(2,:),v(2,1)],'linewidth',2)
% midPoints=[(v(:,1:4)+v(:,2:5))/2,(v(:,5)+v(:,1))/2];
% quiver(midPoints(1,:),midPoints(2,:),u(1,:),u(2,:),0.5,'linewidth',2)
% axis equal
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVec=size(v,2);
u=zeros(2,numVec);
for i=1:numVec
    d=v(:,mod(i,numVec)+1)-v(:,i);
    u(:,i)=[d(2);-d(1)]/norm(d);
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
