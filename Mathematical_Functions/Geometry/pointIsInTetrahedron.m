function boolVal=pointIsInTetrahedron(p,V)
%%POINTISINTETRAHEDRON Given a tetrahedron defined by its vertices,
%       determine whether a specified 3D point is in the tetrahedron. A
%       point that is on the surface of the tetrahedron is considered to be
%       in the tetrahedron.
%
%INPUTS: p A 3XnumPts set of points to test for being inside the
%          tetrahedron.
%        V A 3X4 matrix of the vertices of the tetrahedron.
%
%OUTPUTS: boolVal A numPtsX1 vector of boolean values, where boolVal(i) is
%                 true if the ith point in p is in the tetrahedron and it
%                 is false otherwise.
%
%The algorithm is given in Chapter 13.4.1 of [1].
%
%EXAMPLE:
%Given a tetrahedron, we generate points on a grid inside and outside the
%tetrahedron. We then pass the points to pointIsInTetrahedron. Those points
%that are in the tetrahedron are then displayed in red.
% V=[1,-3, -1,  2;
%    2, 1, -1,  2;
%    3, 2,  0,  1];
% figure(1)
% clf
% hold on
% scatter3(V(1,:),V(2,:),V(3,:))
% for k1=1:3
%     for k2=(k1+1):4 
%         plot3([V(1,k1),V(1,k2)],[V(2,k1),V(2,k2)],[V(3,k1),V(3,k2)],'-k')
%     end
% end
% numPts=50;
% pts=linspace(-4,4,numPts);
% [X,Y,Z]=meshgrid(pts,pts,pts);
% xyzPts=[X(:).';Y(:).';Z(:).'];
% sel=pointIsInTetrahedron(xyzPts,V);
% scatter3(xyzPts(1,sel),xyzPts(2,sel),xyzPts(3,sel),'.r')
% view(28,17)
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=bsxfun(@minus,V(:,2:4),V(:,1));
if(det(M)<0)
    %The ordering must produce a positive determinant. Swapping two columns
    %will change the sign.
    temp=V(:,4);
    V(:,4)=V(:,3);
    V(:,3)=temp;
end

n=zeros(3,4);
n(:,1)=cross(V(:,2)-V(:,4),V(:,3)-V(:,4));
n(:,2)=cross(V(:,1)-V(:,3),V(:,4)-V(:,3));
n(:,3)=cross(V(:,4)-V(:,2),V(:,1)-V(:,2));
n(:,4)=cross(V(:,3)-V(:,1),V(:,2)-V(:,1));

numPts=size(p,2);
boolVal=true(numPts,1);
for curPt=1:numPts
    val=(dot(n(:,1),p(:,curPt)-V(:,4))>0)||(dot(n(:,2),p(:,curPt)-V(:,3))>0)||(dot(n(:,3),p(:,curPt)-V(:,2))>0)||(dot(n(:,4),p(:,curPt)-V(:,1))>0);
    boolVal(curPt)=~val;
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
