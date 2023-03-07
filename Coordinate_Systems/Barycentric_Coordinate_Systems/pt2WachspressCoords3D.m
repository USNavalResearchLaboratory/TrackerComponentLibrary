function [phi,dphi]=pt2WachspressCoords3D(x,v,fAdj,numAdj,un,normalize)
%%PT2WACHSPRESSCOORDS3D Convert a point in a convex polyhedron to 3D
%       Wachspress coordinates. This is a type of barycentric coordinate
%       system that is described in [1] and in the 3D case in Section 8.1
%       in [2]. These coordinates can be useful for Barycentric mapping and
%       for interpolation.
%
%INPUTS: x The 3XN set of points to convert to 3D Wachspress coordinates.
%        v The 3XnumVert set of vertices in the polyhedron.
%     fAdj The maxAdjXnumVert set of indices of adjacent faces to each
%          vertex. The number of advancent faces per vertex is specified by
%          the the numAdj input. The faces for each vertex should be in
%          some type of counterclockwise order when looking at the vertex
%          from outside the polyhedron. For triangular faces, the
%          sortTriangFacesAroundVertices function can be used to get this
%          input and the next input in the correct order.
%   numAdj A numVertX1 list of the number of adjacent faces for each
%          vertex (the number of things in each row of fAdj).
%       un A 3XnumFaces set of unit normals to each face in the polyhedron.
%          The indices in fAdj select the normal for each face.
% normalize If true, an intermediate normalization step is performed. This
%          can reduce finite precision errors when handling badly scaled
%          data, but it also slows things down. The default if omitted or
%          an empty matrix is passed is false.
%
%OUTPUTS: phi A numVertXN set of points converted into 3D Wachspress
%             coordinates. These are typically only valid for points inside
%             the polyhedron.
%        dphi A numBasesX3XnumPts set of derivatives of the elements of phi
%             (rows) taken with respect to the elements of x (columns) for
%             each input measurement (3rd dimension).
%
%EXAMPLE:
%Create a convex polyhedron and consider two points inside it. The outline
%of the polyhedron is displayed along with normals to each of the faces
%(here all triangles), so that on can see that the ordering of the vertices
%in the triangles is correct (counterclockwise when looking at the
%polyhedron, so the technique for getting normals makes them all point
%out). Then, we convert the points into 3D Wachpress coordinates and
%consider the relative error of the reverse conversion. Also, the gradient
%at those points is evaluated and the absolute error compared to finite
%differencing is considered. The derivatives of phi with respect to x are
%considered as well as a second example where once phi have been fit they
%are used to interpolate with a different basis (here, randomly rotated
%basis vectors). The errors are on the order of what one might expect due
%to finite precision limitations.
%List of vertices.
% v=[-3,1,1,-3,1,1,-3,-3,1/2;
%     0,0,0, 0,1,1, 1, 1,1/2;
%     0,0,1, 1,0,1, 1, 0,2];
% numVert=size(v,2);
% x=[[-1;0.5;0.75],[0.5;0.5;1.5]];
% %How the vertices combine to form triangular faces. All vertices are
% %ordered in counterclockwise order when viewed from outside the polyhedron.
% vertInFaces=[8,1,1,1,2,6,5,8,8,4,4,3,6,7;
%              5,5,2,3,5,3,8,7,1,7,3,6,7,4;
%              1,2,3,4,6,2,6,6,4,8,9,9,9,9];
% numFaces=size(vertInFaces,2);
% %Get all of the normal vectors. They are all pointing out of the
% %polyhedron due to the ordering of the vertices.
% un=zeros(3,numFaces);
% for k=1:numFaces
%     i1=vertInFaces(1,k);
%     i2=vertInFaces(2,k);
%     i3=vertInFaces(3,k);
%     un(:,k)=-cross(v(:,i1)-v(:,i2),v(:,i3)-v(:,i2));
%     un(:,k)=un(:,k)/norm(un(:,k));%Normalize
% end      
% 
% %Get a table of which faces are around which vertex, sorted in
% %counterclockwise order around each vertex.
% [fAdj,numAdj]=sortTriangFacesAroundVertices(vertInFaces);
% 
% figure(1)
% clf
% hold on
% %Draw the edges in the faces.
% for k=1:numFaces
%     idx=[vertInFaces(:,k);vertInFaces(1,k)];
%     plot3(v(1,idx),v(2,idx),v(3,idx),'-b','linewidth',2)
% end
% %Draw the normal vectors leaving the centroids of all of the faces in
% %magenta. Also, draw dashed lines from each of the vertices to the centroid
% %points and mark the centroid points in red.
% for k=1:numFaces
%     centPt=mean(v(:,vertInFaces(:,k)),2);
%     scatter3(centPt(1),centPt(2),centPt(3),200,'.r');
%     for curVert=1:3
%         vC=v(:,vertInFaces(curVert,k));
%         plot3([vC(1),centPt(1)],[vC(2),centPt(2)],[vC(3),centPt(3)],'--c')
%     end
%     %Now, plot the normal to the face.
%     quiver3(centPt(1),centPt(2),centPt(3),un(1,k),un(2,k),un(3,k),0.5,'-m','linewidth',2)
% end
% %Scatter the vertices.
% scatter3(v(1,:),v(2,:),v(3,:),400,'.k')
% view(25,9)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% 
% [phi,dphi]=pt2WachspressCoords3D(x,v,fAdj,numAdj,un);
% xBack=barycentricCoords2Pt(phi,v);
% RelErrBack=max(max(abs((xBack-x)./x)))
% 
% f=@(x)pt2WachspressCoords3D(x,v,fAdj,numAdj,un);
% dPhiNumDiff=zeros(numVert,3,2);
% dPhiNumDiff(:,:,1)=numDiff(x(:,1),f,numVert);
% dPhiNumDiff(:,:,2)=numDiff(x(:,2),f,numVert);
% AbsErrDeriv=max(abs((dPhiNumDiff(:)-dphi(:))))
%
% R=randRotMat(3);
% vAlt=R*v;%A different basis set for which the phi are used to interpolate.
% h=@(x)barycentricCoords2Pt(pt2WachspressCoords3D(x,v,fAdj,numAdj,un),vAlt);
% hDeriv=zeros(3,3,2);
% hDeriv(:,:,1)=barycentricCoords2Pt(dphi(:,:,1),vAlt);
% hDeriv(:,:,2)=barycentricCoords2Pt(dphi(:,:,2),vAlt);
% dhNumDiff=zeros(3,3,2);
% dhNumDiff(:,:,1)=numDiff(x(:,1),h,3);
% dhNumDiff(:,:,2)=numDiff(x(:,2),h,3);
% AbsErrDeriv=max(abs((dhNumDiff(:)-hDeriv(:))))
%
%REFERENCES:
%[1] M. S. Floater, A. Gillette, and N. Sukumar, "Gradient bounds for
%    Wachspress coordinates on polytopes," SIAM Journal on Numerical
%    Analysis, vol. 52, no. 1, pp. 515-532, 2014.
%[2] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(normalize))
    normalize=false;
end

numPts=size(x,2);
numVert=size(v,2);
maxFaces=size(fAdj,1);

if(normalize)
    %Normalize to improve finite precision accuracy.
    minVals=min(v,[],2);
    maxVals=max(v,[],2);
    spanVals=maxVals-minVals;
    v=bsxfun(@rdivide,v-minVals,spanVals);
    x=bsxfun(@rdivide,x-minVals,spanVals);
    un=bsxfun(@times,un,spanVals);
    %Renormalize the scaled normal vectors.
    un=bsxfun(@rdivide,un,sqrt(sum(un.*un,1)));
end

%Temporary space used in the loops.
wiv=zeros(numVert-2,1);
Rl=zeros(3,maxFaces-2);
R=zeros(numVert,3);
p=zeros(3,maxFaces);

%For the return values.
phi=zeros(numVert,numPts);
dphi=zeros(numVert,3,numPts);
for curPt=1:numPts
    for curVert=1:numVert
        numFacesCur=numAdj(curVert);
        for curFace=1:numFacesCur
            h=dot(v(:,curVert)-x(:,curPt),un(:,fAdj(curFace,curVert)));
            p(:,curFace)=un(:,fAdj(curFace,curVert))/h;
        end

        for k=1:(numFacesCur-2)
            wiv(k)=det([p(:,k),p(:,k+1),p(:,numFacesCur)]); 
            Rl(:,k)=p(:,k)+p(:,k+1)+p(:,numFacesCur);
        end
        
        if(numFacesCur==3)
            %Triangular face. This formulation gets rid of a possible 0/0
            %issue if one of the weights is small.
            phi(curVert,curPt)=wiv(1);
            R(curVert,:)=Rl(:,1)';
        else
            phi(curVert,curPt)=sum(wiv(1:(numFacesCur-2)));
            R(curVert,:)=(Rl(:,1:(numFacesCur-2))*wiv(1:(numFacesCur-2)))'/phi(curVert,curPt);
        end
    end
    phi(:,curPt)=phi(:,curPt)/sum(phi(:,curPt));
    
    phiRSum=sum(bsxfun(@times,phi(:,curPt),R),1);
    dphi(:,:,curPt)=phi(:,curPt).*bsxfun(@minus,R,phiRSum);

    if(normalize)
        %Undo the normalization
        dphi(:,:,curPt)=bsxfun(@rdivide,dphi(:,:,curPt),spanVals');
    end
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
