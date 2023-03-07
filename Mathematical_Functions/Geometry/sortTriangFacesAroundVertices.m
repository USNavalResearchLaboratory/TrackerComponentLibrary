function [fAdj,numAdj]=sortTriangFacesAroundVertices(vertInFaces,fAdj,numAdj)
%%SORTTRIANGFACESAROUNDVERTICES Consider a closed polyhedron that is made
%       up of triangles whose vertices are sorted in counterclockwise order
%       when looking from the outside the polyhedron at the triangle, given
%       a list of which vertices are in which triangle, obtain a list of
%       which faces (triangles) are adjacent to each vertex with the faces
%       specified in counterclockwise order.
%
%INPUTS: vertInFaces A 3XnumFaces set of indices of vertices in all of the
%                    triangles in the polyhedron (all faces are triangles).
%                    The vertex locations do not need to be provided; just
%                    the indices of the triangles.
%       fAdj, numAdj If available, fAdj is a maxAdjFacesXnumVert matrix,
%                    where for the kth vertex, fAdj(1:numAdj(k),k) is the
%                    set of indices of the faces that are adjacent to that
%                    vertex. It effectively selects a column of
%                    vertInFaces. numAdj is a length numVertices vector. If
%                    these two inputs are omitted, they are reconstructed
%                    from vertInFaces.
%
%OUTPUTS: fAdj,numAdj These have the same definition as fAdj and numAdj on
%                     the input but the faces are sorted in
%                     counterclockwise order around each vertex.
%
%Having the ordering of faces around vertices recorded in a
%counterclockwise order is important in algorithms such as in [1] for
%finding 3D barycentric coordinates.
%
%In a completely closed polyhedron, one can start at the vertex of interest
%and then follow the triangles around in the correct order. Assume that all
%of the trianges are given with indices in counterclockwise order. One can
%rotates the indices of each triangles adjacent to the vertex such that the
%vertex of interest is in the first index. One can then trace the triangles
%around in the correct order by following the indices. The second vertex of
%the first triangle will also be the last index traced. To trace around,
%consider the example below:
%     4-------3
%    / \     / \
%   /   \   /   \
%  /     \ /     \
% 5-------1-------2
%  \      |      /
%   \     |     /
%    \    |    /
%     \   |   /
%      \  |  /
%       \ | /
%        \|/
%         6
%
%Starting with face 1-2-3, the second triangle must share vertex 3. Thus,
%one must search for a triangle whose second index is 3. That leads to
%1-3-4. Similarly, the third index of that triangle will be 4 and the
%second index of the next triangle must be 4, so tht leads to triangle
%1-4-5. Thus, one can continue tracing triangles aorund the vertex until
%the common vertex is number 2, in which case one knows that the final
%triangle has been found.
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2022 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVert=max(vertInFaces(:));

if(nargin<2||isempty(fAdj))
    %Determine the maximum number of adjacent vertices.
    numAdj=zeros(numVert,1);
    for k=1:numVert
        idx=find(any(vertInFaces==k,1));
        numAdj(k)=length(idx);
    end

    %Store the adjacent vertices, but not in any particular order.
    maxAdj=max(numAdj);
    fAdj=zeros(maxAdj,numVert);
    for k=1:numVert
        idx=find(any(vertInFaces==k,1));
        fAdj(1:numAdj(k),k)=idx;
    end
end

for k=1:numVert
    numFaces=numAdj(k);
    %Rotate the vertices in each triangle until vertex k is first.
    for curFace=1:numFaces
        fIdx=fAdj(curFace,k);
        
        if(vertInFaces(2,fIdx)==k)
            vertInFaces(:,fIdx)=vertInFaces([2;3;1],fIdx);
        elseif(vertInFaces(3,fIdx)==k)
            vertInFaces(:,fIdx)=vertInFaces([3;1;2],fIdx);
        end%Otherwise, assume vertInFaces(1,fIdx)==k and no rotation is
           %needed.
    end
    newOrder=zeros(numFaces,1);
    newOrder(1)=1;%The first face remains unchanged.
    assignedFaces=false(numFaces,1);
    assignedFaces(1)=true;
    
    curAddedIdx=1;
    fIdx=fAdj(1,k);
    endNextIdx=vertInFaces(2,fIdx);
    curNextIdx=vertInFaces(3,fIdx);
    while(curNextIdx~=endNextIdx)
        for curFace=2:numFaces
            if(assignedFaces(curFace))
                continue;
            end
            fIdx=fAdj(curFace,k);
            
            if(vertInFaces(2,fIdx)==curNextIdx)
                %Found the next face.
                break;
            end
        end
        
        curAddedIdx=curAddedIdx+1;
        newOrder(curAddedIdx)=curFace;
        assignedFaces(curFace)=true;
        curNextIdx=vertInFaces(3,fIdx);
    end
    fAdj(1:numFaces,k)=fAdj(newOrder,k);
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
