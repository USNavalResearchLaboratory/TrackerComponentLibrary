function triIdxList=createTriangIdxMeshGrid(numRows,numCols,mixedMesh)
%%CREATETRIANGIDXMESHGRID Create a mesh of triangles over a regular
%       square grid. The values returned are just integers starting from 1
%       in rows and columns and can be used as indices in a grid of
%       deformed points that should be triangulated (assuming that the
%       deformation isn't such that triangles would intersect).
%
%INPUTS: numRows The number of rows of vertices in the grid. numRows>=2.
%        numCols The number of columns of vertices in the grid. numCols>=2.
%      mixedMesh If this is false, then every square is broken into two
%                triangles in the same way. If this is true, then
%                neighboring squares have opposite decompositions. The
%                default if omitted or an empty matrix is passed is true. 
%
%OUTPUTS: triIdxList A 2X3X2X(numRows-1)X(numCols-1) matrix of the
%                    triangles. The final two indices select which square
%                    is selected (in rows and columns). Each square
%                    contains two triangles (the third index). For each
%                    triangle triIdxList(:,vertexNum,triIdx) lists the two
%                    coordinates of each vertex. Vertices in vertexNum are
%                    listed in counterclockwise order.
%
%EXAMPLE:
%This creates and plots each type of mesh.
% numRows=5;
% numCols=11;
% numTriang=2*(numRows-1)*(numCols-1);
% 
% triIdxList=createTriangIdxMeshGrid(numRows,numCols,true);
% figure(1)
% clf
% hold on
% for curTri=1:numTriang
%     plot([triIdxList(1,:,curTri),triIdxList(1,1,curTri)],[triIdxList(2,:,curTri),triIdxList(2,1,curTri)],'-b')
% end
% axis equal
% title('Mixed Mesh')
% 
% triIdxList=createTriangIdxMeshGrid(numRows,numCols,false);
% figure(2)
% clf
% hold on
% for curTri=1:numTriang
%     plot([triIdxList(1,:,curTri),triIdxList(1,1,curTri)],[triIdxList(2,:,curTri),triIdxList(2,1,curTri)],'-b')
% end
% axis equal
% title('Non-Mixed Mesh')
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(mixedMesh))
    mixedMesh=true;
end

if(numRows<2)
    error('There must be at least 2 rows.')
end

if(numCols<2)
    error('There must be at least 2 columns.')
end

triIdxList=zeros(2,3,2,numRows-1,numCols-1);
if(mixedMesh)
    for curCol=1:(numCols-1)
        curCol1=curCol+1;
        for curRow=1:(numRows-1)
            curRow1=curRow+1;

            %The triangles in neighboring squares are mirrored so as to form a
            %better triangulation.
            if(mod(curRow+curCol,2))
                triIdxList(:,:,1,curRow,curCol)=[curRow,curRow1,curRow1;
                                                 curCol,curCol1,curCol];
                triIdxList(:,:,2,curRow,curCol)=[curRow,curRow,curRow1;
                                                 curCol,curCol1,curCol1];
            else
                triIdxList(:,:,1,curRow,curCol)=[curRow,curRow,curRow1;
                                                 curCol,curCol1,curCol];
                triIdxList(:,:,2,curRow,curCol)=[curRow,curRow1,curRow1;
                                                 curCol1,curCol1,curCol];
            end
        end
    end
else
    for curRow=1:(numRows-1)
        curRow1=curRow+1;
        for curCol=1:(numCols-1)
            curCol1=curCol+1;
            triIdxList(:,:,1,curRow,curCol)=[curRow,curRow1,curRow1;
                                             curCol,curCol1,curCol];
            triIdxList(:,:,2,curRow,curCol)=[curRow,curRow,curRow1;
                                             curCol,curCol1,curCol1];
        end
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
