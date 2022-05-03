function tetraIdxList=createTetrahedralMeshGrid(numX,numY,numZ,mixedMesh)
%%CREATETETRAHEDRALMESHGRID Create a mesh of tetrahedra over a regular
%       grid of cubes in 3D. The values returned are just integers starting
%       from 1 in (x,y,z) and can be used as indices in a more general
%       grid.
%
%INPUTS: numX, numY, numZ The number of points in the x, y and z dimensions
%              All of these must be >=2.
%      mixedMesh If this is false, then every cube is broken into 5
%                tetrahedra triangles in the same way. If this is true,
%                then neighboring cubes have opposite decompositions. The
%                default if omitted or an empty matrix is passed is true. 
%
%OUTPUTS: tetraIdxList A 3X4X5(numX-1)X(numY-1)X(numZ-1) matrix, where the
%            last three indices select which cube is being tetrahedralized.
%            The third index selects each of the 5 tetrahedra in each cube,
%            the first of which is the one in the center of the cubature
%            that doesn't have an outside face. tetraIdxList(:,i,k) is the
%            ith vertex of the kth tetrahedron.
%
%EXAMPLE:
%This creates and plots each type of mesh. The inner tetrahedron of each
%box is colored differently to make the difference between the patterns
%more visible.
% numX=5;
% numY=2;
% numZ=3;
% numCubes=(numX-1)*(numY-1)*(numZ-1);
% for curPlot=1:2
%     if(curPlot==1)
%         tetraIdxList=createTetrahedralMeshGrid(numX,numY,numZ,true);
%     else
%         tetraIdxList=createTetrahedralMeshGrid(numX,numY,numZ,false);
%     end
%     
%     figure(curPlot)
%     clf
%     hold on
%     tetraIdxSel=1;
%     for k=1:numCubes
%         for tetraIdx=1:5
%             if(tetraIdx==tetraIdxSel)
%                 style='-b';
%                 lineW=4;
%             else
%                 style='-g';
%                 lineW=1;
%             end
%             for i1=1:3
%                 for i2=(i1+1):4
%                     plot3([tetraIdxList(1,i1,tetraIdx,k),tetraIdxList(1,i2,tetraIdx,k)],[tetraIdxList(2,i1,tetraIdx,k),tetraIdxList(2,i2,tetraIdx,k)],[tetraIdxList(3,i1,tetraIdx,k),tetraIdxList(3,i2,tetraIdx,k)],style,'linewidth',lineW)
%                 end
%             end
%         end
%     end
%     axis equal
%     view(-34,35)
%     if(curPlot==1)
%         title('Mixed Mesh')
%     else
%         title('Non-Mixed Mesh')
%     end
% end
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(mixedMesh))
    mixedMesh=true; 
end

if(numX<2)
    error('numX must be >=2.')
end

if(numY<2)
    error('numY must be >=2.')
end

if(numZ<2)
    error('numZ must be >=2.')
end

%These are offsets for the vertices of a 3D box. 
a=[0;0;0];
b=[1;0;0];
c=[1;1;0];
d=[0;1;0];
e=[0;0;1];
f=[1;0;1];
g=[1;1;1];
h=[0;1;1];

%Arrangements of vertices for the two types of tetrahedra.
tetrahedra=zeros(3,4,5,2);
%Arrangement of vertices for the first mesh of 5 tetrahedra in the cube.
tetrahedra(:,:,1,1)=[a,f,c,h];
tetrahedra(:,:,2,1)=[a,b,c,f];
tetrahedra(:,:,3,1)=[a,h,d,c];
tetrahedra(:,:,4,1)=[a,f,h,e];
tetrahedra(:,:,5,1)=[f,g,c,h];

%Arrangement of vertices for the second mesh of 5 tetrahedra in the cube.
tetrahedra(:,:,1,2)=[b,g,d,e];
tetrahedra(:,:,2,2)=[b,e,d,a];
tetrahedra(:,:,3,2)=[b,f,g,e];
tetrahedra(:,:,4,2)=[e,g,d,h];
tetrahedra(:,:,5,2)=[b,g,c,d];

%Five tetrahedrons per 8-point cube.
tetraIdxList=zeros(3,4,5,numX-1,numY-1,numZ-1);

if(mixedMesh)
    for curZ=1:(numZ-1)
        for curY=1:(numY-1)
            for curX=1:(numX-1)
                selIdx=mod(curX+curY+curZ,2)+1;
                tetraIdxList(:,:,:,curX,curY,curZ)=bsxfun(@plus,[curX;curY;curZ],tetrahedra(:,:,:,selIdx));
            end
        end
    end
else
    for curZ=1:(numZ-1)
        for curY=1:(numY-1)
            for curX=1:(numX-1)
                tetraIdxList(:,:,:,curX,curY,curZ)=bsxfun(@plus,[curX;curY;curZ],tetrahedra(:,:,:,1));
            end
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
