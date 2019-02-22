function [xyPoints,adjMat]=getShaped2DLattice(dims,arrayType,dxdy,useTriGrid)
%%GETSHAPED2DLATTICE Obtain a lattice of points in a given shape. This
%               function is useful for obtaining element positions in
%               shapes that are often used for planar arrays. The points
%               can be in a rectangular or a triangular lattice.
%
%INPUTS: dims A 2X1 vector indicating the number of elements in the x and y
%             directions for array types other than hexagonal. For
%             hexagonal arrays, these are the number of elements at the top
%             and the number of elements in the middle.
%   arrayType A string indicating the type of shape of points to generate.
%             Shaped are filled unless indicated. Possible values are:
%             'rectangular' For a dims(1)Xdims(2) rectangular grid.
%             'elliptical' In this instance, dims/2 form the semi-major and
%                          semi-minor axes of an axis-aligned ellipse.
%              'hexagonal' Here dims are the number of elements at the top
%                          of the hexagon and the number of elements in the
%                          middle. The bottom has the same number of
%                          elements as the top.
%               'circular' (The default if omitted or an empty matrix is
%                          passed) The radius of the circle is
%                          min(dims.*dxdy).
%          'unif-circular' This returns an array that is circular and not
%                          filled. The elements are uniformly spaced.
%                          dims(1)>1 is the number of elements and dxdy(1)
%                          is the spacing between elements.
%       dxdy A 2X1 vector with the spacing between elements. The default if
%            omitted or an empty matrix is passed is [0.5;0.5], which is
%            appropriate for lambda/2 spacing in a narrowband antenna array
%            when distances are given in units of lambda.
% useTriGrid A boolean variable indicating whether a triangular lattice
%            should be used for the points. The default if omitted or an
%            empty matrix is passed is true. This is only used with filled
%            shapes.
%
%OUTPUTS: xyPoints A 2XnumPoints set of the points in the plane forming the
%                  selected shape with the given spacing and lattice.
%           adjMat An adjacency matrix. adjMat(i,j) is 1 if element i is
%                  adjacent to element j. Otherwise it is zero. The matrix
%                  is symmetric. adj(i,i) are always 1.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(arrayType))
    arrayType='circular';
end

if(nargin<3||isempty(dxdy))
   dxdy=[0.5;0.5];
end

if(nargin<4||isempty(useTriGrid))
    useTriGrid=true;
end

numPerRow=dims(1);

if(nargout>1)%if adjMat is desired
    dxdyOrig=dxdy;
    %Adjacency is simpler to determine if the spacing fo the points is
    %uniform. The final result can be scaled to the correct size.
    dxdy=[1;1];
end
dx=dxdy(1);

switch(arrayType)
    case 'rectangular'
        numPerCol=dims(2);
        dy=dxdy(2);
        if(useTriGrid)
            xLocsEven=0:dx:((numPerRow-1)*dx);
            stepSize=xLocsEven(2)-xLocsEven(1);
            xLocsOdd=xLocsEven+stepSize/2;
            yLocs=0:dy:((numPerCol-1)*dy);
            
            xyPoints=zeros(2,numPerRow*numPerCol);
            curStart=1;
            for curY=1:numPerCol
                y=yLocs(curY);

                if(mod(curY,2)==0)
                    xyPoints(:,curStart:(curStart+numPerRow-1))=[xLocsEven;y*ones(1,numPerRow)];
                else
                    xyPoints(:,curStart:(curStart+numPerRow-1))=[xLocsOdd;y*ones(1,numPerRow)];
                end
                curStart=curStart+numPerRow;
            end
        else
            xLocs=0:dx:((numPerRow-1)*dx);
            yLocs=0:dy:((numPerCol-1)*dy);
            [Y,X]=meshgrid(yLocs,xLocs);

            xyPoints=[X(:)';Y(:)'];
        end

        %Center about the origin.
        xyPoints=bsxfun(@minus,xyPoints,mean(xyPoints,2));
    case 'circular'
        numPerCol=dims(2);
        dy=dxdy(2);
        xyPoints=getShaped2DLattice(dims,'rectangular',dxdy,useTriGrid);
        
        R=min(numPerRow*dx/2,numPerCol*dy/2);
        sel=sum(xyPoints.^2,1)<=R^2;
        xyPoints=xyPoints(:,sel);
    case 'elliptical'
        numPerCol=dims(2);
        dy=dxdy(2);
        xyPoints=getShaped2DLattice(dims,'rectangular',dxdy,useTriGrid);
        
        Rx=numPerRow*dx/2;
        Ry=numPerCol*dy/2;
        
        sel=(xyPoints(1,:)/Rx).^2+(xyPoints(2,:)/Ry).^2<=1;
        xyPoints=xyPoints(:,sel);
    case 'unif-circular'
        numPoints=dims(1);
        theta=2*pi/numPoints;
        
        xyPoints=zeros(2,numPoints);
        
        %The radius.
        r=(dx/2)/sin(theta/2);
        
        %The first point is on the x axis.
        xyPoints(1,1)=r;
        M=rotMat2D(theta);
        
        for curPoint=2:numPoints
           xyPoints(:,curPoint)=M*xyPoints(:,curPoint-1);
        end
    case 'hexagonal'
        dy=dxdy(2);
        
        numBottom=min(dims);
        numMiddle=max(dims);
        
        if(useTriGrid)
            numElsInRows=[numBottom:numMiddle, (numMiddle-1):-1:numBottom];
            numEls=sum(numElsInRows);
            xyPoints=zeros(2,numEls);

            endVals=cumsum(numElsInRows);
            startVals=endVals-numElsInRows+1;

            rowIdx=1;
            for k=(numBottom-numMiddle):(numMiddle-numBottom)
                idx=startVals(rowIdx):endVals(rowIdx);
                xyPoints(1,idx)=((-(numElsInRows(rowIdx)-1)/2):((numElsInRows(rowIdx)-1)/2))*dx;
                xyPoints(2,idx)=k*dy;

                rowIdx=rowIdx+1;
            end
        else
            if((mod(numBottom,2)==0||mod(numMiddle,2)==0)&&~(mod(numBottom,2)==0&&mod(numMiddle,2)==0))
                numElsInRows=[numBottom:2:(numMiddle-1),(numMiddle-1):-2:numBottom];
            else
                numElsInRows=[numBottom:2:numMiddle, (numMiddle-2):-2:numBottom];
            end
            numEls=sum(numElsInRows);
            xyPoints=zeros(2,numEls);

            endVals=cumsum(numElsInRows);
            startVals=endVals-numElsInRows+1;
            
            rowIdx=1;
            for k=(numBottom-numMiddle):2:(numMiddle-numBottom)
                idx=startVals(rowIdx):endVals(rowIdx);
                
                xyPoints(1,idx)=((-(numElsInRows(rowIdx)-1)/2):((numElsInRows(rowIdx)-1)/2))*dx;
                xyPoints(2,idx)=(k/2)*dy;
                
                rowIdx=rowIdx+1;
            end 
        end
    otherwise
        error('Unknown array type specified.')
end

if(nargout>1)
    %If the adjacency matrix is desired.
    
    %Adjacency is determined via brute force. With a triagnular grid, in 
    %one row, elements are 1 wavelength apart. In the rows above and below,
    %the element is 1 wavelength up and 1/2 wavelength over. Thus, if we
    %just tested for things <=sqrt(5)/2 (=sqrt((1/2)^2+1^2)) then that
    %would be sufficient. With finite, precision, we test for a number that
    %is slightly higher, but is less than sqrt(2), which is the distance to
    %the closest non-adjacent element. A similar thing is done for a square
    %grid.
    
    if(useTriGrid)
        delta=1.26;
    else
        delta=1.7;
    end
    
    numEls=size(xyPoints,2);
    
    adjMat=zeros(numEls,numEls);
    for curEl1=1:(numEls-1)
        adjMat(curEl1,curEl1)=1;
        for curEl2=(curEl1+1):numEls
            dist=norm(xyPoints(:,curEl1)-xyPoints(:,curEl2));
            if(dist<=delta)
                adjMat(curEl1,curEl2)=1;
                adjMat(curEl2,curEl1)=1;
            end
        end
    end

    %Now, scale the points/
    xyPoints(1,:)=xyPoints(1,:)*dxdyOrig(1);
    xyPoints(2,:)=xyPoints(2,:)*dxdyOrig(2);
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
