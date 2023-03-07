function boolVal=pointIsInTriangle(point2D,triVertices)
%%POINTISINTRIANGLE Determine whether a specified 2D points is within a
%           triangle.
%
%INPUTS: point2D A 2XnumPts set of 2D points.
%    triVertices A 2X3 set of vertices of the triangle in any order.
%
%OUTPUTS: bookVal A 1XnumPts set of boolean values that are true if the
%                 point is in the triangle or is on an edge of the
%                 triangle.
%
%The function pt2TriangCoords2D is used to convert the points to 2D
%triangular barycentric coordinates. If the point is within the triangle,
%then all of the cooridnates will be positive. Otherwise, some coordinates
%will be negative.
%
%EXAMPLE:
%In this example, the first 3 points are in or on the triangle and the
%final two points are outside the trangle. Thus, boolPoints=[1,1,1,0,0];
% v=[0,1,0.5;
%    0,0,0.5];
% points=[0.25,0,1,10,2;
%         0.25,0,0,-1,0];
% boolVals=pointIsInTriangle(points,v)
%
%September 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPts=size(point2D,2);
    phi=pt2TriangCoords2D(point2D,triVertices);
    
    boolVal=false(1,numPts);
    for k=1:numPts
        boolVal(k)=all(phi(:,k)>=0);
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
