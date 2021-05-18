function A=signedPolygonArea(vertices)
%%SIGNEDPOLYGONAREA Calculate the signed area of a non-self-intersecting
%                   2D polygon given its vertices. The magnitude of the
%                   area is the typical definition of area that one might
%                   expect. The sign of the area is positive if the points
%                   are arranged in a counterclockwise order around the
%                   polygon and negative if the points are in a clockwise
%                   order.
%
%INPUTS: vertices A 2XN list of the N vertices of the polygon in order
%                 around the polygon. Duplicate vertices are allowed. It
%                 does not matter whether or not the first vertex is
%                 repeated at the end.
%
%OUTPUTS: A The signed area of the polygon.
%
%The formula for computing the signed area of a non-self-intersecting
%polygon is taken from [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Polygon Area." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/PolygonArea.html
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVertices=size(vertices,2);

A=0;
prevVertex=vertices(:,end);
for curIdx=1:numVertices
    curVertex=vertices(:,curIdx);
    A=A+det([prevVertex,curVertex]);
    prevVertex=curVertex;
end
A=A/2;

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
