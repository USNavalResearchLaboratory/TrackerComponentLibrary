function clippedPolygon=clipPolygonSH2D(polygon2Clip,convexClipPolygon,makeFirstEqualLast)
%%CLIPPOLYGONSH2D Clip a non-self-intersecting polygon to a region
%                 specified by a convex polygon using a 2D implementation
%                 of the Sutherland-Hodgman algorithm. If the polygon being
%                 clipped was concave at the vertices outside of the
%                 clipping region, the resulting clipped polygon will have
%                 overlapping edges on the boundary of the region rather
%                 than being split into multiple separate polygons.
%
%INPUTS: polygon2Clip A 2XN list of N>=3 vertices making up the polygon of
%                     the form [x;y]. It does not matter whether vertices
%                     are repeated.
%   convexClipPolygon A 2XNClip list of NClip>=3 vertices making up the
%                     convex clipping polygon. Vertices should not be
%                     repeated. The first vertex must not be repeated on
%                     the end. The vertices must be in a counterclockwise
%                     order.
%  makeFirstEqualLast If this is true, then the last element in the clipped
%                     polygon will be guaranteed to equal the first
%                     element, thus fully closing the polygon. The default
%                     if omitted or an empty matrix is passed is false.
%
%OUTPUTS: clippedPolygon The polygon polygon2Clip clipped to the convex
%                        region specified by convexClipPolygon. If
%                        polygon2Clip is completely does not intersect with
%                        the clipping region, an empty matrix is returned.
%                        If the clipping region is engulfed by the polygon,
%                        then the polygon will be the clipping region.
%
%The Sutherland-Hodgman algorithm is originally from [1], though the
%authors focus a lot of its attention on clipping polyhedra to planes. A
%quick search online will yield many result explaining how it works in 2D
%in a more intuitive manner. The original paper also includes a technique
%for dividing a polygon into multiple parts when it is split by the
%clipping region rather than just having coincident edges on the edge of
%the clipping region. 
%
%The Sutherland-Hodgman algorithm requires that the vertices in the
%clipping polygon be in counterclockwise order.
%
%The maximum number of vertices on the clipped polygon is if every other
%vertex of the polygon to clip is outside of the clipping region, as each
%exit/ entry to the clipping region adds two vertices. Thus, the maximum
%number of edges is 3/2 times the original number.
%
%REFERENCES:
%[1] I. E. Sutherland and G. W. Hodgman, "Reentrant polygon clipping,"
%    Communications of the ACM, vol. 17, no. 1, pp. 32-42, Jan. 1974.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(makeFirstEqualLast))
        makeFirstEqualLast=true;
    end

    %The clipping polygon must have its vertices going in a
    %counterclockwise order. Check whether that is the case.
    if(signedPolygonArea(convexClipPolygon)<0)
        error('The vertices of the clipping polygon should be in a counterclockwise order. Reverse the order of the vertices and try again.');
    end

    curPolygon=polygon2Clip;
    numClipVertices=size(convexClipPolygon,2);
    
    %The first clipping edge will be the one from the end to the beginning.
    prevClipVertex=convexClipPolygon(:,end);
    %For each edge, create the reduced polygon by clipping with that edge.  
    for curClip=1:numClipVertices
        curClipVertex=convexClipPolygon(:,curClip);
        %The polygon will be clipped to this edge.
        curClipEdge=curClipVertex-prevClipVertex;
        
        numVertices=size(curPolygon,2);
        
        %The maximum possible number of vertices in the clipped polygon.
        maxClippedVertices=ceil((3/2)*numVertices)+(makeFirstEqualLast>0);
        clippedPolygonNew=zeros(2,maxClippedVertices);
        
        %Number of vertices that have been added to clippedPolygonNew
        numVerticesNew=0;
        prevVertex=curPolygon(:,end);        
        for curV=1:numVertices
            curVertex=curPolygon(:,curV);
            
            %Note that the clip edges are infinitely extended so that
            %vertices outside of the clipping region are involved.
            if(vertexIsInsideClipEdge(curVertex,prevClipVertex,curClipEdge))
                %If the current vertex is inside of the clipping region,
                %then add it, but if the previous vertex was not in the
                %clipping region, then an extra vertex at the edge of the
                %boundary region needs to be added.
                if(~vertexIsInsideClipEdge(prevVertex,prevClipVertex,curClipEdge))
                    numVerticesNew=numVerticesNew+1;
                    clippedPolygonNew(:,numVerticesNew)=twoLineIntersectionPoint2D([prevVertex,curVertex],[curClipVertex,prevClipVertex]);

                end
                
                numVerticesNew=numVerticesNew+1;
                clippedPolygonNew(:,numVerticesNew)=curVertex;
            elseif(vertexIsInsideClipEdge(prevVertex,prevClipVertex,curClipEdge))
                %If the previous vertex was inside of the clipping region
                %and this vertex is not, then add a line segment from the
                %previous vertex to the edge of the clipping region.
                numVerticesNew=numVerticesNew+1;
                clippedPolygonNew(:,numVerticesNew)=twoLineIntersectionPoint2D([prevVertex,curVertex],[curClipVertex,prevClipVertex]);
            end

            prevVertex=curVertex;
        end
        
        %Resize to fit.
        clippedPolygonNew=clippedPolygonNew(:,1:numVerticesNew);
        curPolygon=clippedPolygonNew;
        
        %The object is not in the viewing area at all.
        if(isempty(curPolygon))
            clippedPolygon=[];
            return
        end
        
        prevClipVertex=curClipVertex;
    end
    
    if(makeFirstEqualLast)
        if(any(curPolygon(:,1)~=curPolygon(:,end)))
            curPolygon(:,end+1)=curPolygon(:,1);
        end
    end
    
    clippedPolygon=curPolygon;
end

function isInsideClipEdge=vertexIsInsideClipEdge(vertex,clipVertex1,clipEdge)
%%VERTEXISINSIDECLIPPEDEDGE Determine whether a vertex is inside of the
%clipping region given an edge. This relies on the fact that the clipping
%region is convex and goes in a counterclockwise direction. Thus, the sign
%of the angle with respect to the clipping edge is used to determine
%whether it is inside or outside. This could be done by extending the
%vectors to 3D (padding a zero on the ends) and then using a cross
%product, whereby the z component would tell us the sign, but the value of
%the z-component is the same as the determinant value used here.
%
%INPUTS:  vertex The 2X1 vertex to test whether it is on the correct size
%                of the edge to be in the clipping region.
%    clipVertex1 The first 2X1 vertex of the boundary used to form the
%                edge for clipping.
%       clipEdge The edge vector, which would be clipVertex2-clipVertex1.
%
%The code below is equivalent to
%isInsideClipEdge=det([vertex-clipVertex1,clipEdge])<=0;
%but should be slightly less susceptible to finite precision errors.

    S=[-clipEdge(1)*vertex(2);
        clipEdge(2)*vertex(1);
       -clipEdge(2)*clipVertex1(1);
        clipEdge(1)*clipVertex1(2)];

    isInsideClipEdge=exactSignOfSum(S)<=0;
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
