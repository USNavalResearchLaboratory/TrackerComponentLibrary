function isConvex=polygonIsConvex(vertices)
%%POLYGONISCONVEX Given an array of vertices that form a closed polygon,
%                 determine whether the polygon is convex.
%
%INPUTS: A 2XN vector of N subsequent vertices in the polygon. None of the
%        vertices may repeat (have a zero-length edge) with the exception
%        of the first vertex, which can optionally be repeated on the end.
%
%OUTPUTS: isConvex A boolean variable. True if the polygon is convex, false
%                  is the polygon is not convex.
%
%A polygon is convex if no internal angle is more than 180 degrees. This
%means that the polygon cannot turn back on itself.If it turns back on
%itself, then the direction of the angle between subsequent vectors will
%change. The function turnOrientation is used with all sets of three
%vertices (all pairs of edges) to determine whether the sign ever changes.
%If the sign never changes going around the polygon, then it is convex.
%
%More on convex polygons can be found at [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Convex Polygon." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ConvexPolygon.html
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%We cannot have duplicate adjacent vertices, except the first vertex can
%optionally be repeated at the end. In such an instance, we test for it and
%get rid of it.
    if(all(vertices(:,1)==vertices(:,end)))
        vertices=vertices(:,1:(end-1));
    end

    numVertices=size(vertices,2);

    if(numVertices<3)
        isConvex=false;
        return;
    end
    
    P1=vertices(:,1);
    P2=vertices(:,2);
    P3=vertices(:,3);
    
    sense=turnOrientation(P1,P2,P3);
    P1=P2;
    P2=P3;
    for curVertex=4:numVertices
        P3=vertices(:,curVertex);
        
        curSense=turnOrientation(P1,P2,P3);
        P1=P2;
        P2=P3;
        
        %If the first node was repeated, so a sense had not been
        %established.
        if(sense==0&&curSense~=0)
            sense=curSense;
            continue; 
        end
        
        if(curSense~=0&&curSense~=sense)
            isConvex=false;
            return; 
        end
    end
    
    %Deal with the last two angles, which contain both the end point as
    %well as the points at the beginning.
    P3=vertices(:,1);
    curSense=turnOrientation(P1,P2,P3);
    P1=P2;
    P2=P3;
    if(curSense~=0&&curSense~=sense)
        isConvex=false;
        return; 
    end
    
    %The final point to test, which starts with vertices(:,numVertices);
    P3=vertices(:,2);
    curSense=turnOrientation(P1,P2,P3);
    if(curSense~=0&&curSense~=sense)
        isConvex=false;
        return; 
    end
    
    isConvex=true;
    return;
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
    