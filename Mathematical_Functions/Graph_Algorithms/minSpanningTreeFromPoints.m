function treeEdgeIdx=minSpanningTreeFromPoints(distMat)
%MINSPANNINGTREEFROMPOINTS Given a set of pariwise distances between
%             points, find the minimum spanning tree between the points.
%             This is a set of lines connecting all of the points
%             such that the sum of the lengths of the lines is minimized.
%             Such a tree is useful, for example, when planning how to lay
%             power lines to connect a series of cities. This algorithm is
%             O(n^2) due to the use of the pairwise distances.
%
%INPUTS: distMat An NXN matrix of positive pairwise distances between
%                points. The matrix is assumed to be symmetric.
%
%OUTPUTS: treeEdgeIdx A set of indices such that the ith edge of the
%                     minimum spanning tree is (i, treeEdgeIdx(i)), where i
%                     refers to the ith point and treeEdgeIdx(i) provides
%                     the index of the point to which the edge connects.
%
%The algorithm is based on MINSPT from [1] and is essentially an
%implementation of Prim's algorithm for dense graphs. When given sparse
%graphs, Prim's algorithm using a priority queue, as described in Chapter
%23.2 of [2] is more efficient.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%[2] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of points
N=size(distMat,1);
NM1=N-1;

%Allocate space and initialize
treeEdgeIdx=(-N)*ones(N,1);
treeEdgeIdx(N)=0;

for l=1:NM1
    dMin=Inf;
    for i=1:NM1
        it=treeEdgeIdx(i);
        
        if(it<0&&i~=it)%If point is not already in the tree.
            d=distMat(-it,i);
 
            if(d<dMin)
                dMin=d;
                iMin=i;
            end
        end
    end%Line 41
    
    %Adjoin new edge
    treeEdgeIdx(iMin)=-treeEdgeIdx(iMin);
    
    for i=1:NM1
        it=treeEdgeIdx(i);
        if(it>0)
            continue;
        end
        d1=distMat(i,iMin);
        d2=distMat(i,-it);
 
        if(d1<d2)
           treeEdgeIdx(i)=-iMin; 
        end
    end
end

%The minimum spanning tree has N-1 edges. The final edge is just the last
%point back to itself and does not count.
treeEdgeIdx=treeEdgeIdx(1:NM1);

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
