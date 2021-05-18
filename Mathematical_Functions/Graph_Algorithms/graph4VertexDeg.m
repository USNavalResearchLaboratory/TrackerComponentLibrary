function [edges,didSucceed]=graph4VertexDeg(d)
%%GRAPH4VERTEXDEG Given a list of the degrees of vertices in a
%           non-directional graph, get a set of edges representing a graph
%           that satisfies the degree profile, or indicate if no such graph
%           can exist. Note that a set of edges for a degree profile is not
%           unique. This just returns one set.
%
%INPUTS: d An nX1 or 1Xn vector of the degrees of the n vertices. These are
%          integer values >=0. They represent the number of edges touching
%          the nodes.
%
%OUTPUTS: edges A 2XnumEdge set of edges that satisfies the degree profile.
%               An empty matrix is returned if all edges have zero degree
%               or if no set ofvedges exists that can result in the given
%               degree profile. An edge [i;j] goes between vertex i and
%               vertex j and adds to the degree count for both vertices.
%    didSucceed This is true if an edges exists that can satisfy d.
%               Otherwise, this is false and edges is an empty matrix.
%
%This function implements Algorithm H in Chapter 7 of [1].
%
%EXAMPLE:
% [edges,didSucceed]=graph4VertexDeg([2;1;6;3;2;1;1])
%The algorithm in the example will succeed and the edge set is
%edges=[7,6,2,5,5,1,1,4;
%       3,3,3,4,3,4,3,3];
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[d,idx]=sort([d(:);0],'descend');

%Allocate space
maxEdges=ceil(sum(d)/2);
edges=zeros(2,maxEdges);
c=zeros(d(1),1);

numEdges=0;

%Step H1
k=d(1);
j=0;
while(k>0)
    j=j+1;
    while(k>d(j+1))
        c(k)=j;
        k=k-1;
    end
end
%If all of the d's are zero.
if(j==0)
    edges=[];
    didSucceed=true;
    return;
end

while(1)
    %Step H2, find n.
    n=c(1);
    if(n==0)
        %Shrink to fit and restore the original indexation.
        edges=idx(edges(:,1:numEdges));
        %Restore the original indexation.
        didSucceed=true;
        return;
    elseif(d(1)>=n)
        edges=[];
        didSucceed=false;
        return;
    end

    %Step H3
    i=1;
    t=d(1);
    r=c(t);
    j=d(n);
    
    %Step H4, generate a new edge.
    while(1)
        c(j)=c(j)-1;
        m=c(t);
        numEdges=numEdges+1;
        edges(:,numEdges)=[n;m];
        d(m)=d(m)-1;
        c(t)=m-1;
        j=j-1;
        if(j==0)
            break;
        elseif(m==i)
            i=r+1;
            t=d(i);
            r=c(t);
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
