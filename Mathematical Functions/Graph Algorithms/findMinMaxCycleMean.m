function [lambda,cycle]=findMinMaxCycleMean(adjMat,findMin)
%%FINDMINMAXCYCLEMEAN Given the adjacency (cost) matrix for a graph, find
%                     the  minimum or maximum cycle mean. A cycle is a set
%                     of vertices in a directed graph that form a loop. The
%                     sum of the weights of the edges in a cycle divided by
%                     the number of vertices in the cycle is the cycle
%                     mean. This algorithm find the value of the minimum
%                     or maximum cycle mean and returns a cycle with that
%                     mean. If there are multiple cycles with the same
%                     mean, only one is returned.
%
%INPUTS: adjMat  An adjacency matrix full of costs for the graph.
%                adjMat(i,j) is the cost of going from vertex i to vertex
%                j. If there is no connection between the vertices in the
%                graph, then set that element to infinity. All valid costs
%                are finite, meaning that no costs should be set to -Inf.
%                it is assumed that no verted is connected back onto itself
%                (the diagonal of adjMat is Inf).
%        findMin An optimal parameter specify whether the minimum or the
%                maximum cycle mean is to be found. The default if omitted
%                is false, meaning that the maximum cycle mean is found.
%
%OUTPUTS: lambda The value of the maximum or minimum cycle mean as
%                specified from findMin. If a minimum/ maximum cycle mean
%                does not exist, then +/- infinity is returned, depending
%                on whether a maximization or minimization problem is being
%                performed.
%          cycle The indices of the vertices in the maximum or minimum
%                cycle mean. The first index is the same as the last index,
%                indicating that a full revolution of the cycle was made.
%                If a minimum/ maximum cycle mean does not exist, then an
%                empty matrix is returned.
%
%The algorithm for finding the minimum/ maximum cycle mean is that of Karp
%in [1]. The implementation is very similar to the description of Karp's
%algorithm given in [2]. The theoretically more efficient algorithms of [2]
%were not used as they are best suited for travering a sparse adjacency
%matrix and sparsity is not used by this function. Moreover, in [2], it is
%noted that Karp's algorithm is the fastest polynomial-time algorithm in
%the worst case.
%
%Consider the simple example
% adjMat=[Inf   -20   Inf;
%         Inf   Inf    10;
%           5   100   Inf];
% [lambda,cycle]=findMinMaxCycleMean(adjMat,false)
%One will get a maximum cycle mean of 55 and the corresponding cycle is
%cycle=[2;3;2].
%
%REFERENCES:
%[1] R. M. Karp, "A characterization of the minimum cycle mean in a
%    digraph," Discrete Mathematics, vol. 23, no. 3, pp. 309-311, 1978.
%[2] A. Dasdan and R. K. Gupta, "Faster maximum and minimum mean
%    cycle algorithms for system-performance analysis," IEEE Transactions
%    on Computer-Aided Design of Integrated Circuits and Systems, vol. 17,
%    no. 10, pp. 889-899, Oct. 1998.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    findMin=false;
end

%To perform minimization, just flip the sign of everything, but keep
%forbidden assignments forbidden.
if(findMin)
    sel=(adjMat==Inf);
    adjMat=-adjMat;
    %Keep forbidden assignments forbidden in the flipped cost matrix.
    adjMat(sel)=Inf;
end

%First, break the graph into strongly connected subgraphs
[setStartList,setLengths,setIdxList,subgraphList]=findStronglyConnectedSubgraphs(adjMat);

numSubgraphs=length(subgraphList);

%Next, run the algorithm on each strongly-connected subgraph. The
%solution is the maximum of all of the cycles
maxLambda=-Inf;
maxCycle=[];%This matters if there is no maximum cycle.
for curSubgraph=1:numSubgraphs
    adjMatCur=subgraphList{curSubgraph};
    
    [lambda,cycle]=maxCycleMeanKarp(adjMatCur);
    
    if(lambda>maxLambda)
        maxLambda=lambda;
        
        %Map the indices of the cycle from the subgraph back to those in
        %the original graph.
        idxMap=setIdxList(setStartList(curSubgraph):(setStartList(curSubgraph)+setLengths(curSubgraph)-1));
        maxCycle=idxMap(cycle);
    end
end

cycle=maxCycle;
if(findMin)
    lambda=-maxLambda;
else
    lambda=maxLambda;
end

end


function [lambda,cycle]=maxCycleMeanKarp(adjMat)
%This is Karp's algorithm from [1], as modified for a maximization in [2].
%It is here as a subroutine to make things more readable above as it must
%be called once for every strongly connected subgraph. 

n=length(adjMat);

%Arbitrarily set the source to the first node.
s=1;

%Head
D=-inf(n+1,n);
D(0+1,s)=0;%Distance to source is zero.
piMat=zeros(n+1,n);%Allocate sopace for the predecessor list.

%Body
for k=1:n%For each level
    for v=1:n%For each node
        %For each predecessor node
        for u=1:n
            %If u is not a predecessor of v
            if(adjMat(u,v)==Inf)
                continue;
            end
            
            temp=D(k-1+1,u)+adjMat(u,v);
            if(D(k+1,v)<temp)
                D(k+1,v)=temp;%max
                piMat(k+1,v)=u;
            end
        end
    end
end

%Now, D(k+1,v) holds the maximum weight of a path of length k from the
%source node to node v. If no such path exists, then D(k+1,v)=-Inf.

M=zeros(n,1);%Allocate space

%Tail
lambda=-Inf;
for v=1:n%For each node
    M(v)=Inf;%The identity for min
    for k=0:(n-1)
        temp=(D(n+1,v)-D(k+1,v))/(n-k);

        if(M(v)>temp)%min
            M(v)=temp;
        end
    end
    
    if(lambda<M(v))
        lambda=M(v);%max
        vOpt=v;
    end
end

%There is at least one maximum cycle mean critical cycle on the path from s
%to vOpt and its length is n-K(vOpt); pi is the predecessor list.

%If lambda is negative infinity, then there does not exist a maximum cycle
%mean critical path.
if(~isfinite(lambda))
    cycle=[];
    return;
end

%We reconstruct the entire path leading to node vOpt. Since we are
%scanning the nodes backwards, path is saved backwards.
path=zeros(n+1,1);
vCur=vOpt;
path(n+1)=vCur;
curPathIdx=1;
for k=n:-1:1
    curPathIdx=curPathIdx+1;
    vCur=piMat(k+1,vCur);
    path(n+1-curPathIdx+1)=vCur;
end

%A cycle exists along the path. We scan the path until we come across a
%node that is repeated. The repeated node form the beginning and ending of
%a path. We can then extract the cycle.

visitedIdx=zeros(n,1);
for curNodeIdx=1:(n+1)
    curNode=path(curNodeIdx);
    %If we have already visited this node, then we have found a complete
    %cycle.
    if(visitedIdx(curNode)~=0)
        cycle=path(visitedIdx(curNode):curNodeIdx);
        break;
    end
    visitedIdx(curNode)=curNodeIdx;
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
