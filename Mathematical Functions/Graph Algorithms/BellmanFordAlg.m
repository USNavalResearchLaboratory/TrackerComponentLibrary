function [retPath,retDist,cycleNodes]=BellmanFordAlg(adjMat,sourceIdx,destIdx)
%%BELLMANFORDALG Use the Bellman-Ford algorithm to find the shortest path
%                from a given source node through a graph represented using
%                an adjacency matrix. Negative edge costs are allowed and
%                if a negative cycle reachable from the source exists in
%                the graph, then that cycle will be returned (and no
%                minimum distance path will exist).
%
%INPUTS: adjMat  An adjacency matrix full of costs for the graph.
%                adjMat(i,j) is the cost of going from vertex i to vertex
%                j. If there is no connection between the vertices in the
%                graph, then set that element to infinity. Cycles
%                accessible from the source involving negative edge weights
%                are detected by the algorithm and mean that no solution
%                exists.
%     sourceIdx  The index of the node (vertex) from which the shortest
%                paths to other vertices is desired.
%        destIdx The index of the destination path. If this parameter is
%                provided, then the path of nodes from the source to the
%                destination is returned along with the path length.
%                Otherwise, information allowing one to reconstruct the
%                whole path is returned if this parameter is omtted.
%       
%OUTPUTS: retPath If destIdx is given, this is the minimum cost sequence of
%                 nodes to get from sourceIdx to destIdx, including the end
%                 nodes; if there is no path to the node, an empty matrix
%                 is returned. Otherwise, this is prevNodes, such that one
%                 can recreate all shortest paths. Starting at index idx,
%                 the next node on the path in the reverse direction from
%                 destimation to source is retPath(idx). By looping, one
%                 can recreate the whole path. If there is no path to a
%                 node, then retPath(idx)=0. If a negative cycle exists
%                 in the path, then an empty matrix is returned. 
%         retDist If destIdx is given, this is the shortest distance from
%                 the source to the destination. If destIdx is not given,
%                 then this is an array where each index holds the shortest
%                 distance from the source to that node. If a negative
%                 cycle exists in the graph, then an empty matrix is
%                 returned.
%      cycleNodes If no negative cost cycles exist in the graph, then this
%                 is an empty matrix. If any negative cycles exist, then
%                 this will be the nodes in order of the first negative
%                 cycle found. The last node in cycleNodes is connected
%                 back to the first node, finishing the loop. A cycle is a
%                 loop in the graph; a negative cost cycle is one where the
%                 sum of the costs in the cycle is negative.
%
%The algorithm implemented is that of Chapter 24.1 of [1], where the
%extraction of a negative cost cycle accessible from the source is added
%when a negative cost cycle is detected.
%
%Ad an example, consider the adjacency matric with two negative cost
%cycles:
% adjMat=[Inf,    -20,    Inf;
%         Inf,    Inf,    10,
%         5,      -15,    Inf];
% [retPath,retDist,cycleNodes]=BellmanFordAlg(adjMat,1,3);
%Looking for a path from 1 to 3, the algorithm would encounter the cycle.
%Thus, the return values would be [] and [] for retPath and retDist and 
%cycleNodes=[3;2], because that is the first cycle found.
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of vertices.
numNodes=length(adjMat);

%The shortest distance yet found from the source to each of the nodes.
dist=inf(numNodes,1);
%The distance from the source to itself is zero. This is the first node to
%be visited.
dist(sourceIdx)=0;
%The previous node in the optimal path from the given node to the source.
prevNodes=zeros(numNodes,1);

for i=1:(numNodes-1)
    for curNode=1:numNodes
        %We have to go over all edges from the current node. If
        %adjMat(curNode,destNode) is infinite, then no edge exists to destNode.
        for destNode=1:numNodes
            w=adjMat(curNode,destNode);%The edge weight
            if(~isfinite(w))
                continue;
            end

            %Perform the relaxation step.
            if(dist(curNode)+w<=dist(destNode))
                dist(destNode)=dist(curNode)+w;
                prevNodes(destNode)=curNode;
            end
        end
    end
end

%Now, we check for a negative weight cycle (reachable from the source). If
%a negative weight cycle is found, then return it. We have to look at all
%edges in the graph.
cycleStart=[];
for curNode=1:numNodes
    for destNode=1:numNodes
        w=adjMat(curNode,destNode);%The edge weight
        if(~isfinite(w))
            continue;
        end
        
        if(dist(curNode)+w<dist(destNode))
            cycleStart=destNode;
            prevNodes(destNode)=curNode;
            break;
        end 
    end
end

%If a negative weight cycle is present, then trace back along it to get the
%cycle to return. However, a problem is that from the starting node, we
%might get lost in a cycle that does not return to the starting node. Thus,
%we have to keep track of when we visited each edge and declare a cycle
%once we revisit any edge, not necessarily the starting edge.

if(~isempty(cycleStart))
    %Allocate the maximum space possible for the cycle.
    cycleNodes=zeros(numNodes,1);
    nodeVisited=zeros(numNodes,1);
    
    curIdx=1;
    cycleNodes(curIdx)=cycleStart;
    nodeVisited(cycleStart)=curIdx;
    while(1)
        nextNode=prevNodes(cycleNodes(curIdx));
        %if we have seen this node before
        if(nodeVisited(nextNode)~=0)
            break;
        end
        cycleNodes(curIdx+1)=prevNodes(cycleNodes(curIdx));
        nodeVisited(nextNode)=curIdx+1;
        
        curIdx=curIdx+1;
    end 
    
    %Shrink the list of cycle nodes down to the size of the number of
    %actual nodes in the cycle and reverse the order so that the nodes are
    %in the correct order.
    cycleNodes=flipud(cycleNodes(nodeVisited(nextNode):curIdx));
    
    %It cannot solve for the optimal path, because of the negative cycle.
    retPath=[];
    retDist=[];
    return;
else
    cycleNodes=[];
end

%If a destination index is given, the recreate the path, if there are no
%cycles.
if(nargin>2)
    %If there is no path to the desired node.
    if(~isfinite(dist(destIdx)))
        retPath=[];
        retDist=Inf;
        return;
    end

    path=zeros(numNodes,1);
    path(1)=destIdx;
    pathStep=1;
    while(path(pathStep)~=sourceIdx)
        path(pathStep+1)=prevNodes(path(pathStep));
        pathStep=pathStep+1;
    end
    retPath=flipud(path(1:pathStep));
    retDist=dist(destIdx);
else
    retPath=prevNodes;
    retDist=dist;
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
