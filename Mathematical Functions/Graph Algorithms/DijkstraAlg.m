function [retPath,retDist]=DijkstraAlg(adjMat,sourceIdx,destIdx,allowNegCosts)
%%DIJKSTRAALG  Use Dijkstra's algorithm to find the shortest path from a
%              given source node through a graph represented using an
%              adjacency matrix. Unlike traditional implementations, the
%              option of allowing possible negative costs exists.
%
%INPUTS: adjMat  An adjacency matrix full of costs for the graph.
%                adjMat(i,j) is the cost of going from vertex i to vertex
%                j. If there is no connection between the vertices in the
%                graph, then set that element to infinity. If the graph
%                might have negative weights along any paths, then
%                allowNegCosts should be true. Cycles accessible from the
%                source involving negative edge weights are detected by the
%                algorithm and mean that no solution exists.
%      sourceIdx The index of the node (vertex) from which the shortest
%                paths to other vertices is desired.
%        destIdx The index of the destination path. If this parameter is
%                provided and is not the empty matrix, then the path of
%                nodes from the source to the destination is returned along
%                with the path length. Otherwise, information allowing one
%                to reconstruct the whole path is returned if this
%                parameter is omtted.
%  allowNegCosts A boolean variable indicating whether an algorithm
%                allowing negative values in adjMat should be used. The
%                default if this parameter is omitted is true. Allowing
%                negative edge costs is slower. If this is set to false and
%                negative edge costs are present, then the algorithm will
%                return incorrect results.
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
%
%The algorithms implemented here are described in Chapter 9.3 of [1].
%Specifically, the algorithm given no negative costs is Chapter 9.3.2 and
%the algorithm for the case of possibly negative costs is 9.3.3.
%
%Solving the shortest augmenting path algorithm  using Dijkstra's algorithm
%is also discussed in Chapter 24.3 of [2].
%
%If in addition to detecting the presence of negative cycles accessible
%from the source, one wishes to find the vertices in the cycle, then the
%function BellmanFordAlg should be used instead of this one.
%
%The algorithm can be demonstrated on the same problem given in Figure 8 of
%[3] for the Viterbi algorithm, as the Dijkstra algorithm is more general
%than the Viterbi algorithm. In this instance, the nodes in Figure 8 are
%numbered increasing for each column, down, from left to right. Thus, the
%adjacency matrix and the shortest path from nodes 1 to 14 are given by the
%commands
% adjMat=inf(14,14);
% adjMat(1,2)=1;
% adjMat(1,3)=1;
% adjMat(2,4)=1;
% adjMat(2,6)=1;
% adjMat(3,5)=2;
% adjMat(3,7)=0;
% adjMat(4,8)=0;
% adjMat(4,10)=2;
% adjMat(5,8)=2;
% adjMat(5,10)=0;
% adjMat(6,9)=1;
% adjMat(6,11)=1;
% adjMat(7,9)=1;
% adjMat(7,11)=1;
% adjMat(8,12)=2;
% adjMat(9,12)=0;
% adjMat(10,13)=1;
% adjMat(11,13)=1;
% adjMat(12,14)=1;
% adjMat(13,14)=1;
% [path,dist]=DijkstraAlg(adjMat,1,14);
%where the answer should be path=[1;3;7;9;12;14]; and dist=3;
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%[2] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%[3] G. D. Forney Jr., "The Viterbi algorithm," Proceedings of the IEEE,
%    vol. 61, no. 3, pp. 268-278, Mar. 1973.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    allowNegCosts=true;
end

%The number of vertices.
numNodes=length(adjMat);

%Create a maximum heap for the indices.
nodeHeap=BinaryHeapOfIndices(numNodes,false);

%The shortest distance yet found from the source to each of the nodes.
dist=inf(numNodes,1);
%The distance from the source to itself is zero. This is the first node to
%be visited.
dist(sourceIdx)=0;
%The previous node in the optimal path from the given node to the source.
prevNodes=zeros(numNodes,1);
%Insert the source node into the heap.
nodeHeap.insert(dist(sourceIdx),sourceIdx);

if(allowNegCosts==false)
    visitedNodes=false(numNodes,1);
    %If all of the values in adjMat are positive.
    while(~nodeHeap.isEmpty())
        %Get the least-cost vertex that has not been visited yet.
        topEl=nodeHeap.deleteTop();
        dist2Min=topEl.key;
        nodeIdxMin=topEl.value;

        %Mark the node as heving been visited.
        visitedNodes(nodeIdxMin)=true;
        %Now, go through the nodes that connect to this one.
        for destNode=1:numNodes
            %If the node has already been visited, then skip it.
            if(visitedNodes(destNode)==true)
                continue;
            end

            alt=dist2Min+adjMat(nodeIdxMin,destNode);
            %If a shorter path to the destNode node exists through nodeIdxMin,
            %than what is already set, then switch to that.
            if(alt<dist(destNode))
                dist(destNode)=alt;
                prevNodes(destNode)=nodeIdxMin;

                %If the node is not already in the heap, the add it.
                %Otherwise, change its value.
                if(nodeHeap.indexIsInHeap(destNode)==false)
                    nodeHeap.insert(dist(destNode),destNode);
                else
                    nodeHeap.changeIndexedKey(dist(destNode),destNode);
                end
            end
        end
    end
else%Some edges can have negative costs.
    %This is to detect negative cycles so the algorithm does not iterate
    %forever.
    numTimesVisited=zeros(numNodes,1);
    
    while(~nodeHeap.isEmpty())
        %Get the least-cost vertex that has not been visited yet.
        topEl=nodeHeap.deleteTop();
        dist2Min=topEl.key;
        nodeIdxMin=topEl.value;

        numTimesVisited(nodeIdxMin)=numTimesVisited(nodeIdxMin)+1;
        
        %If the graph has a negative cycle.
        if(numTimesVisited(nodeIdxMin)>numNodes)
            retPath=[];
            retDist=[];
            return;
        end
        
        %Now, go through the nodes that connect to this one.
        for destNode=1:numNodes
            alt=dist2Min+adjMat(nodeIdxMin,destNode);
            %If a shorter path to the destNode node exists through 
            %nodeIdxMin, than what is already set, then switch to that.
            if(alt<dist(destNode))
                dist(destNode)=alt;
                prevNodes(destNode)=nodeIdxMin;

                %If the node is not already in the heap, the add it. 
                %Otherwise, change its value.
                if(nodeHeap.indexIsInHeap(destNode)==false)
                    nodeHeap.insert(dist(destNode),destNode);
                else
                    nodeHeap.changeIndexedKey(dist(destNode),destNode);
                end
            end
        end
    end
end

%If a destination index is given, the recreate the path.
if(nargin>2&&~isempty(destIdx))
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

%Get rid of the heap.
nodeHeap.delete();

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
