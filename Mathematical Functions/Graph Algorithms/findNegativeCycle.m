function cycleNodes=findNegativeCycle(adjMat)
%FINDNEGATIVECYCLE Given an adjacency matrix of a graph, find a negative
%                  cost cycle or indicate whether none exists. Unlike 
%                  simply running the function BellmanFordAlg, this
%                  algorithm will find a negative cost cycles even if it is
%                  not accessible from a particular source node. A cycle is
%                  a loop in the graph; a negative cost cycle is one where
%                  the sum of the costs in the cycle is negative.
%
%INPUTS: adjMat  An adjacency matrix full of costs for the graph.
%                adjMat(i,j) is the cost of going from vertex i to vertex
%                j. If there is no connection between the vertices in the
%                graph, then set that element to infinity. If the graph
%                might have negative weights along any paths, then
%                allowNegCosts should be true. Cycles accessible from the
%                source involving negative edge weights are detected by the
%                algorithm and mean that no solution exists.
%
%OUTPUTS:cycleNodes If no negative cost cycles exist in the graph, then 
%                this is an empty matrix. If any negative cycles exist, 
%                then this will be the nodes in order of the first negative
%                cycle found (this does not return all negative cycles
%                found. The last node in cycleNodes is corrected back to
%                the first node, finishing the loop.
%
%The algorithm is that of [1].
%
%As an example, consider
% adjMat=[Inf,    -20,    Inf,    Inf,    Inf;
%         Inf,    Inf,    1,      Inf,    Inf;
%         1,      Inf,    Inf,    1,      Inf;
%         Inf,    Inf,    Inf,    Inf,    1;
%         Inf,    Inf,    Inf,    Inf,    Inf];
% cycleNodes=findNegativeCycle(adjMat)
%There is a negative cycle around 1->2->3->1. However, if one used
%BellmanFordAlg starting at node 4, it would never find the cycle as nodes
%1, 2, and 3, are not accessible from node 4. On the other hand, this
%algorithm returns cycleNodes=[3;1;2];
%
%REFERENCES:
%[1] X. Huang, "Negative-weight cycle algorithms," in Proceedings of the
%    International Conference on Foundations of Computer Science, Las
%    Vegas, NV, 26-29 Jun. 2006, pp. 109-115. [Online].
%    Available: http://ww1.ucmss.com/books/LFS/CSREA2006/FCS4906.pdf
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numNodes=length(adjMat);

%Step 1: Augment the adjacency matrix with two new vertices. The first
%new vertex is connected to all existing nodes in the graph. All existing
%nodes in the graph are connected to the second new vertex. Both new
%vertices are appended to the end of adjMatNew.
adjMatNew=[[adjMat,inf(numNodes,1),ones(numNodes,1)];
          [ones(1,numNodes),Inf,Inf];
           inf(1,numNodes+2)];
source=numNodes+1;

%Step 2: Use the Bellman-Ford algorithm to find the cycles, if any. Since
%the extra nodes were placed at the end of the matrix, they will not affect
%the indices returned.
[~,~,cycleNodes]=BellmanFordAlg(adjMatNew,source);

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
