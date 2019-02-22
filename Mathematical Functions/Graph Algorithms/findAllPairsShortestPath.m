function [distMat,pathMatrix]=findAllPairsShortestPath(adjMat)
%%FINDALLPAIRSSHORTESTPATH Given a graph represented by an adjacency
%                   matrix, find the shortest paths between all pairs of
%                   vertices in the graph. The Floyd-Warshall algorithm is
%                   used. If the distance between any node and itself is
%                   negative, then a negative cost cycle involving that
%                   node exists.
%
%INPUTS: adjMat  An adjacency matrix full of costs for the graph.
%                adjMat(i,j) is the cost of going from vertex i to vertex
%                j.
%
%OUTPUTS: distMat A matrix such that distMat(i,j) is the minimum distance
%                 from node i to node j. Values of dist(i,i) that are not
%                 infinite indicate cycles. Values of dist(i,i) that are
%                 negative indicate negative cycles and are not
%                 representative of the actual distance when making a
%                 single revolution of the cycle; they can be much larger.
%      pathMatrix A matrix that can be used to recreate the shortest path
%                 between any nodes. To go from node i to node j, the path
%                 is constructed as
%                     if(pathMatrix(i,j)==0)
%                         path=[];%No path exists between the nodes.
%                     else
%                         path=i;
%                         k=i;
%                         while(1)
%                             k=pathMatrix(k,j);
%                             path=[path;k];
%                             if(k==j)
%                                 break;%path now holds the path.
%                             end
%                         end
%                     end
%                  When there is a cycle, the path begins and ends
%                  with the same node.
%
%The Floyd-Warshall all-pairs shortest path algorithm is given in Chapter
%6.5 of [1] and Chapter 25.2 of [2]. The method of generating the path
%matrix is closer to that of [2] than that of [1]. As noted in [3], if no
%negative cycles are present in the graph, then the magnitude of all of
%the distances in distMat are less than N*max(max(abs(adjMat))), where N is
%the number of nodes (the length of adjMat). However, if a negative cost
%cycle is present, then it is shown that there exist graphs where values
%are 2*6^(N-1)*max(max(abs(adjMat))) during execution of the algorithm (do
%to going around the path multiple times), meaning that overflow issues can
%arise in certain very large scenarios. Of course, when using double
%precision arithmetic, that should not be a problem.
%
%The algorithm can find all nodes that are in cycles, but it will not find
%all cycles, because a single node can be in multiple cycles.
%
%REFERENCES:
%[1] C. H. Papadimitriou and K. Steiglitz, Combinatorial Optimization:
%    Algorithms and Complexity. Englewood Cliffs, NJ: Prentice-Hall Inc.,
%    1982.
%[2] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%[3] S. Hougardy, "The floyd-warshall algorithm on graphs with negative
%    cycles," Information Processing Letters, vol. 110, no. 8-9,
%    pp. 279-281, Apr. 2010.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of vertices.
numNodes=length(adjMat);

%Initialize the distance matrix.
distMat=adjMat;

%Initialize the path matrix
pathMatrix=zeros(numNodes,numNodes);

for i=1:numNodes
    for j=[1:(i-1),(i+1):numNodes]
        if(isfinite(adjMat(i,j)))
            pathMatrix(i,j)=j;
        end
    end
end

%These three loops are the part of the algorithm that finds the shortest
%distances.
for j=1:numNodes
    %The inner two loop is over all (i,k)~=j.
    for i=[1:(j-1),(j+1):numNodes]
        for k=[1:(j-1),(j+1):numNodes]
            %The if-statement performs
            %distMat(i,k)=min(distMat(i,k),distMat(i,j)+distMat(j,k));
            %and updates pathMatrix if
            %distMat(i,k)>distMat(i,j)+distMat(j,k).
            if(distMat(i,k)>distMat(i,j)+distMat(j,k))
                distMat(i,k)=distMat(i,j)+distMat(j,k);
                pathMatrix(i,k)=pathMatrix(i,j);
            end
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
