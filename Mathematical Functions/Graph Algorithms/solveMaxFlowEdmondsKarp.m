function [maxFlow,F]=solveMaxFlowEdmondsKarp(CMat,source,sink)
%%SOLVEMAXFLOWEDMONDSKARP Given a graph of vertices (nodes) with capacity
%                         limitations between nodes, find the maximum flow
%                         that can pass between the source and sink nodes.
%                         The connectivity of the nodes is specified by an
%                         adjacency matrix CMat. The Edmonds-Karp algorithm
%                         is used.
%
%INPUTS: CMat An NXN adjacency matrix for N nodes. CMat(i,j) is the amount
%             of flow that can go from node i to node j. If the nodes are
%             not connected in the graph, then CMat(i,j)=0. The values of
%             the diagonal elements (nodes connected with themselves) are 
%             ignored.
%      source The index of the source node.
%        sink The index of the sink node.
%
%OUTPUTS: maxFlow The maximum amount of flow that can be pushed from the
%                 source to the sink in the graph.
%           F     A matrix indicating how much flow is going through the
%                 nodes in the maximum flow solution. F(i,j) is the amount
%                 of flow going from node i to node j. Note that
%                 F(i,j)=-F(j,i) as amount of flow going from node i to
%                 node j is the negative of that in the opposite direction.
%
%The Edmonds-Karp algorithm, which is a modification of the Ford-Fulkerson
%algorithm, is described in Chapter 26.2 of [1].
%
%The maximum flow problem maximizes sum(F(source,:)) subject to the
%constraints:
%1) Flow Conservation: sum(F(u,:))==0 for all u that are not the source or
%                      sink.
%2) Capacity Constraint: 0<=F(u,v)<=CMat(u,v)
%3) Skew Symmetry: F(u,v)=-F(v,u)
%
%The algorithm can be demonstrated using the example from Figure 26.1 of
%the book:
% source=1;
% sink=6;
% CMat=[0, 16, 13, 0,  0, 0;
%       0,  0, 10, 12,  0, 0;
%       0,  4,  0, 0, 14, 0;
%       0,  0,  9, 0,  0, 20;
%       0,  0,  0, 7,  0, 4;
%       0,  0,  0, 0,  0, 0];
% [maxFlow,F]=solveMaxFlowEdmondsKarp(CMat,source,sink);
%One will find the maximum flow to be 23 and F is 
% F =
%      0    12    11     0     0     0
%    -12     0     0    12     0     0
%    -11     0     0     0    11     0
%      0   -12     0     0    -7    19
%      0     0   -11     7     0     4
%      0     0     0   -19    -4     0
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The total number of vertices in the network.
n=length(CMat);

%Initial flow magnitude is zero.
maxFlow=0;
%This matrix will hold the flow values between vertices that have been
%established. F(v1,v2) is the flow between vertex A and vertex B on the
%established path. 
F=zeros(n,n);

while(1)
    %Given a flow network F having flow f, there is a residual network,
    %which consists of edges whose flows can be increased. The residual
    %capacity of an edge from vertex 1 to vertex 2 is
    %ARes(v1,v2)=CMat(v1,v2)-F(v1,v2). That is the extra amount of flow
    %that might be possible to push through the edge before exceeding the
    %capacity. The residual network is the set of edges such that
    %ARes(v1,v2)>0. The residual network is itself a flow problem with a
    %smaller CMat.
    
    %An augmenting path is a path from source to sink in the residual
    %network. Each edge on the augmenting path allows for more flow without
    %violating capacity constraints.
    
    %The Fork-Fulkerson algorithm finds some augmenting path and then
    %increases the flow on each edge in F by the residual capacity. The
    %Edmonds-Karp algorithm chooses the augmenting path to be the
    %augmenting path that uses the fewest nodes, which can be found using a
    %breadth-first search. They proved that choosing the shortest (in terms
    %of number of nodes) augmenting path lets the algorithm terminate in
    %polynomial time. In comparison, the Ford-Fulkerson algorithm can get
    %into an infinite loop. The algorithm terminates when there is no more
    %augmenting path. The last augmenting path found is thus the maximum
    %flow path.
    
    [addedFlow,augPath]=shortestAugmentingPath(CMat,source,sink,F);
    
    %When the capacity of the returned path is zero, then the problem has
    %been solved.
    if(addedFlow==0)
        break;
    end
    
    maxFlow=maxFlow+addedFlow;
    
    %The total flow on each edge of the path increases by the residual
    %capacity. The path is scanned backwards from sink to source.
    v=sink;
    while(v~=source)
        u=augPath(v);
        F(u,v)=F(u,v)+addedFlow;
        F(v,u)=F(v,u)-addedFlow;
        v=u;
    end
end
end

function [pathFlow,thePath]=shortestAugmentingPath(CMat,source,sink,F)
%Breadth-first search for the shortest augmenting path.

n=length(CMat);
deltaFlow=zeros(n,1);
deltaFlow(source)=Inf;

%Mark all of the vertices as unvisited.
thePath=-1*ones(n,1);
%This ensures that the path never returns to the source. This also means
%that the path will have to be scanned in reverse.
thePath(source)=-2;
%Initialize a queue with the source in it.
theQueue=queue(source);

while(~theQueue.isEmpty)
    %Get the top node in the queue.
    u=theQueue.dequeueTop().Data;
    
    %Search all of the neighboring connections that have a positive flow
    %capacity.
    for v=1:n
        %The path must go somewhere that has not already been visited.
        if(v==u||thePath(v)~=-1)
            continue;
        end
        
        %Find the residual capacity.
        resCap=CMat(u,v)-F(u,v);
        if(resCap>0)
            %Connect u and v.
            thePath(v)=u;
            %The change in flow into this node.
            deltaFlow(v)=min(deltaFlow(u),resCap);
            if(v==sink)
                %If a path connecting the source and the sink has been
                %formed, then we have found the shortest augmenting path. 
                pathFlow=deltaFlow(sink);
                return;
            else
                %Add the node to the queue to possibly revisit later.
                theQueue.addToTop(v);
            end
        end
    end
end

%If we get here, then there is no augmenting path.
pathFlow=0;
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
