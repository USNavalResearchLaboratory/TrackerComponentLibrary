function [F,totalCost,exitCode]=minCostFlow(AMat,CMat,b,maxIter)
%%MINCOSTFLOW Solve the minimum cost flow problem using a strong polynomial
%             time cycle cancelling algorithm that works with non-integer
%             costs. Minimum cost flow problems include transportation
%             problems.
%
%INPUTS: AMat An NXN matrix of costs in a directed graph such that
%             AMat(i,j) is the cost of an edge going from vertex i to
%             vertex j. If no vertex goes from node i to node j, then a
%             cost of 0 should be inserted.
%        CMat An NXN matrix of capacities of the edges in the graph.
%             CMat(i,j) is the capacity of an edge from node i to node j.
%             If no edge exists, then a capacity of zero should be used. It
%             is assumed that all capacities are finite.
%           b An NX1 vector of the supply provided by each node. It is
%             required that sum(b)=0.
%     maxIter An optional parameter specifying the maximum number of
%             iterations that should be performed by the algorithm. If
%             omitted, a maximum of 100+10N iterations is used.
%
%OUTPUTS: F The flow matrix. F(i,j) is the amount of flow going from node i
%           to node j. Note that F(i,j)=-F(j,i). If the algorithm could not
%           obtain an initial feasible solution, then an empty matrix is
%           returned. If the maximum number of iterations is exceeded, then
%           the algorithm will return a feasible solution that is not
%           optimal.
% totalCost The total cost of the flow. This is the quantity being
%           minimized. If the algorithm could not obtain an initial
%           feasible solution, then an empty matrix is returned.
%  exitCode A parameter indicating whether an error occurred or whether the
%           algorithm terminated successfully. Possible values are:
%           0 The algorithm was successful.
%           1 No feasible solution could be found.
%           2 The maximum number of iterations was reached.
%           3 Inputs are invalid. This means that either C contains NaN
%             values, sum(b)~=0 or AMat contains nonfinite values.
%
%The minimum cost flow problems seeks to find a flow F to minimize
%sum(sum(AMat.*F)) subject to the constraints:
%4) Total Vertex Flow: sum(F(i,:))=b(i) for all i.
%2) Capacity Constraint: 0<=F(u,v)<=CMat(u,v)
%3) Skew Symmetry: F(u,v)=-F(v,u)
%Whereas the maximum flow problem, which is solved in the function
%solveMaxFlowEdmondsKarp seeks to maximize the flow between two nodes the
%minimum cost flow problem seeks to minimize the total cost of a fixed
%amount of flow.
%
%The strong polynomial cycle cancelling algorithm of [1] is used. However,
%the epsilon-based convergence acceleration method described in [1] is not
%used. Rather, the algorithm is more similar to the simple cycle
%cancellation algorithm of Chapter 7.2 of [2], except the minimum mean
%cycle is always augmented rather than any negative cost cycle.
%
%The complexity of the non-epsilon accelerated algorithm in [1] is bounded
%as O(N^2*m^3*log(N)), where N is the number of vertices in the graph, and
%m is the number of edges. The true complexity depends on the efficiency of
%the function computeResidualCapacity, which is used to obtain the negative
%cycles.
%
%EXAMPLE: This is example 7.1 in Chapter 7.2 of [2].
% AMat=[0, 4, 1, 0;
%       0, 0, 2, 5;
%       0, 3, 0, 2;
%       0, 0, 0, 0];
% CMat=[0, 2, 2, 0;
%       0, 0, 1, 1;
%       0, 1, 0, 1;
%       0, 0, 0, 0];%Capacity matrix
% b=[2;0;0;-2];
% [F,totalCost,exitCode]=minCostFlow(AMat,CMat,b)
%One will get a flow matrix of 
% F =[ 0     0     2     0;
%      0     0    -1     1;
%     -2     1     0     1;
%      0    -1    -1     0];
%having a total cost of 12.
%
%OTHER EXAMPLES: Transportation problems can be transformed into minimum
%cost flow problems and vice versa, as shown in Chapters 7.4 and 7.5 of
%[2]. Thus, the examples in the comments to the function
%solveTransportationProblem can be taken as additional examples of the
%minimum cost flow problem.
%
%REFERENCES:
%[1] A. V. Goldberg and R. E. Tarjan, "Finding minimum-cost circulations
%    by canceling negative cycles," Journal of the Association for
%    Computing Machinery, vol. 36, no. 4, pp. 873-886, Oct. 1989.
%[2] C. H. Papadimitriou and K. Steiglitz, Combinatorial Optimization:
%    Algorithms and Complexity. Englewood Cliffs, NJ: Prentice-Hall Inc.,
%    1982.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVertex=length(b);

if(nargin<5||isempty(maxIter))
    maxIter=100+10*numVertex;
end

if(any(isnan(CMat(:)))||sum(b)~=0||any(~isfinite(AMat(:))))
    exitCode=3;%The problem is infeasible. 
    F=[];
    totalCost=[];
    return;
end

%Step 1: Find a feasible flow.
%In this step, we must find a feasible flow F. This is a flow that
%satisfies constraints 1-3 but does not necessarily minimize the cost
%sum(sum(AMat.*F)). This can be done by ignoring costs and solving a
%maximum flow problem based only on the constraints. The maximum flow
%problem does not have demand constraints on all of the nodes, rather it
%only has them on a sink and source which are added.
%
%To enforce constraint 1, we have to create a network with an artificial
%sink and source. If b(i) is positive, then we need an arc from the source
%to that node with capacity b(i). If b(i) is negative, then we need an arc
%from the sink to that node with capacity -b(i). If the maximum flow
%problem solved on the modified graph saturates all of the edges coming out
%of the source, then the problem is feasible.

%Augment the capacity matrix.
CMatAugment=[CMat,zeros(numVertex,2);
             zeros(2,numVertex+2)];

source=numVertex+1;
sink=numVertex+2;
totalFlowIn=0;
for curVertex=1:numVertex
    if(b(curVertex)>0)
        CMatAugment(source,curVertex)=b(curVertex);
        totalFlowIn=totalFlowIn+b(curVertex);
    else
        CMatAugment(curVertex,sink)=-b(curVertex);
    end
end

[maxFlow,F]=solveMaxFlowEdmondsKarp(CMatAugment,source,sink);

%Check whether all of the source nodes have been saturated. This is true if
%the total flow in equals the maximum flow found. A comparison to an
%epsilon value is used as it is not clear whether it can be guaranteed
%within finite precision bounds that the two quantities will always be
%equal for feasible problems, though they seem to generally be equal for
%feasible problems.
if(abs(maxFlow-totalFlowIn)>eps(totalFlowIn))
    exitCode=1;%The problem is infeasible. 
    F=[];
    totalCost=[];
    return;
end

%Get rid of the artificial source and sink nodes from the flow matrix.
F=F(1:numVertex,1:numVertex);

%To make the computation of the residual matrix simpler, we will transform
%the cost matrix. This simplifies the computation of costs in the residual
%flow network.
AMatOrig=AMat;
AMat=AMat-AMat';

%Compute the residual flow network, also known as the incremental flow
%network. Definition 7.3 from Chapter 7.2 of [2] explains what the
%incremental flow network is and how the costs are assigned. The residual
%capacity is updated during augmentation of the flow matrix F, so that it
%need not be recomputed each loop.

%Find the capacity of the residual flow network.
residualCapacity=CMat-F;
for curIter=1:maxIter
    %Step 2: Find the minimum mean cost cycle with respect to the costs on
    %        the residual flow network, also known as the incremental flow
    %        network.
    
    %The costs associated with the residual flow network are AMat(i,j) on
    %forward nodes and -AMat(i,j) on backward nodes. Hence the reason why
    %we replaced AMat with AMat-AMat'.
    costMat=Inf(numVertex,numVertex);
    costMat(residualCapacity~=0)=AMat(residualCapacity~=0);
        
    [minCycleMean,cycleVertices]=findMinMaxCycleMean(costMat,true);
    
    if(minCycleMean>=0)
        exitCode=0;
        FP=F;
        FP(F<0)=0;
        totalCost=sum(sum(FP.*AMatOrig));
        return;
    end
    
    %Step 3: The minimum mean cycle is negative. We must cancel the minimum
    %        mean cycle. That is, we augment the graph along the negative
    %        cycle. We adjust the flow around the cycle by as much
    %        as possible without violating capacity constraints, so that
    %        the negative cycle no longer exists. 
    
    %Determine the minimum residual capacity of the cycle.
    numVertexInCycle=length(cycleVertices);
    minCapacity=Inf;
    startVertex=cycleVertices(1);
    for curEndVertex=2:numVertexInCycle
        endVertex=cycleVertices(curEndVertex);
        
        curCapacity=residualCapacity(startVertex,endVertex);
        minCapacity=min(curCapacity,minCapacity);
        
        startVertex=endVertex;
    end
    
    %Push the amount of flow around the negative cycle that will saturate
    %the minimum capacity arc on the minimum mean cycle. Note that the skew
    %symmetry constraint dictates how the flow is to be added (how
    %augmentation is performed).
    startVertex=cycleVertices(1);
    for curEndVertex=2:numVertexInCycle
        endVertex=cycleVertices(curEndVertex);
        
        %Update the flow
        F(startVertex,endVertex)=F(startVertex,endVertex)+minCapacity;
        F(endVertex,startVertex)=F(endVertex,startVertex)-minCapacity;

        %Update the residual flow network. The update is opposite that of
        %F, because the residual flow is CMat-F.
        residualCapacity(startVertex,endVertex)=residualCapacity(startVertex,endVertex)-minCapacity;
        residualCapacity(endVertex,startVertex)=residualCapacity(endVertex,startVertex)+minCapacity;
   
        startVertex=endVertex;
    end
end

%Find the cost of the flow.
FP=F;
FP(F<0)=0;
totalCost=sum(sum(FP.*AMatOrig));

%If we are here, the maximum number of iterations has elapsed. We will
%check whether convergence occurred on the final iteration, and it not,
%then we will return accordingly.
costMat=Inf(numVertex,numVertex);
costMat(residualCapacity~=0)=AMat(residualCapacity~=0);

minCycleMean=findMinMaxCycleMean(costMat,true);
if(minCycleMean>=0)
    exitCode=0;
    return;
end

%If we get here, then the algorithm terminated without cancelling all
%negative cycles. Thus, the algorithm did not converge.
exitCode=2;
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
