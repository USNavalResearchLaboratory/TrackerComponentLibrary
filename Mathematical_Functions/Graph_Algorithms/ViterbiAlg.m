function [minCostPath,minCost]=ViterbiAlg(costMats)
%%VITERBIALG Use the Viterbi algorithm to find the recursive, optimal
%            (minimum cost) solution to a trellis optimization problem.
%            This can be, for example to estimate the state sequence of a
%            discrete-time finite-state Markov process. The Viterbi
%            algorithm is a form of dynamic programming. The algorithm is
%            implemented for a fixed number of states at each step. If the
%            number of states varies, then just put Inf for the cost of
%            the unused states at each step.
%
%INPUTS: costMats A numStatesXnumStatesXnumTransitions hypermatrix of costs
%                 of transitioning between states at each state. This is a
%                 set of adjacency matrices, one for each stage. There are
%                 numStates states at each stage. To go from the kth state
%                 to the k+1th state, the const costMats(i,j,k) is the cost
%                 of transitioning from state i forward to state j.
%                 Infinite costs can be used when transitions are not
%                 allowed. The cumulative cost function is the sum of costs
%                 along the path through the states.
%
%OUTPUTS: minCostPath The indices of the states in the minimum cost path.
%                     This is a numTransitions+1 vector. When multiple
%                     paths have equal cost, this just chooses the first
%                     one.
%             minCost The cost of the minCost path. This is the sum of the
%                     transition costs in costMats.
%
%A good explanation of the Viterbi algorithm is given in [1]. The total
%cost is the sum of all of the edge costs. The example in Figure 8 of the
%paper can be implemented using
% costMats=inf(4,4,5);
% costMats(1,1,1)=1;
% costMats(1,3,1)=1;
% costMats(1,1,2)=1;
% costMats(1,3,2)=1;
% costMats(3,2,2)=2;
% costMats(3,4,2)=0;
% costMats(1,1,3)=0;
% costMats(1,3,3)=2;
% costMats(2,1,3)=2;
% costMats(2,3,3)=0;
% costMats(3,2,3)=1;
% costMats(3,4,3)=1;
% costMats(4,2,3)=1;
% costMats(4,4,3)=1;
% costMats(1,1,4)=2;
% costMats(2,1,4)=0;
% costMats(3,2,4)=1;
% costMats(4,2,4)=1;
% costMats(1,1,5)=1;
% costMats(2,1,5)=1;
% [minCostPath,minCost]=ViterbiAlg(costMats)
%whereby the minimum cost path is [1;3;4;2;1;1] and the minimum cost is 3.
%
%REFERENCES:
%[1] G. D. Forney Jr., "The Viterbi algorithm," Proceedings of the IEEE,
%    vol. 61, no. 3, pp. 268-278, Mar. 1973.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Renormalize each step to keep the values reasonable.
%%If the costMats is a single matrix.
numStates=size(costMats,1);
numTransitions=size(costMats,3);
%This is necessary to decode the least-cost path through the trellis at the
%end of the algorithm. This holds the idex of the minimum node back at each
%step
backwardPath=zeros(numStates,numTransitions);
curCumCosts=zeros(numStates,1);
prevCumCosts=zeros(numStates,1);
cumCostOffset=0;

for curTrans=1:numTransitions
    %Each destination state only takes the minimum cost input from the past
    %states.
    overallMinCost=Inf;
    for curDestState=1:numStates
        [minCost,minIdx]=min(prevCumCosts+costMats(:,curDestState,curTrans));
        curCumCosts(curDestState)=minCost;
        backwardPath(curDestState,curTrans)=minIdx;
        
        if(minCost<overallMinCost)
            overallMinCost=minCost;
        end
    end
    
    %To minimize finite precision problems if the cumulative path costs
    %have large magnitudes compared to the individual path costs, we will
    %subtract off the smallest value from all other values. The cumulative
    %effects of these subtractions needs to be taken into account to get
    %the correct value for the minimim cost path the return.
    curCumCosts=curCumCosts-overallMinCost;
    cumCostOffset=cumCostOffset+overallMinCost;
    prevCumCosts=curCumCosts;
end

%We now need to extract the minimum cost path.
minCostPath=zeros(numTransitions+1,1);
[minCost,minIdx]=min(curCumCosts);
minCostPath(end)=minIdx;
minCost=minCost+cumCostOffset;

for curTrans=numTransitions:-1:1
    minCostPath(curTrans)=backwardPath(minCostPath(curTrans+1),curTrans);
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
