function [minCostTuples,gain]=assign3DBB(C,maximize,boundType,initMethod,maxIter,epsVal)
%%ASSIGN3DBB Solve the data fusion axial 3D assignment problem using a
%         branch-and-bound algorithm. Such problems are NP-hard and thus
%         cannot be solved in polynomial time, so this function can only
%         solve problems of a moderate size. The optimization problem being
%         solved is minimize (or maximize)
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\sum_{k=1}^{n3}C_{i,j,k}*\rho_{i,j,k}
%         subject to
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\rho_{i,j,k}<=1 for all k
%         \sum_{i=1}^{n1}\sum_{k=1}^{n3}\rho_{i,j,k}<=1 for all j
%         \sum_{j=1}^{n2}\sum_{k=1}^{n3}\rho_{i,j,k} =1 for all i
%         \rho_{i,j,k} = 0 or 1
%         assuming that n1<=n2<=n3, and C is and n1Xn2Xn3 cost matrix.
%         This is equivalent to the optimization problem
%         min (or max) sum_{i} C(i,phi_2(i),phi_3(i))
%         where phi_2 and phi_3 are length n1 arrangements of n2 and n3
%         items over which the minimization is performed. In general, it is
%         not required that n1<=n2<=n3 (for any other ordering, the
%         equality constraint is on the dimension with the fewest
%         elements), so the solution is given as tuples solving the problem
%         min (or max) sum_{i} C(tuple(1,i),tuple(2,i),tuple(3,i))
%         where the tuples satisfy all of the above constraints.
%
%INPUTS: C An n1Xn2Xn3 cost hypermatrix. C cannot contain any NaNs and the
%          largest finite element minus the smallest element is a finite
%          quantity (does not overflow) when performing minimization and
%          where the smallest finite element minus the largest element is
%          finite when performing maximization.  Forbidden assignments can
%          be given costs of +Inf for minimization and -Inf for
%          maximization. During minimization, there should be no elements
%          with -Inf cost and no elements with +Inf cost during
%          maximization.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
% boundType This selects the bound to use for the branch and bound method.
%          These correspond to the method input of the assign3DLB function.
%          The default if this parameter is omitted or an empty matrix is
%          passed is 2.
% initMethod A parameter indicating how the branch-and-bound algorithm
%          should be initialized. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Do not
%            use an initial estimate.
%          1 Use an initial solution (and lower bound) via algorithm 0 of
%            the assign3D function.
%  maxIter If initMethod=1, then this is the maximum number of iterations
%          to allow the initialization routine in assign3D to perform. If
%          initMethod~=1, then this parameter is ignored. The default if
%          omitted or an empty matrix is passed is 10.
%   epsVal If initMethod=1, then this is both the AbsTold and RelTol inputs
%          to assign3D to determine whether the function has converged in
%          terms of the relative duality gap. If initMethod~=1, then this
%          parameter is ignored. The default if omitted or an empty matrix
%          is passed is eps(1).
%
%OUTPUTS: tuples A 3Xn1 matrix where tuples(:,i) is the ith assigned tuple
%                in C as described above.
%           gain The cost value of the optimal assignment found. This is
%                the sum of the values of the assigned elements in C.
%
%The algorithm is similar to the branch-and bound procedure describes in
%[1], which uses the "classical" branching of [2]. There is a total of
%min([n1,n2,n3]) elements to choose out of the entire C matrix. Suppose
%that n1=min([n1,n2,n3]). Then, the tuples chosen will have the first index
%from 1 to n1 and the indices of the remaining two tuple components will
%vary. The branching here is on the tuple chosen for i=1 of C(i,j,k), then
%the next level is for i=2 until branching is done for i=n1. At each level
%of branching there are n2*n3 possibilities, as is illustrated in Fig. 1 of
%[1]. Each time a branch is taken, all tuples of (j,k) sharing any common
%tuples with the branch taken (and any previous branches taken) must be
%removed. At each level, lower bounds for all n2*n3 possible branches are
%computed. Those branches with lower bounds that are less than the best
%found solution are discarded.
%
%REFERENCES:
%[1] W. P. Pierskalla, "The multidimensional assignment problem,"
%    Operations Research, vol. 16, no. 2, pp. 422-431, Mar-Apr. 1968.
%[2] R. E. Burkard and R. Rudolf, "Computational investigations on 3-
%    dimensional axial assignment problems," Belgian Journal of Operations
%    Research, Statistics and Computer Science, vol. 32, no. 1-3, pp.
%    85-98, 1993.
%
%September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(epsVal))
	epsVal=eps(1); 
end

if(nargin<5||isempty(maxIter))
    maxIter=10;
end

if(nargin<4||isempty(initMethod))
	initMethod=0; 
end

if(nargin<3||isempty(boundType))
    boundType=2; 
end

if(nargin<2||isempty(maximize))
	maximize=false; 
end

%For some types of lower bounds, the cost matrix must have all non-negative
%elements for the assignment algorithm to work. This forces all of the
%elements to be positive. The delta is added back in when computing the
%gain in the end.
if(maximize==true)
    CDelta=max(C(:));
    C=-C+CDelta;
else
    CDelta=min(C(:));
    C=C-CDelta;
end

%The algorthm was implemented assuming that n1<=n2<=n3 for the
%dimensionality of C. If a matrix having a different arrangement of indices
%is passed, one could permute the indices to make the assumptions below
%hold.
nVals=size(C);
%Deal with trailing singleton dimensions.
if(length(nVals)==2)
    nVals=[nVals,1];
elseif(length(nVals)==1)
    nVals=[nVals,1,1];
end

if(length(nVals)<3||~(nVals(1)<=nVals(2)&&nVals(2)<=nVals(3)))
    error('It is required that size(C,1)<=size(C,2)<=size(C,3)')
end

n1=nVals(1);
n2=nVals(2);
n3=nVals(3);

%Get an initial estimate.
switch(initMethod)
    case 0%Do not use an initial estimate.
        q=-Inf;
        gain=Inf;
        minCostTuples=[];
    case 1
        subgradMethod=0;
        maximize=false;
        [minCostTuples,gain,q,exitCode]=assign3D(C,false,subgradMethod,maximize,maxIter,epsVal,epsVal);
        if(exitCode==-3||exitCode==-1)
            %If no valid assignment was found, then we just start the
            %algorithm without initialization. This most commonly occurs
            %when the problem is infeasible and no solution will be found
            %in the end anyway.
            q=-Inf;
            gain=Inf;
            minCostTuples=[];
        end
    otherwise
        error('Invalid initMethod specified.');
end

%If the initial approximation returned a suboptimal solution.
if(abs(gain-q)>=epsVal*gain)
    numLevels=n1;%The number of things to assign.

    %There are n1*n2*n3 total tuples. However, we define a level by the
    %value of the first index of the tuple. In the initial level, there are
    %n2*n3 possible tuples (fix the value of the first index of C and the
    %other two are free). In subsequent levels, there are
    %(n2-curLevel+1)*(n3-curLevel+1) possibilities.
    maxTuplesPerLevel=n2*n3;

    %Generate all possible tuples of the second and third coordinates that
    %can be assigned. The assignment of the first coordinate is the level
    %of the recursion below.
    origFreeTuples=genAllTuples([n2-1;n3-1],false)+1;
    
    %The list of tuples that are candidates in the current level.
    freeTuples=zeros(2,maxTuplesPerLevel,numLevels);
    %Save the tuples that could possibly be assigned in the first level.
    freeTuples(:,:,1)=origFreeTuples;
    
    %These values link the modified tuples at each level to the original
    %tuples. These are important, because the ordering of the tuples in the
    %lower levels will be changed.
    freeTupleOrigIdx=zeros(maxTuplesPerLevel,numLevels);
    freeTupleOrigIdx(:,1)=1:maxTuplesPerLevel;

    %This keeps track of how many tuples are being considered at each
    %level after eliminating those that conflict with prior assignments and
    %eliminating those whose lower bounds are too high.
    numFreeTuples=zeros(numLevels,1);
    %In the top level, all possibilities are initially considered.
    numFreeTuples(1)=maxTuplesPerLevel;
    
    %The list of tuples that were candidates when entering the current
    %level. In the current level, some tuples might be removed due to a low
    %cost at that level. However, when going to a higher level, all tuples
    %that do not conflict with assignments at lower levels must be
    %considered. Thus, this stores the value of the tuples in the current
    %level before pruning.
    freeTuplesEntering=zeros(2,maxTuplesPerLevel,numLevels);
    freeTuplesEntering(:,:,1)=origFreeTuples;

    %Like freeTupleOrigIdx, these link the values in freeTuplesEntering to
    %the original set of tuples.
    freeTuplesEnteringOrigIdx=zeros(maxTuplesPerLevel,numLevels);
    freeTuplesEnteringOrigIdx(:,1)=1:maxTuplesPerLevel;

    %This keeps track of how many tuples are being considered at each
    %level after eliminating those that conflict with prior assignments but
    %not eliminating anything where the lower bounds are too high.
    numFreeTuplesEntering=zeros(numLevels,1);
    %In the top level, all possibilities are considered.
    numFreeTuplesEntering(1)=maxTuplesPerLevel;

    %The tuple that is assigned at each level, after the indexation has
    %been changed to relate it to the current cost submatrix at that level.
    %We don't need to save the first index here, because it is just the
    %level number. We don't need to save the assigned tuple in the maximum
    %level, because we use it right away.
    assignedTuples=zeros(2,numLevels-1);
    
    %This links the assigned tuple to the index of the tuple in
    %origFreeTuples. It is simpler to use this than to try to reconstruct
    %the tuple from assignedTuples, where each tuple has been modified
    %based on assignments at lower levels. Unlike assignedTuples, this
    %needs to be numLevels in size.
    assignedTupleOrigIdx=zeros(numLevels,1);

    %The cumulative assigned cost at each level. The final cost is computed
    %at the top level and does not need to be stored in this.
    cumAssignedCost=zeros(numLevels-1,1);

    %Allocate space for the bounded total costs for assigning something at
    %each level. Each level is given by the value in the first index that
    %is assigned. For each first index, there are n2*n3 values. There are
    %no bounds for the deepest level, because the cost matrix can be
    %directly used.
    boundVals=zeros(maxTuplesPerLevel,numLevels-1);

    %This holds what was the minimum cost value before going to a higher
    %level of recursion. Thus, when returning to a particular level, if the
    %minimum cost value changed, one should check the previous bounds and
    %see if any of the candidate tuples should be thrown out and never
    %visited.
    prevMinCost=zeros(numLevels-1,1);

    %Allocate space for the marginal costs at each level (the cost
    %submatrices). We could just allocate this like
    %CostMatCur=zeros(n1,n2,n3,numLevels); and then easily address each
    %level with the final index. However, that wastes a lot of space,
    %because each level removes one element from every dimension. The first
    %matrix is n1Xn2Xn3 in size, the next is (n1-1)*(n2-1)*(n3-1) in size
    %and so on until the final level which is 1X(n2-n1+1)X(n3-n1+1) in
    %size. To save space, we shall allocate the minimum amount of space
    %required, as well as an offset array, so that we know where the
    %elements of the matrices start and end so that we can later consider
    %each matrix individually.
    costMatCurStartIdx=zeros(numLevels+1,1);
    costMatCurStartIdx(1)=1;
    for curLevel=2:(numLevels+1)
        offset=curLevel-2;
        costMatCurStartIdx(curLevel)=costMatCurStartIdx(curLevel-1)+(n1-offset)*(n2-offset)*(n3-offset);
    end
    %The final value in costMatCurStartIdx is one past the end of the
    %matrix.
    CostMatCur=zeros(costMatCurStartIdx(numLevels+1)-1,1);
    %Assign the marginal cost matrix of the first level.
    CostMatCur(costMatCurStartIdx(1):(costMatCurStartIdx(2)-1))=C(:);

    %Enter the recursion (unrolled).
    curLevel=1;
    increasingLevelNumber=true;
    while(curLevel>0)
        if(curLevel==numLevels)
            %If we have reached the top level, the final assignment is just
            %choosing the most likely free tuple in the current cost
            %matrix.

            minVal=Inf;
            minIdx=1;
            for curTuple=1:numFreeTuples(curLevel)
                i=1;
                j=freeTuples(1,curTuple,curLevel);
                k=freeTuples(2,curTuple,curLevel);
                %The dimensionality of the cost matrix at the current
                %level.
                dims=[n1-curLevel+1;n2-curLevel+1;n3-curLevel+1];
                %The index within the cost matrix at the current level.
                idx=sub2ind(dims,i,j,k);
                curCost=CostMatCur(costMatCurStartIdx(curLevel)+idx-1);

                if(curCost<minVal)
                    minVal=curCost;
                    minIdx=curTuple;
                end
            end
            if(numLevels>1)
                minVal=cumAssignedCost(curLevel-1)+minVal;
            end

            if(minVal<gain)%If a new global minimum has been found.
                gain=minVal;
                
                assignedTupleOrigIdx(curLevel)=freeTupleOrigIdx(minIdx,curLevel);
                
                minCostTuples=zeros(3,numLevels);
                for i=1:numLevels
                    minCostTuples(:,i)=[i;origFreeTuples(:,assignedTupleOrigIdx(i))];
                end
            end

            increasingLevelNumber=false;
            curLevel=curLevel-1;
            continue;
        elseif(increasingLevelNumber==true)
            %We have just entered this level, so we have to first compute
            %the lower bounds for the current level.

            curTuple=1;
            while(curTuple<=numFreeTuples(curLevel))
                i=1;
                j=freeTuples(1,curTuple,curLevel);
                k=freeTuples(2,curTuple,curLevel);
                
                %The i,j,k indices of the tuple are with respect to indices
                %in the shrunken cost matrix CostMatCur(i,j,k,curLevel),
                %not to the set of tuples with respect to the global cost
                %matrix C. Since we are picking off rows one at a time, i
                %will always be 1.
                
                %The dimensionality of the cost matrix at the current
                %level.
                dims=[n1-curLevel+1;n2-curLevel+1;n3-curLevel+1];
                %The index within the cost matrix at the current level.
                idx=sub2ind(dims,i,j,k);
                curCost=CostMatCur(costMatCurStartIdx(curLevel)+idx-1);

                %If the assignment is not allowed.
                if(~isfinite(curCost))
                    curTuple=curTuple+1;
                    continue;
                end

                %Create the cost matrix at the next level up by removing
                %the appropriate entry from each dimension of the current
                %cost matrix. First, we extract the matrix for the current
                %level.
                CCur=reshape(CostMatCur(costMatCurStartIdx(curLevel):(costMatCurStartIdx(curLevel+1)-1)),[n1-curLevel+1,n2-curLevel+1,n3-curLevel+1]);
                %The indices of the cost matrix one level up.
                spanNext=costMatCurStartIdx(curLevel+1):(costMatCurStartIdx(curLevel+2)-1);
                CostMatCur(spanNext)=CCur(2:(n1-curLevel+1),[1:(j-1),(j+1):(n2-curLevel+1)],[1:(k-1),(k+1):(n3-curLevel+1)]);
                
                if(curLevel>1)
                    boundVals(curTuple,curLevel)=cumAssignedCost(curLevel-1)+curCost+assign3DLB(reshape(CostMatCur(spanNext),[n1-curLevel,n2-curLevel,n3-curLevel]),boundType);
                else
                    boundVals(curTuple,1)=curCost+assign3DLB(reshape(CostMatCur(spanNext),[n1-curLevel,n2-curLevel,n3-curLevel]),boundType); 
                end

                %If the bound is no better than the current best solution,
                %then remove that tuple as a contender in this level. We
                %just overwrite it with the last tuple, decriment the
                %number of tuples present and then revisit the tuple with
                %the same index (which will now be the last tuple).
                if(boundVals(curTuple,curLevel)>=gain)
                    CostMatCur(costMatCurStartIdx(curLevel)+idx-1)=Inf;
                    
                    freeTuples(:,curTuple,curLevel)=freeTuples(:,numFreeTuples(curLevel),curLevel);
                    freeTupleOrigIdx(curTuple,curLevel)=freeTupleOrigIdx(numFreeTuples(curLevel),curLevel);
                    numFreeTuples(curLevel)=numFreeTuples(curLevel)-1;
                else
                    curTuple=curTuple+1;
                end
            end

            %If there are not enough tuples left to make a full assignment.
            %The assumption here is that numLevels>1, so if we are here,
            %there must be one tuple left to assign after assigning the one
            %in this level.
            if(numFreeTuples(curLevel)==0)
                increasingLevelNumber=false;
                curLevel=curLevel-1;
                continue;
            end

            %Sort the tuples by cost.
            [boundVals(1:numFreeTuples(curLevel),curLevel),idx]=sort(boundVals(1:numFreeTuples(curLevel),curLevel),'descend');
            freeTuples(:,1:numFreeTuples(curLevel),curLevel)=freeTuples(:,idx,curLevel);
            freeTupleOrigIdx(1:numFreeTuples(curLevel),curLevel)=freeTupleOrigIdx(idx,curLevel);    

            %Record the current minimum cost value.
            prevMinCost(curLevel)=gain;
        end

        %If we are here, then either increasingLevelNumber=false and we are
        %backtracking, or increasingLevelNumber=true, and we have just
        %entered this level.
        if(increasingLevelNumber==false)
            %If we are backtracking, then we remove the last tuple visited
            %in this level from consideration.
            i=1;
            j=assignedTuples(1,curLevel);
            k=assignedTuples(2,curLevel);
            
            %The dimensionality of the cost matrix at the current
            %level.
            dims=[n1-curLevel+1;n2-curLevel+1;n3-curLevel+1];
            %The index within the cost matrix at the current level.
            idx=sub2ind(dims,i,j,k);
            %freeTuples(:,curLevel,numFreeTuples(curLevel),curLevel)
            %contains the last assigned tuple for this level. Here, we
            %eliminate it from further consideration by setting
            %CostMatCur(i,j,k,curLevel) to Inf and reducing the count of
            %the tuple by one.
            CostMatCur(costMatCurStartIdx(curLevel)+idx-1)=Inf;
            numFreeTuples(curLevel)=numFreeTuples(curLevel)-1;
            
            %If the last tuple visited leaves too few tuples left for a
            %full assignment, then go up another level.
            if(numFreeTuples(curLevel)==0)
                increasingLevelNumber=false;
                curLevel=curLevel-1;
                continue;
            end

            %Also, if the minimum bound changed, remove all tuples whose
            %bounds are larger than the new minimum.
            if(prevMinCost(curLevel)~=gain)
                prevMinCost(curLevel)=gain;

                %Find all of the tuples where the bound is no better than
                %the current best solution and remove them as contenders in
                %this level.

                %Remember that when entering this level, the largest bounds
                %are at the beginning of the boundVals array, because it
                %has been sorted in decreasing order. Thus, we have to find
                %the first index where the bound is <gain.

                idx=find(boundVals(1:numFreeTuples(curLevel),curLevel)<gain,1);

                %If this gets rid of all tuples.
                if(isempty(idx))
                    curLevel=curLevel-1;
                    continue;
                elseif(idx>1)
                    %If here, remove the tuples with bounds that are too
                    %large.
                    sel1=1:(numFreeTuples(curLevel)-idx+1);
                    sel2=idx:numFreeTuples(curLevel);
                    numFreeTuples(curLevel)=numFreeTuples(curLevel)-idx+1;

                    boundVals(sel1,curLevel)=boundVals(sel2,curLevel);

                    %Set the entries in the cost matrix associated with the
                    %tuples that are being removed to Inf, so they won't be
                    %assigned.
                    for curTuple=1:(idx-1)
                        i=1;
                        j=freeTuples(1,curTuple,curLevel);
                        k=freeTuples(2,curTuple,curLevel);
                        %The dimensionality of the cost matrix at the
                        %current level.
                        dims=[n1-curLevel+1;n2-curLevel+1;n3-curLevel+1];
                        %The index within the cost matrix at the current
                        %level.
                        idx=sub2ind(dims,i,j,k);
                        CostMatCur(costMatCurStartIdx(curLevel)+idx-1)=Inf;
                    end

                    freeTuples(:,sel1,curLevel)=freeTuples(:,sel2,curLevel);
                    freeTupleOrigIdx(sel1,curLevel)=freeTupleOrigIdx(sel2,curLevel);
                end
            end

            increasingLevelNumber=true;
        end

        %The tuple having the lowest cost is visited first.
        %Assign the tuple in this current level.
        assignedTuples(:,curLevel)=freeTuples(:,numFreeTuples(curLevel),curLevel);
        assignedTupleOrigIdx(curLevel)=freeTupleOrigIdx(numFreeTuples(curLevel),curLevel);

        %Add the cost of the assigned tuple to the cumulative assigned
        %cost.
        i=1;
        j=assignedTuples(1,curLevel);
        k=assignedTuples(2,curLevel);
        %The dimensionality of the cost matrix at the current
        %level.
        dims=[n1-curLevel+1;n2-curLevel+1;n3-curLevel+1];
        %The index within the cost matrix at the current level.
        idx=sub2ind(dims,i,j,k);
        curCost=CostMatCur(costMatCurStartIdx(curLevel)+idx-1);

        if(curLevel>1)
            cumAssignedCost(curLevel)=cumAssignedCost(curLevel-1)+curCost;
        else
            cumAssignedCost(curLevel)=curCost;
        end

        %Construct the cost matrix for the next level. This means removing
        %the entire row, column, etc. in each dimension that contains the
        %assigned tuple.
        CCur=reshape(CostMatCur(costMatCurStartIdx(curLevel):(costMatCurStartIdx(curLevel+1)-1)),[n1-curLevel+1,n2-curLevel+1,n3-curLevel+1]);
        %The indices of the cost matrix one level up.
        spanNext=costMatCurStartIdx(curLevel+1):(costMatCurStartIdx(curLevel+2)-1);
        CostMatCur(spanNext)=CCur(2:(n1-curLevel+1),[1:(j-1),(j+1):(n2-curLevel+1)],[1:(k-1),(k+1):(n3-curLevel+1)]);

        %Copy the tuples that can be considered for the next and subsequent
        %levels into freeTuples(:,:,curLevel+1) and
        %freeTuplesEntering(:,:,curLevel+1).
        %This draws tuples for the next level down from
        %freeTuplesEntering(:,:,curLevel), because freeTuples(:,:,curLevel)
        %may have had tuples removed only for curLevel due to bounding. At
        %the same time as the tuples are added, the coordinates of all
        %valid tuples are adjusted to deal with the removed elements of
        %CostMatCur.
        
        numFreeTuples(curLevel+1)=0;
        for curTuple=1:numFreeTuplesEntering(curLevel)
            if(any(freeTuplesEntering(:,curTuple,curLevel)==assignedTuples(:,curLevel)))
                continue;
            else
                numFreeTuples(curLevel+1)=numFreeTuples(curLevel+1)+1;

                %Free tuples are indexed with respect to the shrunken
                %matrix of costs in CostMatCur(:,:,:,curLevel+1). This
                %means that for index values > the assigned index, one must
                %be subtracted.
                freeTuples(:,numFreeTuples(curLevel+1),curLevel+1)=freeTuplesEntering(:,curTuple,curLevel)-(freeTuplesEntering(:,curTuple,curLevel)>assignedTuples(:,curLevel));
                freeTupleOrigIdx(numFreeTuples(curLevel+1),curLevel+1)=freeTuplesEnteringOrigIdx(curTuple,curLevel);
            end
        end
        numFreeTuplesEntering(curLevel+1)=numFreeTuples(curLevel+1);
        freeTuplesEntering(:,1:numFreeTuples(curLevel+1),curLevel+1)=freeTuples(:,1:numFreeTuples(curLevel+1),curLevel+1);
        freeTuplesEnteringOrigIdx(1:numFreeTuples(curLevel+1),curLevel+1)=freeTupleOrigIdx(1:numFreeTuples(curLevel+1),curLevel+1);

        %Go to the next level.
        curLevel=curLevel+1;
    end
end

%Adjust the gain for the initial offset of the cost matrix.
if(maximize==true)
    gain=-gain+CDelta*n1;
else
    gain=gain+CDelta*n1;
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
