function [curSetPart,recurSetData]=getNextSetMPartition(param1,m,startIdx)
%%GETNEXTSETPARTITION Given a set of n labeled objects, generate all
%                ways of partitioning the objects into m subsets, given
%                that the order of the objects in the subsets does not
%                matter. Sets are numbered starting at 0. Whereas
%                genAllSetPartitions go through all partition sizes, this
%                function is limited only to partition sizes of m.
%
%INPUTS: param1 If one wishes to get the first length-m set partition in
%               the sequence, this this should be n, the number of objects
%               in the total set and the second input, m, must be given. If
%               one wishes to get the next set in the sequence, then param1
%               is the recurData returned by the function and no other
%               inputs can be given.
%             m If this is given, then the first m-part set partition in
%               the sequence is produced. This is the number of parts. This
%               should not be provided when subsequent set partitions in
%               the sequence are desired.
%      startIdx The starting index of the set partitions. If omitted or an
%               empty matrix is passed, the default of 1 is used. This
%               should only be passed when the first combination is
%               desired. It should not be passed when recurData is passed
%               as param1 to get subsequent length-m set partitions.
%
%OUTPUTS: curSetPart The current nX1 length-m set partition with the first
%                    set given index startIdx.
%          recurData A data structure that can be passed back to this
%                    function to get subsequent binary combinations.
%
%There is a total of StirlingNumber2(n,m) length-m set partitions.
%
%This is a callback-form of Algorithm 1 in genAllSetMPartitions. The
%algorithm works by using getNextMPartition to go through all possible
%sizes of the set partitions and then going from the largest to the
%smallest sized partition group, use getNextCapacitatedCombo to get the
%actual combinations, noting that once a combination of elements has been
%assigned to one partition group, those indices are not available to
%subsequent ones.
%
%In Matlab, additional speed savings could be obtained by implementing this
%in a handle class, because we cannot directly overwrite the elements of
%the recurData structure. Rather, being passed by value, copies are made
%and new copies have to be formed for the next return value.
%
%This would be more efficient if implement as a class, because we could
%make an iterator and use pointers rather than copying the data between
%calls.
%
%This function implements the algorithm of Ehrlich in [1], which is
%explained in [2]. For the m=2 case, the function genAllBinCombinations is
%used and a special condition for m=1 has been added. A modification to
%the algorithm to avoid negative indexation has been made.
%
%EXAMPLE:
%Here, we show that this produces the same results as the equivalent
%algorithm in genAllSetMPartitions.
% n=6;
% m=4;
% startIdx=1;
% allMSets=genAllSetMPartitions(n,m,startIdx,1);
% 
% numMSetPars=StirlingNumber2(n,m);
% allMRecurSets=zeros(n,numMSetPars);
% [curSetPart,recurData]=getNextSetMPartition(n,m,startIdx);
% allMRecurSets(:,1)=curSetPart;
% curIdx=2;
% 
% [curSetPart,recurData]=getNextSetMPartition(recurData);
% while(~isempty(curSetPart))
%     allMRecurSets(:,curIdx)=curSetPart;
%     curIdx=curIdx+1;
%     [curSetPart,recurData]=getNextSetMPartition(recurData);
% end
% assert(all(allMRecurSets(:)==allMSets(:)))
%The assertion throws no error, because the results are the same.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If the first set partition is desired.
if(nargin>1)
    if(nargin<3||isempty(startIdx))
        startIdx=1; 
    end
    
    n=param1;
    
    idx2Assign=zeros(n,m);
    numIdx2Assign=zeros(m,1);
    %All indices can be assigned at the start.

    curSetPartition=zeros(n,1);

    %To store the data used by getNextCapacitatedCombo.
    allRecurData=cell(m,1);

    %All indices are available for the first partition set. This does not
    %change.
    idx2Assign(:,1)=1:n;
    numIdx2Assign(1)=n;

    curPartition=getNextMPartition(n,m);
        
    %Set the recursion data that does not change.
    recurSetData.n=n;
    recurSetData.m=m;
    recurSetData.startIdx=startIdx;
    
    skipOuterLoop=false;
else
    %Get the recursion set data that does not change.
    recurSetData=param1;
    n=recurSetData.n;
    m=recurSetData.m;
    startIdx=recurSetData.startIdx;
    
    %Get the recursion set data that changes.
    idx2Assign=recurSetData.idx2Assign;
    numIdx2Assign=recurSetData.numIdx2Assign;
    curSetPartition=recurSetData.curSetPartition;
    allRecurData=recurSetData.allRecurData;
    curPartition=recurSetData.curPartition;
    
    %Get the variables in the loop.
    partSizes=recurSetData.partSizes;
    numReps=recurSetData.numReps;
    numUniqueParts=recurSetData.numUniqueParts;
    startPartNum=recurSetData.startPartNum;
    backtracking=recurSetData.backtracking;
    curPartSet=recurSetData.curPartSet;
    
    skipOuterLoop=true;
end

while(1)
    if(skipOuterLoop==false)
        %If we have gone past the final combination.
        if(isempty(curPartition))
            curSetPart=[];
            return
        end
        
        %curPartition contains values in decreasing order. Starting at the
        %largest partition and going to the smallest, we form all combinations
        %of assigning indices to the partitions. Partitions are numbered from
        %1.

        %First, we get the unique partition sizes and the number of repeated
        %partition sizes. 
        [partSizes,numReps]=runLenEncode(curPartition);
        numUniqueParts=length(partSizes);

        startPartNum=cumsum([0;numReps(1:(numUniqueParts-1))]);

        %Recursively go through all of the capacitated combinations of
        %assigning the n elements to the repeated partitions. 
        curPartSet=1;
        backtracking=false;
        skipOuterLoop=true;
    end
    while(curPartSet>0)
        if(backtracking==false)%If we just entered this level.
            [theCombo,allRecurData{curPartSet}]=getNextCapacitatedCombo(partSizes(curPartSet),numReps(curPartSet),numIdx2Assign(curPartSet));

            assignCapacitatedCombo(theCombo);

            if(curPartSet==numUniqueParts)
                %Save the current set partition.
                curSetPart=curSetPartition+(startIdx-1);

                %Repeat the current level but get into backtracking mode,
                %because after finishing this level, it will be necessary
                %to backtrack.
                backtracking=true;
                
                saveRecurData();
                return
            end
            %Go to the next level
            curPartSet=curPartSet+1;
            continue;
        else
            %If backtracking, we need to go to the next capacitated combo
            %in this level.
            [theCombo,allRecurData{curPartSet}]=getNextCapacitatedCombo(allRecurData{curPartSet});

            if(isempty(theCombo))
                %If we have passed the last capacitated combo at this
                %level.
                curPartSet=curPartSet-1;
                continue;%Keep backtracking.
            end

            assignCapacitatedCombo(theCombo);

            if(curPartSet==numUniqueParts)
                %Save the current set partition.
                curSetPart=curSetPartition+(startIdx-1);
                
                saveRecurData();
                return;
            end

            %Go to the next level
            curPartSet=curPartSet+1;
            backtracking=false;
            continue;
        end
    end

    skipOuterLoop=false;
    curPartition=getNextMPartition(curPartition);  
end

    function saveRecurData()
        %Assign the variables that change:
        recurSetData.idx2Assign=idx2Assign;
        recurSetData.numIdx2Assign=numIdx2Assign;
        recurSetData.curSetPartition=curSetPartition;
        recurSetData.allRecurData=allRecurData;
        recurSetData.curPartition=curPartition;
    
        %Get the variables in the loop.
        recurSetData.partSizes=partSizes;
        recurSetData.numReps=numReps;
        recurSetData.numUniqueParts=numUniqueParts;
        recurSetData.startPartNum=startPartNum;
        recurSetData.backtracking=backtracking;
        recurSetData.curPartSet=curPartSet;
    end

    function assignCapacitatedCombo(theCombo)
        totalBins=size(theCombo,1);
        totalReps=size(theCombo,2);

        %Initially, for the next level, all indices are avaiable to
        %assign.
        if(curPartSet<numUniqueParts)
            %Indices for the next level are removed as they are
            %assigned.
            idx2Assign(1:numIdx2Assign(curPartSet),curPartSet+1)=idx2Assign(1:numIdx2Assign(curPartSet),curPartSet);
            numIdx2Assign(curPartSet+1)=numIdx2Assign(curPartSet);
        end
        
        %Assign the values in each of the numReps partitions starting from
        %startPartNum(curPartSet) for the partition numbering.
        partNumOffset=startPartNum(curPartSet);

        for curRep=1:totalReps
            for curBin=1:totalBins
                partNum=partNumOffset+curRep;
                curBaseIdx=theCombo(curBin,curRep);

                curIdx=idx2Assign(curBaseIdx,curPartSet);
                %Assign the index.
                curSetPartition(curIdx)=partNum;
            end
        end
        
        %Remove the assigned indices from the next level.
        if(curPartSet<numUniqueParts)
            %As opposed to explicitely deleting things, we just overwrite
            %the element with the last element in the array and then reduce
            %the indication of the array size.
            idx2Del=sort(theCombo(:),'descend');
            for k=1:(totalReps*totalBins)
                curBaseIdx=idx2Del(k);

                idx2Assign(curBaseIdx,curPartSet+1)=idx2Assign(numIdx2Assign(curPartSet+1),curPartSet+1);
                numIdx2Assign(curPartSet+1)=numIdx2Assign(curPartSet+1)-1;
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
