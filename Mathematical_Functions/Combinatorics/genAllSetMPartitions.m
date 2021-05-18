function allSets=genAllSetMPartitions(n,m,startIdx,algorithm)
%%GENALLSETMPARTITIONS Given a set of n labeled objects, generate all
%                ways of partitioning the objects into m subsets, given
%                that the order of the objects in the subsets does not
%                matter. Sets are numbered starting at 0. Whereas
%                genAllSetPartitions go through all partition sizes, this
%                function is limited only to partition sizes of m.
%
%INPUTS: n The number of objects to partition.
%        m The number of subsets to form; 1<=m<=n.
% startIdx The starting index of the set partitions. If omitted or an empty
%          matrix is passed, the default of 1 is used.
% algorithm An optional parameter indicating which algorithm to use to
%          generate the subset sequence. Possible values are:
%          0 Use the algorithm given by Knuth's solution to problem 17 of
%            Chapter 7.2.1.5 of [1], which is his realization of the
%            recursive algorithm of Ruskey in [2]. A special condition to
%            handle the m=1 case has been added.
%          1 (The default if omitted or an empty matrix is passed) Use
%            getNextMPartition to go through all possible sizes of the set
%            partitions and then going from the largest to the smallest
%            sized partition group, use getNextCapacitatedCombo to get the
%            actual combinations, noting that once a combination of
%            elements has been assigned to one partition group, those
%            indices are not available to subsequent ones.
%
%OUTPUTS: allSets The nXnumSetPars collection of all length-m set
%                 partitions. allSets(:,i) is the ith set partition and the
%                 n elements given the number of the set to which that
%                 element is assigned, starting with startIdx.
%
%There is a total of StirlingNumber2(n,m) length-m set partitions.
%
%Algorithm 0 uses recursion. The number of recursion can be quite high. For
%example, for n=10 and m=4, the algorithm goes down 6753 levels. However,
%algorithm 0 is typically faster than algorithm 1.
%
%Scanning from the first unique partition encounted to the last in each of
%the solutions from algorithm 0, the partition numbers increase. In
%algorithm 1, the same partition numbers are used, but the first partition
%encountered is not always the lowest indexed one. For example, in
%algorithms 0, allSets(1,i) for all i is startIdx. However, that is not
%always the case with algorithm 1.
%
%EXAMPLE:
% allSet0=genAllSetMPartitions(5,4,0,0)
% allSet1=genAllSetMPartitions(5,4,0,1)
%One will get:
%allSet0=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%         0, 1, 1, 1, 1, 1, 1, 1, 1, 1;
%         1, 1, 0, 2, 2, 2, 2, 2, 2, 2;
%         2, 2, 2, 2, 1, 0, 3, 3, 3, 3;
%         3, 3, 3, 3, 3, 3, 3, 2, 1, 0];
%allSet1=[0, 0, 0, 0, 1, 1, 1, 1, 1, 1;
%         0, 2, 2, 2, 0, 0, 0, 2, 2, 2;
%         3, 0, 3, 3, 0, 3, 3, 0, 0, 3;
%         1, 1, 0, 1, 2, 0, 2, 0, 3, 0;
%         2, 3, 1, 0, 3, 2, 0, 3, 0, 0];
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%[2] F. Ruskey, "Simple combinatorial Gray codes constructed by reversin
%    sublists," in International Symposium on Algorithms and Computation,
%    Hong Kong, pp.201-208, Dec. 1993.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(m>n)
    error('m must be <=n.')
end

if(nargin<4||isempty(algorithm))
    algorithm=0; 
end

if(nargin<3||isempty(startIdx))
    startIdx=1;
end

switch(algorithm)
    case 0
        allSets=genAllSetMPartitionsKnuth(n,m)+startIdx;
    case 1
        allSets=genAllSetMPartitionsCrouse(n,m)+(startIdx-1);
    otherwise
        error('Unknown algorithm specified.')
end
end

function allSets=genAllSetMPartitionsKnuth(n,m)
%%GENALLSETMPARTITIONSEHRLICH Generate all r-subset set paritions of n
%                   elements using Knuth's solution to problem 17 of
%                   Chapter 7.2.1.5 of [1], which is his realization of the
%                   recursive algorithm of Ruskey in [2]. A special
%                   condition to handle the m=1 case has been added. The
%                   subsets are numbered starting with 0.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%[2] F. Ruskey, "Simple combinatorial Gray codes constructed by reversin
%    sublists," in International Symposium on Algorithms and Computation,
%    Hong Kong, pp.201-208, Dec. 1993.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

totalNumSets=StirlingNumber2(n,m);

allSets=zeros(n,totalNumSets);

%The initial set partition.
a=[zeros((n-m),1);(0:(m-1)).'];
curPartIdx=1;

if(m==1)
    %The special case where everything is in one set.
    allSets=a;
    return;
end

f(m,n,0);
    function f(mu,nu,sigma)
    %%F The recursive "forward" function.
    %
    %July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
        if(mu==2)
            allSets(:,curPartIdx)=a;
            curPartIdx=curPartIdx+1;
        else
            f(mu-1,nu-1,mod(mu+sigma,2));
        end
        
        if(nu==mu+1)
            a(mu)=mu-1;

            allSets(:,curPartIdx)=a;
            curPartIdx=curPartIdx+1;

            while(a(nu)>0)
                a(nu)=a(nu)-1;
                allSets(:,curPartIdx)=a;
                curPartIdx=curPartIdx+1;
            end
        elseif(nu>mu+1)
            if(mod(mu+sigma,2))%Odd
                a(nu-1)=mu-1; 
            else%Even
                a(mu)=mu-1;
            end

            if(mod(a(nu)+sigma,2))%Odd
                b(mu,nu-1,0);
            else%Even
                f(mu,nu-1,0);
            end
            
            while(a(nu)>0)
                a(nu)=a(nu)-1;
                if(mod(a(nu)+sigma,2))%Odd
                    b(mu,nu-1,0);
                else%Even
                    f(mu,nu-1,0);
                end
            end
        end
    end

    function b(mu,nu,sigma)
    %%B The recursive "backward" function.
    %
    %July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        if(nu==mu+1)
            while(a(nu)~=mu-1)
                allSets(:,curPartIdx)=a;
                curPartIdx=curPartIdx+1;
                
                a(nu)=a(nu)+1;
            end

            allSets(:,curPartIdx)=a;
            curPartIdx=curPartIdx+1;
            
            a(mu)=0;
        elseif(nu>mu+1)
            if(mod(a(nu)+sigma,2))%Odd
                f(mu,nu-1,0);
            else%even
                b(mu,nu-1,0);
            end
            
            while(a(nu)<mu-1)
                a(nu)=a(nu)+1;
                if(mod(a(nu)+sigma,2))%Odd
                    f(mu,nu-1,0);
                else%even
                    b(mu,nu-1,0);
                end
            end
            
            if(mod(mu+sigma,2))%Odd
                a(nu-1)=0;
            else%Even
                a(mu)=0;
            end
        end
        
        if(mu==2)
            allSets(:,curPartIdx)=a;
            curPartIdx=curPartIdx+1;
        else
            b(mu-1,nu-1,mod(mu+sigma,2));
        end
    end
end

function setMPartList=genAllSetMPartitionsCrouse(n,m)
%%GENALLSETMPARTITIONS Generate all r-subset set paritions of n
%                   elements using an algorithm based on m-partitions and
%                   capacitated combinations. The subsets are numbered
%                   starting with 1.
%
%This algorithm uses getNextMPartition to go through all possible sizes of
%the set partitions. Going from the largest partition size to the smallest,
%we go through all combinations of assigning the n indices to the
%partition. However, since the partitions of equal size are unordered, we
%group them together and go through them using getNextCapacitatedCombo.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

totalSetMPartitions=StirlingNumber2(n,m);
setMPartList=zeros(n,totalSetMPartitions);
curSetMPartition=1;

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
while(~isempty(curPartition))
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
    while(curPartSet>0)
        if(backtracking==false)%If we just entered this level.
            [theCombo,allRecurData{curPartSet}]=getNextCapacitatedCombo(partSizes(curPartSet),numReps(curPartSet),numIdx2Assign(curPartSet));

            assignCapacitatedCombo(theCombo);

            if(curPartSet==numUniqueParts)
                %Save the current set partition.
                setMPartList(:,curSetMPartition)=curSetPartition;
                curSetMPartition=curSetMPartition+1;
                
                %Repeat the current level but get into backtracking mode,
                %because after finishing this level, it will be necessary
                %to backtrack.
                backtracking=true;
                continue
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
                setMPartList(:,curSetMPartition)=curSetPartition;
                curSetMPartition=curSetMPartition+1;
                continue;
            end

            %Go to the next level
            curPartSet=curPartSet+1;
            backtracking=false;
            continue;
        end
    end

    curPartition=getNextMPartition(curPartition);  
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
