function partners=assignStableRoommates(prefLists)
%%ASSIGNSTABLEROOMMATES Solve the basic form of the stable roommates 
%        problem. In this problem, an even number of individuals rank their
%        preferences for each other. The algorithm then finds a stable 
%        pairing of individuals, if one exists. An assignment is stable if
%        there is no pair of (non-assigned) members where BOTH prefer each 
%        other over their current assignment. This function finds one
%        solutions; though multiple solutions may exist. If no solutions
%        exist, an empty matrix is returned.
%
%INPUTS: prefLists An NX(N-1) matrix where prefLists(n,i) is the ith
%                  preferred roommate for person n. N must be an even
%                  number. The "roommate" is just a number from 1 to N that
%                  cannot equal n. Alternatively, an NXN matrix can be
%                  passed, where the last column is equal to the row
%                  number. This is allowed as this algorithm just augments
%                  any NX(N-1) matrix to have that format before proceeding
%                  (it plays a role in identifying infeasible assignments).
%
%OUTPUTS: partners An NX1 vector listing which each partner of the N
%                  individuals is assigned to share a room in a stable
%                  assignment (thus, partners is a permutation vector). If
%                  no stable assignment exists, then an empty matrix is
%                  returned.
%
%The O(N^2) algorithm of [1] is used to solve the problem.
%   
%Example with a solution (from [1]):
% prefLists=[4,6,2,5,3;
%            6,3,5,1,4;
%            4,5,1,6,2;
%            2,6,5,1,3;
%            4,2,3,6,1;
%            5,1,4,2,3];
% partners=assignStableRoommates(prefLists)
%The solution should be partners=[6;3;2;5;4;1] for pairs of (1,6), (2,3),
%and (4,5)
%
%Example without a solution (from [1]):
% prefLists=[2,6,4,3,5;
%            3,5,1,6,4;
%            1,6,2,5,4;
%            5,2,3,6,1;
%            6,1,3,4,2;
%            4,2,5,1,3];
% partners=assignStableRoommates(prefLists)
%The solution in this instance is an empty matrix, because no stable
%solution exists.
%   
%Example with multiple solutions (from [1]). Only 1 solution is found by
%this function, though.
% prefLists=[2,5,4,6,7,8,3;
%            3,6,1,7,8,5,4;
%            4,7,2,8,5,6,1;
%            1,8,3,5,6,7,2;
%            6,1,8,2,3,4,7;
%            7,2,5,3,4,1,8;
%            8,3,6,4,1,2,5;
%            5,4,7,1,2,3,6];
% partners=assignStableRoommates(prefLists)
%The solution that will be found by this algorithm is
%partners=[4;3;2;1;6;5;8;7]. However, the other two possible solutions are
%[1,2,3,4,8,7,6,5] and [5,6,7,8,1,2,3,4].
%  
%REFERENCES:
%[1] R. W. Irving, "An Efficient Algorithm for the 'Stable Rommates'
%    Problem," Journal of Algorithms, vol. 6, no. 4, pp. 577-595, Dec.
%    1985.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    numPeople=size(prefLists,1);
    
    %If the preference lists do not already contain "self" (sentinel) as
    %the last column, then add it.
    if(size(prefLists,2)<numPeople)
        prefLists=[prefLists,(1:numPeople)'];
    end
    
    firstUnmatched=1;
    firstInCycle=1;
    
    %Create the ranking matrix to quickly compare preferences and
    %initialize leftmost and rightmost matrices
    ranking=zeros(numPeople,numPeople);
    %Extra entries are sentinels for procedure
    %find1stwithMoteThan1Potential
    leftmost=ones(numPeople+1,1);
    rightmost=numPeople*ones(numPeople+1,1);
    for curPerson=1:numPeople
        for rank=1:numPeople
            ranking(curPerson,prefLists(curPerson,rank))=rank;
        end
    end
    
    [solutionPossible,leftmost,rightmost]=phase1Reduce(prefLists,ranking,leftmost,rightmost);
    
    %Proper initialization unnecessary
    second(1:numPeople)=leftmost(1:numPeople)+1;
    
    solutionFound=false;
    tail=zeros(numPeople,1);%Initially, nothing in tail (all zeros)
    cycle=zeros(numPeople,1);
    while(solutionPossible&&solutionFound==false)
        firstUnmatched=find1stwithMoteThan1Potential(firstUnmatched,leftmost,rightmost);
        if(firstUnmatched>numPeople)
            solutionFound=true;
            break;
        else
            [cycle,firstInCycle,lastInCycle,firstUnmatched,second,tail]=seekCycle(prefLists,ranking,cycle,firstInCycle,firstUnmatched,second,tail,rightmost);
            
            [solutionPossible,leftmost,second,rightmost]=phase2Reduce(prefLists,ranking,cycle,leftmost,second,rightmost,firstInCycle,lastInCycle);
        end
    end
    
    if(solutionFound)
        partners=zeros(numPeople,1);
        
        for curPerson=1:numPeople
            partners(curPerson)=prefLists(curPerson,leftmost(curPerson));
        end
    else
        partners=[];%No solution found.
    end
end

function [solutionPossible,leftmost,rightmost]=phase1Reduce(prefLists,ranking,leftmost,rightmost)
    numPeople=size(prefLists,1);
    
    %Initially, no one proposed to (all zeros).
    setProposedTo=zeros(numPeople,1);
    for curPerson=1:numPeople
        proposer=curPerson;
        while(1)
            %Best potential partner
            nextChoice=prefLists(proposer,leftmost(proposer));
            %nextChoice holds current
            current=prefLists(nextChoice,rightmost(nextChoice));
            
            while(ranking(nextChoice,proposer)>ranking(nextChoice,current))
                %Proposer is rejected by nextChoice 
                leftmost(proposer)=leftmost(proposer)+1;
                nextChoice=prefLists(proposer,leftmost(proposer));
                current=prefLists(nextChoice,rightmost(nextChoice));
            end
            %nextChoice holds proposer
            rightmost(nextChoice)=ranking(nextChoice,proposer);
            %and rejects current
            proposer=current;
            
            if(setProposedTo(nextChoice)==0)
                break;
            end
        end
        setProposedTo(nextChoice)=1;
    end
    
    solutionPossible=(proposer==nextChoice);
end

function firstUnmatched=find1stwithMoteThan1Potential(firstUnmatched,leftmost,rightmost)
%%Find the first person with more than one potential partner.
    while(leftmost(firstUnmatched)==rightmost(firstUnmatched))
        firstUnmatched=firstUnmatched+1;
    end
end

function [cycle,firstInCycle,lastInCycle,firstUnmatched,second,tail]=seekCycle(prefLists,ranking,cycle,firstInCycle,firstUnmatched,second,tail,rightmost)
    numPeople=size(prefLists,1);
    
    if(firstInCycle>1)
        %Last person in previous tail.
        person=cycle(firstInCycle-1);
        %His second choice may have to be updated
        posInCycle=firstInCycle-1;
        cycleSet=tail;
    else
        cycleSet=zeros(numPeople,1);%Nothing in the cycle set
        posInCycle=1;
        person=firstUnmatched;
    end
    
    %Generate sequence
    while(1)
        cycleSet(person)=1;%Add person to the cycle set.
        
        cycle(posInCycle)=person;
        posInCycle=posInCycle+1;
        posInList=second(person);
        
        while(1)%Update second choice for current person
            nextChoice=prefLists(person,posInList);
            posInList=posInList+1;
            
            if(ranking(nextChoice,person)<=rightmost(nextChoice))
               break;
            end; 
        end
        
        second(person)=posInList-1;
        person=prefLists(nextChoice,rightmost(nextChoice));
        
        %If sequence starts to cycle.
        if(cycleSet(person)==1)
            break;
        end
    end
    lastInCycle=posInCycle-1;
    tail=cycleSet;
    
    %Work back to the beginning of the cycle
    while(1)
       posInCycle=posInCycle-1;
       tail(cycle(posInCycle))=0;
        
       if(cycle(posInCycle)==person)
           break;
       end
    end
    firstInCycle=posInCycle;
end

function [solutionPossible,leftmost,second,rightmost]=phase2Reduce(prefLists,ranking,cycle,leftmost,second,rightmost,firstInCycle,lastInCycle)    
    for rank=firstInCycle:lastInCycle
        %Allow next person in cycle to be rejected
        proposer=cycle(rank);
        leftmost(proposer)=second(proposer);
        %Proper update unnaecessary at this stage
        second(proposer)=leftmost(proposer)+1;
        
        nextChoice=prefLists(proposer,leftmost(proposer));
        rightmost(nextChoice)=ranking(nextChoice,proposer);
        %nextChoice holds proposer.
    end
    
    rank=firstInCycle;
    solutionPossible=true;
    while(rank<=lastInCycle&&solutionPossible)
        %Check that no one has run out of potential partners
        proposer=cycle(rank);
        solutionPossible=(leftmost(proposer)<=rightmost(proposer));
        rank=rank+1;
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
