function [house4Agent,agent4House]=stableHouseAlloc(prefLists,houseList)
%%STABLEHOUSEALLOC This solves the stable house allocation problem with
%                  seniority. Given a set of agents who rank a collection
%                  of houses, we want to find a stable allocation of houses
%                  to agents. That is, no pair of agents would both prefer
%                  to swap houses versus keeping the one they are in.
%                  However, the notion of seniority comes in when an agent
%                  has a house and wants to trade up. In such an instance,
%                  the same type of assignment will take place, but where
%                  the agents with seniority are guaranteed not to house a
%                  house they rank lower than the current one they are
%                  occupying. This is a type of 2D assignment problem.
%                  
%INPUTS: prefLists A numAgentsXnumPrefs list of where
%                  prefLists(curAgent,curRank) is the number of the house
%                  (starting from 0) ranked in position curRank by the
%                  agent. Agents do not have to rank all houses. However,
%                  unranked houses are assumed to be undesired and agents
%                  not raking enough houses can be left without a house
%                  despite availability. All agents must rank at least one
%                  house. If agents do not all rank the same number of
%                  houses, all houses after the number ranked by the agent
%                  must be set to zero. Houses are specified by index
%                  starting from 1 to numHouses.
%        houseList A numHousesX1 or 1XnumHouses vector where
%                  houseList(numHouse) is the index of the agent currently
%                  residing in the house (and thus posessing seniority so
%                  that he is guaranteed not to get a house worse than his
%                  current house). Unoccupied houses are set to 0.
%
%OUTPUTS: house4Agent A numAgentsX1 vector where house4Agent(curAgent) is
%                     the index of the house to which the agent is
%                     assigned, and 0 otherwise.
%         agent4House A numHousesX1 vector where agent4House(curHouse)
%                     provides the index of the agent assigned to the house
%                     or 0 if the house is unassgined.
%
%This implements the top trading cycles algorithm as described in [1]. The
%algorithm has been slightly modified to handle preference lists of
%different lengths.
%
%Example from [1]:
% prefLists=[2,6,5,1,4,3,7;
%            7,1,6,5,4,3,2;
%            2,1,4,7,3,6,5;
%            2,4,3,6,1,7,5;
%            4,3,7,1,2,5,6];
% houseList=[1;
%            2;
%            3;
%            4;
%            0;
%            0;
%            0];
% [house4Agent,agent4House]=stableHouseAlloc(prefLists,houseList)
%And the optimal assignment of houses to agents is house4Agent=[2;7;1;4;3]
%
%Consider an example where there are not enough houses. In such an
%instance, we can see those with seniority getting preference:
% prefLists=[1,2,3,4;
%            2,3,4,1;
%            4,3,2,1;
%            1,2,3,4;
%            3,1,2,4];
% houseList=[1;3;2;4];
% [house4Agent,agent4House]=stableHouseAlloc(prefLists,houseList)
%Here, the stable assignment is house4Agent=[1;3;4;2;0]; The last agent
%cound not get a house, because he did not already have seniority.
%
%The effects of seniority can be seen with two more example using the same
%prefLists as above:
% houseList=[2;3;0;1];
% [house4Agent1,agent4House1]=stableHouseAlloc(prefLists,houseList)
% houseList=[2;3;5;1];
% [house4Agent2,agent4House2]=stableHouseAlloc(prefLists,houseList)
%In the first example, the stable assignment is house4Agent1=[1;2;4;3;0],
%leaving the last agent out again. However, in the second example, the
%stable assignment is house4Agent =[1;2;4;0;3]; as the last agent has
%seniority and cannot be booted out.
%
%The algorithm can also handle short preference lists, which might leave
%extra agents unassigned. In the previous case:
% prefLists=[1,2,3;
%            2,3,4;
%            4,3,2;
%            1,2,3;
%            3,1,2];
% houseList=[1;3;2;4];
% [house4Agent,agent4House]=stableHouseAlloc(prefLists,houseList)
%Here, we get an assignment of house4Agent=[1;3;4;2;0], which is the same
%as before. However, if we make the preference lists too short and they are
%bit conflicting, then some agents who would normally have seniority can
%get stuck with nothing:
% prefLists=[1;
%            2;
%            4;
%            1;
%            1];
% houseList=[1;3;2;4];
% [house4Agent,agent4House]=stableHouseAlloc(prefLists,houseList)
%The stable assignment here is house4Agent=[1;2;4;0]. In this instance,
%there are enough houses, but the last two agents only want house 1, which
%must go to the first agent due to seniority, so they get nothing.
%
%REFERENCES:
%[1] A. Abdulkadiro?lu�lu and T S�nmez, "House Allocation with Existing
%    Tenants," Journal of Economic Theory, vol. 88, no. 2, pp. 233-260,
%    Oct. 1999.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
numAgents=size(prefLists,1);
numPrefs=size(prefLists,2);
numHouses=size(houseList,1);

%houseList(houseIdx) indicates which agent is currently occupying the
%house. We will make an inverse list so that invHouseList(agentIdx), so
%that we can determine which, if any house an agent is currently occupying.
invHouseList=zeros(numAgents,1);
for curHouse=1:numHouses
    if(houseList(curHouse)~=0)
        invHouseList(houseList(curHouse))=curHouse;
    end
end

%prefLists(agent,rank)=house gives the house of the given rank for the
%given agent. We will make a set of inverse pref lists such that
%invPrefLists(agent,house)=rank. Houses that are not in a preference list
%get assigned 0 in the inverse list.
invPrefLists=zeros(numAgents,numHouses);
for curAgent=1:numAgents
    for curRank=1:numPrefs
        curHouse=prefLists(curAgent,curRank);
        %This is needed in case some agents do not rank all of the houses
        %and thus put in zeros.
        if(curHouse~=0)
            invPrefLists(curAgent,curHouse)=curRank;
        end
    end
end

%Rather than deleting agents and houses from prefLists as they are assigned,
%which would involve computationally expensive matrix manipulations, we
%just adjust linked lists. The first pair of arrays is for the agent (rows
%of the matrix).
nextAgentList=[2:(numAgents),0];
prevAgentList=0:(numAgents-1);

%Next, we have matrices for the columns of the matrix. Matrices are
%necessary, as things  in different positions of each list need to be
%deleted.
nextPrefList=repmat([2:(numPrefs),0],[numAgents,1]);
prevPrefList=repmat(0:(numPrefs-1),[numAgents,1]);

%Now, see if any of the preference lists are short(contain zeros). If so,
%then that will be an early end to the list, which should be marked in
%nextPrefList
for curAgent=1:numAgents
    foundIdx=find(prefLists(curAgent,:)==0,1);
    
    if(~isempty(foundIdx))
        nextPrefList(curAgent,foundIdx-1)=0;
    end
end

firstAgent=1;
firstPrefList=ones(numAgents,1);

%Allocate space for return variables. Also, if agents or houses are
%unassigned in the end, then these correctly hold zeros.
agent4House=zeros(numHouses,1);
house4Agent=zeros(numAgents,1);

cycle=findCycle(prefLists,houseList,firstAgent,firstPrefList);
while(~isempty(cycle))
    numInCycle=size(cycle,2);
    
    for curInCycle=1:numInCycle
        theAgent=cycle(1,curInCycle);
        theHouse=cycle(2,curInCycle);
        
        %If we have an agent that is not assigned to a house, because his
        %list is too short, then he is just removed and not assinged to
        %anything (assigned to 0). In such an instance, numInCycle should
        %equal one.
        if(theHouse~=0)
            agent4House(theHouse)=theAgent;
            house4Agent(theAgent)=theHouse;
        end
        
        %Remove the agent from the preference list.        
        %Update the doubly linked list.
        if(prevAgentList(theAgent)~=0)
            nextAgentList(prevAgentList(theAgent))=nextAgentList(theAgent);
        else%If it was the first thing in the list
            firstAgent=nextAgentList(theAgent);
        end
        if(nextAgentList(theAgent)~=0)
            prevAgentList(nextAgentList(theAgent))=prevAgentList(theAgent);
        end
        
        %Mark any previous house the agent might have been occupying as
        %vacant. Note that it does not matter if this is the house the
        %agent is assigned in the cycle.
        if(invHouseList(theAgent)~=0)
            houseList(invHouseList(theAgent))=0;
        end
        
        %Remove the assigned house from the preference lists of all of the other agents, assuming
        %that the agent was assigned to a house and is not just eliminated for having too short a preference list.
        if(theHouse~=0)
            curAgent=firstAgent;
            while(curAgent~=0)
                curRank=invPrefLists(curAgent,theHouse);
            
                %Eliminate the house from the preference list of the
                %current agent.
                if(curRank~=0)
                    %Now, update the rest of the doubly linked list.
                    if(prevPrefList(curAgent,curRank)~=0)
                        nextPrefList(curAgent,prevPrefList(curAgent,curRank))=nextPrefList(curAgent,curRank);
                    else%If it was the first thing in the list.
                        firstPrefList(curAgent)=nextPrefList(curAgent,curRank);
                    end
                    if(nextPrefList(curAgent,curRank)~=0)
                        prevPrefList(curAgent,nextPrefList(curAgent,curRank))=prevPrefList(curAgent,curRank);
                    end
                end
                
                curAgent=nextAgentList(curAgent);
            end
        end
    end
    cycle=findCycle(prefLists,houseList,firstAgent,firstPrefList);
end
end

function cycle=findCycle(prefLists,houseList,firstAgent,firstPrefList)
    numAgents=size(prefLists,1);
    numHouses=size(houseList,1);

    if(firstAgent==0)
       cycle=[];
       return;
    end
    
    %The maximum possible cycle length.
    maxObj=min(numAgents,numHouses);
    
    %Allocate maximum possible space for the cycle...
    cycle=zeros(2,maxObj);
    
    cycleLength=0;
    
    %Start with the first agent that has not been removed form the selection process
    curAgent=firstAgent;
    while(1)
        %If the agent's list was too short, so there are no more assignable houses that are
        %acceptable to him, then return him as the cycle with zero for the assignment
        %indicating that he is to be removed as unassignable.
        if(firstPrefList(curAgent)==0)
            cycle=[curAgent;0];
            return;
        end
        
        %Index of the house top-ranked by the current agent.
        houseIdx=prefLists(curAgent,firstPrefList(curAgent));
        
        %Add the current assignment to the cycle.
        cycleLength=cycleLength+1;
        cycle(1,cycleLength)=curAgent;
        cycle(2,cycleLength)=houseIdx;
        
        %Get the next agent index. If the house is occupied, then this points to the owner.
        %Otherwise, this is zero and the next agent is the lowest index agent.
        if(houseList(houseIdx)==0)
            %Vacant house, so we return to the first one in the list.
            cycle=cycle(:,1:cycleLength);
            return;
        else
            curAgent=houseList(houseIdx);%Occupied house
        end

        %If a cycle is encountered, the return the cycle.
        cycleStartIdx=find(curAgent==cycle(1,:));
        
        if(~isempty(cycleStartIdx))
            cycle=cycle(:,cycleStartIdx:cycleLength);
            return;
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
