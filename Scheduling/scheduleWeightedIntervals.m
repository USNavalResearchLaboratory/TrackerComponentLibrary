function [weightVal,jobsSelected]=scheduleWeightedIntervals(intervals,weights)
%%SCHEDULEWEIGHTEDINTERVALS Given a set of jobs having fixed time intervals
%                   in which one would like to execute them as well as
%                   positive weights associated with each job, determine
%                   the subset of jobs that are mutually compatible (that
%                   do not overlap in time) that maximize the total weight.
%
%INPUTS: intervals A 2XnumIntervals matrix of intervals such that
%                  intervals(:,i) is the ith interval where
%                  intervals(1,i)<=intervals(2,i).
%          weights A numIntervalsX1 or 1XnumIntervals vector of the weights
%                  associated with the intervals. Note that all entries of
%                  weights must be >=0 and finite.
%
%OUTPUTS: weightVal The scalar total weight of the subset of jobs having
%                   maximum weight.
%      jobsSelected A numSelectedX1 set of the intervals that form the
%                   maximum weight subset of jobs to be executed that do
%                   note overlap. If intervals is an empty matrix, then
%                   jobsSelected is an empty matrix.
%
%The algorithm to solve the problem is based on dynamic programming and is
%described in Chapters 6.1 and 6.2 of [1].
%
%EXAMPLE:
% intervals=[1,3,0,4,3,5,6,8;
%            4,5,6,7,8,9,10,11];
% weights=[1;2;3;4;5;6;7;8];
% [weightVal,jobsSelected]=scheduleWeightedIntervals(intervals,weights)
%The optimal weight value is 13 and the optimal assignment of jobs is
%jobsSelected=[8;4;1]. That is actually the same as the maximum subset of
%jobs obtained using scheduleIntervals without considering the weights.
%However,the weights matter. For example in 
% intervals=[1,3,0,4,3,5,6,8;
%            4,5,6,7,8,9,10,11];
% weights=[1;2e8;3;4;5;6;7;8];
% [weightVal,jobsSelected]=scheduleWeightedIntervals(intervals,weights)
%The high weight of the second job forces it to be in any optimal
%assignment. Thus, jobsSelected=[8;2];
%
%Also, one might want to check that the algorithm works when the intervals
%are not sorted in order of increasing finishing time. Thus, 
% intervals=[3,8, 1,3,6, 5,0,4;
%            5,11,4,8,10,9,6,7];
% weights=[2;8;1;5;7;6;3;4];
% [weightVal,jobsSelected]=scheduleWeightedIntervals(intervals,weights)
%Here, one again gets a weight value of 13, but the jobsSelected=[2;8;3].
%These are the same intervals selected in the  first example.
%
%REFERENCES:
%[1] J. Kleinberg and É. Tardos, Algorithm Design, 1st ed. Pearson, 2005.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(intervals))
    weightVal=[];
    jobsSelected=[];
    return
end

numIntervals=size(intervals,2);

%Intervals/weights are sorted by finishing time
[~,idx]=sort(intervals(2,:),'ascend');
intervals=intervals(:,idx);
weights=weights(idx);

%This will be the largest index of the sorted intervals such p(i) gives the
%largest index i<j such that job i is compatible with j. p(1) is just left
%set to 0.
p=zeros(numIntervals,1);
for curP=2:numIntervals
    %Perform a binary search for the value or the next lowest value.
    %we are finding the rightmost interval such that the finishing time of
    %that interval is <= the start time of interval curP.
    [~, foundIdx]=binSearch(intervals(2,1:(curP-1)),intervals(1,curP),1);
    %If the found index is compatible
    if(intervals(1,curP)>=intervals(2,foundIdx))
        p(curP)=foundIdx;
    else
        p(curP)=0;
    end
end

%Given the p values, we can use dynamic programming to solve the problem.
opt=zeros(numIntervals+1,1);
opt(0+1)=0;

for curWeight=1:numIntervals
    if(p(curWeight)==0)
        opt(curWeight+1)=max([weights(curWeight),opt(curWeight-1+1)]);
    else
        opt(curWeight+1)=max([weights(curWeight)+opt(p(curWeight)+1),opt(curWeight-1+1)]);
    end
end

weightVal=opt(end);

%If the actual optimal interval assignment is desired.
if(nargout>1)
    %Allocate the maximum possible amount of space.
    jobsSelected=zeros(numIntervals,1);
    numJobsAdded=0;
    
    curWeight=numIntervals;
    
    while(1)
        if(p(curWeight)==0)
            if(weights(curWeight)>opt(curWeight-1+1))
                numJobsAdded=numJobsAdded+1;
                jobsSelected(numJobsAdded)=idx(curWeight);
                break;
            else
                curWeight=curWeight-1;
            end
        else
            if(weights(curWeight)+opt(p(curWeight)+1)>opt(curWeight-1+1))
                numJobsAdded=numJobsAdded+1;
                jobsSelected(numJobsAdded)=idx(curWeight);
                curWeight=p(curWeight);
            else
                curWeight=curWeight-1;
            end
        end
    end
    
    %Shrink to the actual number selected
    jobsSelected=jobsSelected(1:numJobsAdded);
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
