function jobsSelected=scheduleIntervals(intervals)
%%SCHEDULEINTERVALS Given a set of jobs having fixed time intervals in
%           which one would like to execute them, determine the maximum
%           subset of jobs that are mutually compatible (that do not
%           overlap in time).
%
%INPUTS: intervals A 2XnumIntervals matrix of intervals such that
%                  intervals(:,i) is the ith interval where
%                  intervals(1,i)<=intervals(2,i).
%
%OUTPUTS: jobsSelected A numSelectedX1 set of the intervals that form the
%                  maximum subset of jobs to be executed that do not
%                  overlap. If intervals is an empty matrix, then
%                  jobsSelected is an empty matrix.
%
%A simple greedy algorithm is optimal. The algorithm is described in
%Chapter 4.1 of [1].  
%
%EXAMPLE:
% intervals=[0,1,3,3,4,5,6,8;
%            6,4,5,8,7,9,10,11];
% jobsSelected=scheduleIntervals(intervals);
%The maximum subset of jobs that can be scheduled is jobsSelected=[2;5;8].
%
%REFERENCES:
%[1] J. Kleinberg and É. Tardos, Algorithm Design, 1st ed. Pearson, 2005.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(intervals))
   jobsSelected=[];
   return;
end

numIntervals=size(intervals,2);

%Intervals are sorted by finishing time
[~,idx]=sort(intervals(2,:),'ascend');
intervals=intervals(:,idx);

%Allocate the maximum possible amount of space
jobsSelected=zeros(numIntervals,1);

jobsSelected(1)=idx(1);

%The time when the last job selected finished.
lastFinish=intervals(2,1);
numJobsSelected=1;

for curInterval=2:numIntervals
    if(intervals(1,curInterval)>=lastFinish)
        numJobsSelected=numJobsSelected+1;
        jobsSelected(numJobsSelected)=idx(curInterval);
        lastFinish=intervals(2,curInterval);
    end
end

%Shrink to fit the actual number of intervals assigned.
jobsSelected=jobsSelected(1:numJobsSelected);

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
