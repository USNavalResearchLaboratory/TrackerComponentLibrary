function [jobs,tStarts]=scheduleMinLatenessDense(deadlines,durations,tStart)
%%SCHEDULEMINLATENESSDENSE Given an ordered list of scalar deadlines for
%                   certain tasks and an ordered list of durations of those
%                   tasks, create a schedule of when to execute the tasks
%                   so as to minimize the lateness of the tasks. Tasks are
%                   scheduled densely (no space between them).
%
%INPUTS: deadlines A numJobsX1 or 1XnumJobs array of deadline when each job
%                  must be completed. These could be, for example, seconds
%                  after some epoch.
%        durations A numJobsX1 or 1XnumJobs array of how long it takes to
%                  complete each job.
%           tStart An optional parameter indicating the scalar time when
%                  scheduling the jobs can commence. If this parameter is
%                  omitted or an empty matrix is passed, then tStart=0 is
%                  used.
%
%OUTPUTS: jobs A numJobsX1 list indicating the order in which the jobs
%              should be executed. If deadlines is empty, then an empty
%              matrix is returned.
%      tStarts A numJobsX1 list of when each job should begin being
%              executed. If deadlines is empty, then an empty matrix is
%              returned.
%
%The algorithm is a very simple greedy algorithm described in Chapter 4.3
%of [1]. The optimal order is just obtained by sorting the deadlines.
%
%EXAMPLE:
% deadlines=[14;6;8;9;9;15];
% durations =[3;3;2;1;4;2];
% tStart=0;
% [jobs,tStarts]=scheduleMinLatenessDense(deadlines,durations,tStart)
%The optimal ordering of jobs is [2;3;4;5;1;6]; and the corresponding start
%times are [0;3;5;6;10;13].
%
%REFERENCES:
%[1] J. Kleinberg and É. Tardos, Algorithm Design, 1st ed. Pearson, 2005.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(deadlines))
    jobs=[];
    tStarts=[];
    return;
end

%Times start from 0 if no initial time is given.
if(nargin<3||isempty(tStart))
    tStart=0;
end

numJobs=length(deadlines);

%The optimal order is just greedy --by earliest deadline. 
[~,jobs]=sort(deadlines,'ascend');
durations=durations(jobs);

tStarts=zeros(numJobs,1);

tCur=tStart;
for curJob=1:numJobs
    tStarts(curJob)=tCur;
    tCur=tCur+durations(curJob);
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
