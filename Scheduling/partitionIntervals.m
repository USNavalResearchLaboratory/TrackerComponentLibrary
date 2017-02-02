function partition4Interval=partitionIntervals(intervals)
%%PARTITIONINTERVALS Given a set of possible overlkapping intervals, find a
%               way of arranging them into a minimum number of slots such
%               that no two intervals overlap. This could be used, for
%               example, to schedules classes into rooms when professors
%               set the times of the classes (assuming all rooms are large
%               enough).
%
%INPUTS: intervals A 2XnumIntervals matrix of intervals such that
%                  intervals(:,i) is the ith interval where
%                  intervals(1,i)<=intervals(2,i).
%
%OUTPUTS: partition4Interval A numIntervalsX1 vector where partitions(i)
%                  gives the integer partition to which the ith interval
%                  belongs.If intervals is an empty matrix, then
%                  partition4Interval will be an empty matrix too.
%
%The algorithm is a simple greedy algorithm as described in Chapter 4.1 of
%[1].
%
%EXAMPLE:
% intervals=[9,   9,   9,    11,  11, 13,  13,  14,  15,  15;
%            10.5,12.5,10.5, 12.5,14, 14.5,14.5,16.5,16.5,16.5];
% partition4Interval=partitionIntervals(intervals)
%The optimal set of partitions is partition4Interval=[1;2;3;1;3;1;2;3;1;2].
%
%REFERENCES:
%[1] J. Kleinberg and É. Tardos, Algorithm Design, 1st ed. Pearson, 2005.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(intervals))
    partition4Interval=[];
    return;
end

numIntervals=size(intervals,2);

%Sort intervals in order of increasing start time.
[~,idx]=sort(intervals(1,:),'ascend');
intervals=intervals(:,idx);
%revIdx is to undo the ordering.
[~,revIdx]=sort(idx,'ascend');

%There can be at most numIntervals partitions.
partition4Interval=zeros(numIntervals,1);
finishTimes=zeros(numIntervals,1);

%The first interval goes into the first partition
partition4Interval(1)=1;
finishTimes(1)=intervals(2,1);%The finish time
numPartitions=1;
for curInterval=2:numIntervals
    %Find the first partition that the current interval is compatible with.
    foundPart=0;
    for curPart=1:numPartitions
        if(finishTimes(curPart)<=intervals(1,curInterval))
            foundPart=curPart;
            break;
        end
    end
    
    %Must add another partition as it does not fit in any of the existing
    %ones.
    if(foundPart==0)
        numPartitions=numPartitions+1;
        foundPart=numPartitions;
    end
    %Schedule the interval in the found partition
    partition4Interval(curInterval)=foundPart;
    finishTimes(foundPart)=intervals(2,curInterval);
end

partition4Interval=partition4Interval(revIdx);

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
