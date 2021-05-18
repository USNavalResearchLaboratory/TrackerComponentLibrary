function [cycles,numInCycle,n]=getPermutationCycles(perm)
%%GETPERMUTATIONCYCLES Given a permutation, get the cycles that make up the
%                      permutation. These can be viewed as instructions on
%                      how to swap elements to make a permutation. See
%                      details below.
%
%INPUTS: perm An NX1 or 1XN permutation fo the values 1:n.
%      cycles A maxCycleLengthXnumCycles matrix where the values of each
%             cycle are given as column vectors. Since cycles can be of
%             differing lengths, the number of rows is the maximum cycle
%             length.
%    numInCycle A numCyclesX1 vector that holds the length of each cycle in
%             cycles. For the ith cycle, only elements
%             cycles(1:numInCycle(i),i) are used in the cycle. The other
%             elements of cycles are zero.
%           n This just repreats the length of perm.
%
%An explanation of permutation cycles is in Chapter 1.3.3 of [1]. The basic
%idea is that each permutation moves a values from one position to another.
%An example in [1] is the permutation [[3,4,6,2,5,1]. One can write this
%below the indices of the positions as
% 1,2,3,4,5,6
% 3,4,6,2,5,1
%This can be viewed as instructions to swap elements. We start, 1 goes to
%3. Well, 3 goes to 6 and 6 goes to 1. Following a chain of those
%instructions is a cycle. Collecting all such chains collects all such
%cycles. Thus, this function is by following the cycles, marking each
%visited node (by flipping the sign in perm), so that we know when we have
%completed a full cycle (we encounter a negative value).
%
%If one has more cycles than are needed to minimally represent a
%permutation, then one can string along a great many cycles, run
%multiplyPermCycles and then run this function to get a simpler cycle
%representation.
%
%EXAMPLE:
%The aforementioned example above is
% [cycles,numInCycle,n]=getPermutationCycles([3,4,6,2,5,1]);
%One will get three cycles, one of them singleton (meaning that the value
%stayed in its original position). The results are:
% cycles=[1, 2, 5;
%         3, 4, 0;
%         6, 0, 0];
% numInCycle=[3,2,1];
%So, there are three cycles, which are (1,3,6), (2,4) and (5). The function
%multiplyPermCycles will provide a permutation for a given set of cycles.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(perm);

cycles=zeros(n,n);
numInCycle=zeros(n,1);

maxCycleLength=0;
numCycles=0;

while(1)
    curIdx=1;
    while(curIdx<=n)
        %We have found an unscanned cycle.
        if(perm(curIdx)>0)
            break;
        end
        curIdx=curIdx+1;
    end

    if(curIdx>n)
        %All cycles have been found. 
        %Size to fit.
        cycles=cycles(1:maxCycleLength,1:numCycles);
        numInCycle=numInCycle(1:numCycles);
        return;
    end

    numCycles=numCycles+1;
    curCycleLen=0;
    %We have the start of the cycle. We just need to trace through it.
    while(perm(curIdx)>0)
        curCycleLen=curCycleLen+1;
        cycles(curCycleLen,numCycles)=curIdx;
        perm(curIdx)=-perm(curIdx);%Mark as visited.
        curIdx=-perm(curIdx);
    end
    numInCycle(numCycles)=curCycleLen;
    maxCycleLength=max(numInCycle,maxCycleLength);
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
