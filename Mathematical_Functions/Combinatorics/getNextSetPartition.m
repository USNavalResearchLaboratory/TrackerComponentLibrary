function [partition,recurVals]=getNextSetPartition(param1,recurVals)
%%GETNEXTSETPARTITION Get the next way of partitioning a set of n unique
%                  items in lexicographic order. The partition is a
%                  length-n vector q that specifies which partition each
%                  item belongs to. The order of the items in the
%                  partitions is not important and the order of the
%                  partitions themselves is not important. If only n is
%                  passed, then the values for the first partition in the
%                  sequence are returned.
%
%INPUTS: param1 If only this parameter is passed, and nothing else, then
%               this is
%               n The number of items in the set to be partitioned. Given
%                 only this input, the function will return the values in
%                 the first partition.
%               If two inputs are provided, then param 1 is
%               partition The current set partition that should be updated
%                         to get the next set partition.
% recurVals A structure containing elements p and nc, which are necessary
%          to get the next item in the sequence (these values would be
%          supplied by the previous iteration). The elements of the
%          structure are
%          p A return value for the current partition; p(i) is the number
%            of elements in the ith class of the output partition for i=1
%            to nc.
%         nc A return value for the current partition; the number of
%            classes in the output partition (the number of subsets).
%
%OUTPUTS: partition The updated partition. If the final partition was
%                   passed, then partition will be an empty matrix. If only
%                   n was passed, then partition will be the first
%                   partition. Paritions are numbered in increasing order.
%                   That is, the first occurence of partition numbered i in
%                   partition is always before the first occurences of
%                   partitions numbered j>1.
%         recurVals The updated value in recurValues the function
%                   nextSetPartition can be called again to get the next
%                   set partition if partition was not empty. If
%                   recurVals.nc==n, then the returned partition is the
%                   final partition and a subsequent call will return an
%                   empty  matrix.
%
%The algorithm is based on NEXEQU in Chapter 11 of [1]. There is a total of
%BellNumber(n) set partitions for a particular n.
%
%Given a set (1,2,3), the possible partitions and the corresponding values
%of q are
%(1,2,3)        partition=[1;1;1]
%(1,2)(3)       partition=[1;1;2]
%(1,3)(2)       partition=[1;2;1]
%(2,3)(1)       partition=[1;2;2]
%(1)(2)(3)      partition=[1;2;3]
%This function always returns q so that the partition numbers increase.
%That is, there are no jumps. For example, q=[2;1;2] is equivalent to
%q=[1;2;1], but this function will only return the latter.
%
%%The function is either called as
%[q,recurVals]=getNextSetPartition(n);
%to get the first set parition or as
%[q,recurVals]=getNextSetPartition(q,recurVals);
%to get subsequent set partitions.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If only n was passed, get the first set partition.
if(nargin==1)
    n=param1;
    recurVals.p(1)=n;
    recurVals.nc=1;
    partition=ones(n,1);
    return;
end

partition=param1;
n=length(partition);

p=recurVals.p;
nc=recurVals.nc;

%If the final set partition was passed, return the empty matrix.
if(nc==n)
   partition=[];
   return;
end

%Step B.
m=n;

%Step C
while(1)
    L=partition(m);
    if(p(L)~=1)
        break;
    end
    partition(m)=1;
    m=m-1;
end

%Step D
nc=nc+m-n;
p(1)=p(1)+n-m;
if(L==nc)
    nc=nc+1;
    p(nc)=0;
end

%Step E
partition(m)=L+1;
p(L)=p(L)-1;
p(L+1)=p(L+1)+1;

recurVals.p=p;
recurVals.nc=nc;

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
