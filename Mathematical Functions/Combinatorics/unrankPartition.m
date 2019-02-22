function thePartition=unrankPartition(rank,n)
%%UNRANKPARTITION Given the rank of a partiion of the integer n, when
%                 considering the partitions given in reverse lexicographic
%                 order, find the partition of the given rank. A partition
%                 of the integer n is a set of integers that sum to n.
%
%INPUTS: rank The integer rank of the partition. rank>0.
%           n The integer being partitioned.
%
%OUTPUTS: thePartition A numElementsX1  partition of n of the given rank.
%
%Let p(n,m) be the number of partitions of n with up to m parts. It is
%known that p(n,m)=p(n,m-1)+p(n-m,m). This relation is given in Chapter
%7.2.1.4 of [1]. The rank obtained here is derived using that recurrence
%relation.
%
%EXAMPLE:
%Here, we demonstrate that the unrankPartition function is consistent with
%the rankPartition function.
% n=20;
% [thePartition,theData]=getNextPartition(n);
% while(~isempty(thePartition))
%     assert(all(unrankPartition(rankPartition(thePartition),n)==thePartition))
%     [thePartition,theData]=getNextPartition(n,theData);
% end
%There should be no error indicating that the unranked partitions are all
%the same as the ranked partitions.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

thePartition=zeros(n,1);%We do not know how many parts to have...

rank=numberOfPartitions(n)-rank;

if(rank<0)
    error('The rank given is more than the total number of partitions of n')
end

curIdx=0;
nLeft=n;
while(rank>0)
    parts=1;
    numParts=numMPartitions(nLeft,parts,true);
    while(numParts<rank)
        parts=parts+1;
        numParts=numMPartitions(nLeft,parts,true);
    end
    if(numParts>rank)
        parts=parts-1;
        numParts=numMPartitions(nLeft,parts,true);
    end
    
    curIdx=curIdx+1;
    thePartition(curIdx)=parts+1;
    rank=rank-numParts;
    nLeft=nLeft-thePartition(curIdx);
end

while(nLeft>0)
    curIdx=curIdx+1;
    thePartition(curIdx)=1;
    nLeft=nLeft-1;
end

%Shrink to fit.
thePartition=thePartition(1:curIdx);

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
