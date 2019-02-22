function theRank=rankPartition(thePartition)
%%RANKPARTITION Rank a given partition of the integer n when considering
%               the partitions given in reverse lexicographic order. A
%               partition of the integer n is a set of integers that sum to
%               n.
%
%INPUTS: thePartition An sX1 or 1Xs partition. The values start at 1.
%
%OUTPUTS: the Rank. The integer rank of the partition starting from 1.
%
%Let p(n,m) be the number of partitions of n with up to m parts. It is
%known that p(n,m)=p(n,m-1)+p(n-m,m). This relation is given in Chapter
%7.2.1.4 of [1]. The rank generated is derived using that recurrence
%relation.
%
%EXAMPLE:
%The function getNextPartition generates partitions in reverse
%lexicographic order. Here, we demonstrate that the partitions generated
%are ranked correctly. We partition the integer 6.
% n=6;
% [thePartition,theData]=getNextPartition(n);
% while(~isempty(thePartition))
%     theRank=rankPartition(thePartition)
%     [thePartition,theData]=getNextPartition(n,theData);
% end
%One will see that theRank goes from 1 to 11 in order.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    m=length(thePartition);
    n=sum(thePartition);
    
    thePartition=sort(thePartition,'descend');

    theRank=0;
    for idx=1:m
        curParts=thePartition(idx)-1;
        if(curParts>0)
            theRank=theRank+numMPartitions(sum(thePartition(idx:m)),curParts,true);
        end
    end
    theRank=numberOfPartitions(n)-theRank;
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
