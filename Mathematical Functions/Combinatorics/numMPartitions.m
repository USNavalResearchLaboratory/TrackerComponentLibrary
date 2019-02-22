function val=numMPartitions(n,m,getUpToM)
%%NUMMPARTITIONS Get the number of ways of partitioning the integer n into
%                m parts, where order does not matter (unlike compositions,
%                where order matters). That is the number of ways of
%                getting m nonzero integers that sum to n. The amount of
%                memory for this implementation is proportional to n*m and
%                thus the solution becomes slow for very large n and m
%                values.
%
%INPUTS: n The positive n>0 integer to be partitioned.
%        m The number of parts into which the integer is to be partitioned.
% getUpToM This parameter changes what is returned. If false (the default
%          if omitted or an empty matrix is passed), the number of
%          ways of partitioning n into m parts is returned. If true, then
%          the total number of ways of partitioning n into m or fewer parts
%          is returned.
%
%OUTPUTS: val The number of ways of partitioning the integer n into m
%             parts or if getUpToM=true, then this is the total number of
%             ways of partitioning n into m or fewer parts.
%
%The counting method is based on the recurrence relation in Chapter 7.2.1.4
%of [1]. The algorithm could be implemented to use less memory. However,
%that would require shifting rows in a table up (to reuse memory) and would
%make the code a lot more complicated.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(getUpToM))
    getUpToM=false;
end

P=eye(n+1,m+1);
P(:,1+1)=1;

for curN=2:n
    for curK=1:min(curN,m)
        P(curN+1,curK+1)=P(curN-1+1,curK-1+1)+P(curN-curK+1,curK+1);
    end
end
if(getUpToM)
    val=sum(P(end,:));
else
    val=P(end,end);
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
