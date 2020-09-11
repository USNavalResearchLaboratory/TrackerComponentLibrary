function a=getNextMPartition(param1,m)
%%GETNEXTMPARTITION Get the next partition of the integer n into m
%                   parts. This is m nonzero items such that the sum is n.
%                   The partitions are visited in reverse lexicographic
%                   order (colex order). The function can be called with
%                   two parameters to get the first n,m partition, or with
%                   1 parameter (that being the previous partition) to
%                   get the next partition. The total number of
%                   partitions can be found using the function
%                   numMPartitions(n,m). As combinations are to
%                   permutations, partitions are to compositions: The order
%                   of the terms in an m-partitions do not matter, whereas
%                   they matter in a composition.
%
%INPUTS: param1 If two parameters are passed to the function, then param1
%            is n, the integer to be partitioned. Otherwise, param1 is a,
%            the current partition of the chosen integer n into m parts.
%            This means that a must be length m and its elements must sum
%            to n. The elements must be ordered in decreasing order. The
%            first partition in the series must be a(1)=n-m+1; a(2:m)=1;
%            Note that n>=m>=2 if a is passed, because a size 1 partition
%            only has one value (the initial value).
%          m This parameter is not passed if the next m-partition in the
%            sequence is desired. However, if one wishes to obtain the
%            first value in the sequence, then m is the positive integer
%            number of parts into which the integer n is split, n>=m;
%
%OUTPUTS: a The next partition in the series or the empty matrix if the
%           parition given was the last partition in the series. If the
%           function was called getNextMPartition(n,m), then a is the first
%           parition.
%
%The algorithm is Algorithm H in Chapter 7.2.1.4 of [1].
%
%The function can either be called as
%a=getNextMPartition(n,m);
%to get the first parition of the integer n into m parts or as
%a=getNextMPartition(a);
%to get the next partition.
%
%Often, one might prefer all partitions of m integers that sum to a
%constant value n INCLUDING zeros. This function can be used to produce
%such a modified set of partitions by initializing it as
%a=getNextMPartition(n+m,m);
%and for every composition generated, including the first, the modified
%composition starting from zero is given by a-1.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Return the first partition in the sequence if two parameters are passed.
if(nargin>1)
    n=param1;
    a(1)=n-m+1;
    a(2:m)=1;
    return;
end

a=param1;

m=length(a);

%Deal with the special case of a 1D partition. In this instance, there is
%only 1 possible solution, so if one is calling this function again, it
%will be past the end of all possible solutions.
if(m==1)
    a=[];
    return;
end

while(a(2)<a(1)-1)
    %Step H3 in Knuth
    a(1)=a(1)-1;
    a(2)=a(2)+1;
    return;
end

%Step H4 in Knuth, with an added check to avoid reading past the end of the
%array if m=2.
j=3;
if(j>m)
    a=[];
    return;
end
s=a(1)+a(2)-1;
while(a(j)>=a(1)-1)
    s=s+a(j);
    j=j+1;
    %If we have passed the last partition
    %(Termination part of step H5)
    if(j>m)
       a=[];
       return;
    end
end

%Step H5 in Knuth.
x=a(j)+1;
a(j)=x;
j=j-1;

%Step H6 in Knuth
while(j>1)
    a(j)=x;
    s=s-x;
    j=j-1;
    a(1)=s;
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
