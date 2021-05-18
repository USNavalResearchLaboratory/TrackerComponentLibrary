function partList=genAllSetPartitions(n)
%%GENALLSETPARTITIONS Generate all ways of partitioning a set of n unique
%                  items. A partition is a length-n vector that specifies
%                  which partition each item belongs to. The order of the
%                  items in the partitions is not important and the order
%                  of the partitions themselves is not important.
%
%INPUTS: n The integer number of unique items in the set to partition; n>=1.
%
%OUTPUTS: partList The nXnumPart list of set partitions. The indexation
%                  begins at 1. for the ith partition partList(k,i) is the
%                  index of the set to which item i belongs.
%
%This function implements Algorithm H of Section 7.2.1.5 of [1]. The number
%of set partitions is given by BellNumber(n).
%
%EXAMPLE:
%A set of 3 objects is partitioned as 
% partList=genAllSetPartitions(3)
%One will get
% partList=[1,1,1,1,1;
%           1,1,2,2,2;
%           1,2,1,2,3];
%These corresponds to the partitions of a set of elements (1,2,3) as:
%(1,2,3)
%(1,2)(3)
%(1,3)(2)
%(1)(2,3)
%(1)(2)(3)
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numSetPartitions=BellNumber(n);

partList=zeros(n,numSetPartitions);

%Step H1, allocate and initialize.
a=zeros(n,1);
b=ones(n-1,1);
m=1;

for curPart=1:numSetPartitions
    %Step H2
    partList(:,curPart)=a;
    
    if(a(n)~=m)
        %Step H3
        a(n)=a(n)+1;
        continue;
    end
    
    %Step H4
    j=n-1;
    while(a(j)==b(j))
        j=j-1;
    end
    
    %Step H5
    if(j==1)
        partList=partList+1;%Change indexation to 1.
        return;
    end
    a(j)=a(j)+1;
    
    %Step H6
    m=b(j)+(a(j)==b(j));
    j=j+1;
    while(j<n)
        a(j)=0;
        b(j)=m;
        j=j+1;
    end
    a(n)=0;
end

%It only gets here if n=1.
partList=partList+1;

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
