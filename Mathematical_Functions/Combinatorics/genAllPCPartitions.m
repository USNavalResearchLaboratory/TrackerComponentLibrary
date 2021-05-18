function [cList,lList]=genAllPCPartitions(n)
%%GENALLPCPARTITIONS Generate all partitions of the integer n in part-count
%               format. That is, find all n such that
%               n=c(1)+2*c(2)+...n*c(n). The number of parts is sum(c)
%               (some elements of c will be zero). All element of c and
%               >=0.
%
%INPUTS: n The integer to partition. n>0.
%
%OUTPUTS: cList An nXnumPart list of the values of each part for each of
%               the numPart partitions.
%         lList An (n+1)XnumPart link table such that when considering the ith
%               partition, if c=cList(:,i) and l=lList(:,i), then if there
%               are t nonzero values in c, having indices k1,k2,...,k(t-1),
%               kt, with k1<k2<...<kt, then l(0+1)=k1, l(k1+1)=k2,
%               l(k(t-1)+1)=kt, l(kt+1)=0.
%               
%This function implements Problem 5 of Section 7.2.1.4 of [1]. The total
%number of partitions is numPart=numberOfPartitions(n).
%
%EXAMPLE:
%The partitions for n=5 are:
% [cList,lList]=genAllPCPartitions(5)
%where one will see that
%cList=[5,3,1,2,0,1,0;
%       0,1,2,0,1,0,0;
%       0,0,0,1,1,0,0;
%       0,0,0,0,0,1,0;
%       0,0,0,0,0,0,1];
%and 
%lList=1,1,1,1,2,1,5;
%      0,2,2,3,3,4,4;
%      0,0,0,0,3,3,3;
%      0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0];
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPart=numberOfPartitions(n);

cList=zeros(n,numPart);
lList=zeros(n+1,numPart);

%Allocate space. Indexation is from 0 in [1].
c=zeros(n+1,1);
l=zeros(n+1,1);

%Step C1
c(0+1)=1;
c(1+1)=n;
%All other c are zero.

l(0+1)=1;
%l(1) is already 0.

for curPart=1:numPart
    skipC6=false;
    %Step C2, visit the partition.
    cList(:,curPart)=c(2:(n+1));
    lList(:,curPart)=l;

    %Step C3, branch.
    j=l(0+1);
    k=l(j+1);
    if(c(j+1)~=1)
        if(j>1)
            %Step C5
            c(1+1)=j*(c(j+1)-1)-1;
            skipC6=true;
        else
            %Step C4
            c(1+1)=c(1+1)-2;
            c(2+1)=c(2+1)+1;
            l((c(1+1)>0)+1)=2;
            if(k~=2)
                l(2+1)=k;
            end
            continue;
        end
    end
    
    %Step C6
    if(skipC6==false)
        if(k==0)
            return;
        end

        c(j+1)=0;
        c(1+1)=k*(c(k+1)-1)+j-1;
        j=k;
        k=l(k+1);
    end

    %Step C7
    if(c(1+1)>0)
        l(0+1)=1;
        l(1+1)=j+1;
    else
        l(0+1)=j+1;
    end
    
    c(j+1)=0;
    c(j+1+1)=c(j+1+1)+1;
    if(k~=j+1)
        l(j+1+1)=k; 
    end
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
