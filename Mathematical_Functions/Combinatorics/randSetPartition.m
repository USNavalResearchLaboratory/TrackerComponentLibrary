function [q,nc]=randSetPartition(n,sortEntries)
%%RANDSETPARTITION Generate a random partition of a set of n items.  The
%                  partition is a length-n vector q that specifies which
%                  partition each item belongs to. nc is the number of
%                  subsets in the partition q. The order of the items in
%                  the partitions is not important and the order of the
%                  partitions themselves is not important. See the function
%                  getNextSetPartition for partition generation.
%
%INPUTS: n The number of items that should be randomly partitioned.
% sortEntries An optional parameter specifying whether the partition
%          numbers should come in order of occurrence. For example,
%          q=[2;1;2] and q=[1;2;1] are the same partition of n=3. If
%          sortEntries is true, then only q=[1;2;1] would be returned out
%          of the two. The default if this parameter is omitted is false.
%
%OUTPUTS: q A random set partition of n. This is a length n vector with
%           values that can range from 1 to n.
%        nc The number of subsets in the random set partition.
%
%The algorithm is RANEQU in Chapter 12 of [1].
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    sortEntries=false;
end

%Allocate space.
b=zeros(n,1);
q=zeros(n,1);

b(1)=1;

for L=1:(n-1)
    sumVal=1/L;
    L1=L-1;
    if(L1~=0)
        for k=1:L1
           sumVal=(sumVal+b(k))/(L-k); 
        end
    end
    b(L+1)=(sumVal+b(L))/(L+1);
end

m=n;
L=0;
z=-1;
while(1)
    if(z<0)
        z=m*b(m)*rand(1);
        k=0;
        L=L+1;
    end
    q(m)=L;
    m=m-1;

    if(m==0)
        break;
    end

    z=z-b(m);
    k=k+1;
    z=z*k;
end
p=randperm(n);
q=q(p);

nc=L;

%If one wants the indices of the sets to start at 1 and go up by 1 when a
%new set is encountered. This is nice for display, but is usually not
%needed in general.
if(sortEntries)
    val=0;
    for idx=1:n
        if(q(idx)>val+1)
            val=val+1;
            
            theVal=q(idx);

            selIsVal=q(:)==val;
            selNotVal=q(:)==q(idx);

            q(selNotVal)=0;
            q(selIsVal)=theVal;
            q(selNotVal)=val;
        elseif(q(idx)==val+1)
            val=val+1;
        end
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
