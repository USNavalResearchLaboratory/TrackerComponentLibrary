function [rank,WT]=rankSetPartition(q,WT)
%%RANKSETPARTITION Return the set partition with the given rank in the
%           lexicographic ordering of all set parititions having length(n).
%           The partition is a length-n vector q that specifies which
%           partition each item belongs to. The order of the items in the
%           partitions is not important.
%
%INPUTS: q The nX1 set partition vector that should be ranked. Paritions
%          are numbered in increasing order starting from 1. That is, the
%          first occurence of partition numbered i in q is always before
%          the first occurences of partitions numbered j>1.
%       WT The is, optionally, the Williamson table. It is the table T in
%          [1]. The table is a function of length(q). If omitted, then the
%          table will be computed and returned so that futue calls to the
%          function can use 
%
%OUTPUTS: rank The rank of the set partition in lexicographic order.
%              0<=rank<BellNumber(n).
%           WT The Williamson table for a length-n set partition. Passing
%              this to rankSetPartition on subsequent calls will speed up
%              the algorithm.
%
%This function implements the ranking algorithm of [1]. See
%getNextSetPartition for more information on set partitions.
%
%REFERENCES:
%[1] S. G. Williamson, "Ranking algorithms for lists of partitions," SIAM
%    Journal on Computing, vol. 5, no. 4, pp. 602-617, 1976.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(q);

if(nargin<3||isempty(WT))
    %Compute the Williamson triangle
    WT=zeros(n,n);
    WT(1,1:n)=1;
    for v=2:n
        for mu=1:(n-v+1)
            WT(v,mu)=mu*WT(v-1,mu)+WT(v-1,mu+1);
        end
    end
end

q=q-1;

rank=0;

m=cummax(q);

for t=0:(n-2)
    idx=n-t;
    
    rank=rank+q(idx)*WT(t+1,m(idx-1)+1);
end

rank=rank+1;

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
