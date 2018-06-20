function [R,WT]=unrankSetPartition(n,rank,WT)
%%UNRANKSETPARTITION Return the set partition of n items of a given rank in
%                    lexicographic order.
%
%INPUTS: n The number of elements in the set partition n>=1.
%     rank The rank of the item in the set partition desired.
%          0<=rank<BellNumber(n) 
%       WT The is, optionally, the Williamson table. It is the table T in
%          [1]. The table is a function of length(q). If omitted, then the
%          table will be computed and returned so that futue calls to the
%          function can use 
%
%OUTPUTS: q The nX1 set partition vector corresponding to rank. Paritions
%           are numbered in increasing order starting from 1. That is, the
%           first occurence of partition numbered i in q is always before
%           the first occurences of partitions numbered j>1.
%        WT The Williamson table for a length-n set partition. Passing this
%           to rankSetPartition on subsequent calls will speed up the
%           algorithm.
%
%This function implements the unranking algorithm of [1]. The algorithm is
%described in a slightly more direct manner in [2], along with a parallel
%variant, which is not implemented here.
%
%REFERENCES:
%[1] S. G. Williamson, "Ranking algorithms for lists of partitions," SIAM
%    Journal on Computing, vol. 5, no. 4, pp. 602-617, 1976.
%[2] Z.Kokosinski,"A parallel dynamic programming algorithm for unranking
%    set partitions," Technical Transactions Automatic Control, vol. 11,
%    2013.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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

R=ones(n,1);

i=n-1;
j=1;
m=j;

while(rank>0)
    if(m*WT(i,j)<=rank)
        rank=rank-m*WT(i,j);
        R(n-i+1)=R(n-i+1)+m;
        j=max(j,R(n-i+1));
        m=j;
        i=i-1;
    else
        m=m-1;
    end
end


end