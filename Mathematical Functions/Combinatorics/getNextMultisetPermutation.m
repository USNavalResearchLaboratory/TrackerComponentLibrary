function a=getNextMultisetPermutation(a,isInit)
%%GETNEXTMULTISETPERMUTATION Get the next permutation in a multiset in
%                   lexicographic order. If the input sequence has no
%                   repeats, then one is considering normal permutations,
%                   not multiset ones.
%
%INPUTS: a The multiset for which the next multiset permutation in the
%          sequence is desired. If this is just the multiset and one wants
%          to get the first permutation in the sequence, then isInit should
%          be true. Otherwise, false.
%   isInit This indicates whether a is an unordered multiset and one wishes
%          to get the first multiset permutation in lexicographic order.
%          The default if this parameter is omitted or an empty matrix is
%          passed is false, which means that this function will produce the
%          next permutation.
%
%OUTPUTS: a The nX1 vector that is next multiset permutation in
%           lexicographic order, unless isInit=true, in which case this is
%           the first value. If the input a was the last value in the
%           sequence and isInit=false, then this function returns an empty
%           matrix.
%
%This function implements Algorithm L of Chapter 7.2.1 of [1]. 
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(isInit))
   isInit=false; 
end

%If this is the first input in the sequence
if(isInit==true)
    a=sort(a(:),'ascend');
    return;
end

n=length(a);

%If a scalar is passed, then it is the only permutation.
if(n==1)
    a=[];
    return;
end

%Step L2, find j.
j=n-1;
while(a(j)>=a(j+1))
    j=j-1;
    %If the algorithm has gone through all multiset permutations.
    if(j==0)
        a=[];
        return;
    end
end

%Step L3
l=n;
while(a(j)>=a(l))
    l=l-1;
end
temp=a(j);
a(j)=a(l);
a(l)=temp;

%Step L4
a((j+1):n)=a(n:-1:(j+1));

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
