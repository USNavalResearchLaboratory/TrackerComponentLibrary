function perm=unrankPermutation(rank,n)
%%UNRANKPERMUTATION Return the permutation of the given rank in the
%                   lexicographic ordering of permutations consisting of
%                   numElements elements, where the where the first element
%                   of perm is the most significant. This requires that
%                   factorial(numElements) is small enough to be evaluated
%                   without a loss of precision.
%
%INPUTS: rank The order of the desired permutation of [1;2;3;...;n] in
%             lexicographic order. Note that 0<=rank<=(n!-1).
%           n The number of elements in the desired permutation.
%
%OUTPUTS: perm The permutation having the given lexicographic rank and
%              number of elements (having values 1 to n). If a rank equal
%              to or greater than the  total number of unique permutations
%              is given, then an empty matrix is returned.
%
%The algorithm is from [1].
%
%REFERENCES:
%[1] J. Liebehenschel, "Ranking and unranking of lexicographically ordered
%    words: An average-case analysis," Journal of Automata, Languages and
%    Combinatorics, vol. 26, no. 4, pp. 227-268, Jul. 2004.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(rank>=factorial(n))
        perm=[];
        return;
    end
    
    perm=zeros(n,1);%Allocate space for the result.
    a=1:n;
    
    for ii=0:(n-1)
       val=factorial(n-ii-1);
       k=floor(rank/val);
       rank=mod(rank,val);
       jj=0;
       while(k>=0)
           if(a(jj+1)~=0)
              k=k-1; 
           end
           jj=jj+1;
       end
       perm(ii+1)=a(jj-1+1)-1;
       a(jj-1+1)=0;
    end
    
    %The algorithm returns values from 0 to n-1. This changes it to 1 to n.
    perm=perm+1;
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
