function rank=rankPermutation(perm)
%%UNRANKPERMUTATION Return the ranked position of a permutation in
%                   the lexicographic ordering, where the first element of
%                   perm is the most significant.
%
%INPUTS: perm An n-length vector consisting of integers from 1 to n.
%
%OUTPUTS: rank The rank of the permutation in a lexicographic ordering of
%              permutations, counting from zero.
%
%The algorithm is from [1].
%
%REFERENCES:
%[1] J. Liebehenschel, "Ranking and unranking of lexicographically ordered
%    words: An average-case analysis," Journal of Automata, Languages and
%    Combinatorics, vol. 26, no. 4 pp. 227-268, Jul. 2004.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(perm);

%The algorithm is meant for permutation indices starting from 0. This
%adjusts perm for that.
perm=perm-1;

rank=0;
for ii=0:(n-2)
    k=perm(ii+1); 
    for jj=0:(ii-1)
        if(perm(jj+1)<perm(ii+1))
           k=k-1;
        end
    end
    rank=rank+k*factorial(n-ii-1);
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
