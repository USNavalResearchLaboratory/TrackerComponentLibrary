function numPerms=numMultisetPermutations(E)
%%NUMMULTISETPERMUTATIONS Determine the number of unique permutations of
%               the elements in a given multiset. A multiset is a set of
%               elements where some elements are repeated.
%
%INPUTS: E An nX1 or 1Xn vector of the elements in the multiset, including
%          repetitions. The elements of E must be finite; E cannot contain
%          any NaNs.
%
%OUTPUTS: numPerms The number of multiset permutations of E that are
%                  possible.
%
%The number of multiset permutations is given by the multinomial of a
%vector of the number of times each unique element is repeated.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

E=sort(E,'descend');
E=E(:);
%First, we will determine the number of times each unique element in the
%sorted E is repeated. This assumes that the values in E are all
%finite. 
numReps=diff(find(diff([Inf;E;-Inf])));
%The inner diff([Inf;E;-Inf]) marks where repetitions end. The outer
%find finds the indices of the nonzero elements. The differences between
%those indices is the number of elements in each repetitions. This is
%efficient as it avoids loops in Matlab, but it is an inefficient method of
%determining the multiplicity of the unique elements in E.

numPerms=multinomial(numReps);

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
