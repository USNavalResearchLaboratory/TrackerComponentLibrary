function val=numArrangements(n,m)
%%NUMARRANGEMENTS Get the total number of way of putting n items into m
%                 spaces. That is the total number of arrangements of n
%                 items into m parts. If m>n, then 0 is returned.
%
%INPUTS: n The total number of items from which one draws.
%        m The number of slots into which these items can be placed.
%
%OUTPUTS: val The number of arrangements.
%
%Arragements are generating by going through all possible combinations of
%which of the n elements are in the m slots and then permuting the ordering
%of the elements into the slots. Thus, the total number of possibilities
%the the product of the number of combinations and the number of
%permutation of what was chosen in each combination.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=binomial(n,m)*factorial(m);

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
