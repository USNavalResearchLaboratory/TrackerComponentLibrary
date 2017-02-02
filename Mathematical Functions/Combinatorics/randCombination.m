function combo=randCombination(n,m)
%%RANDCOMBINATION  Generate a random combination of m elements drawn from a
%            set of n elements. This will not work if the total number of
%            possible combinations (binomial(n,m)) is so large as to
%            overflow.
%
%INPUTS:    n       The number of items from which m items are chosen for
%                   the combination.
%           m       The number of items chosen.
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order. The lowest item is indexed zero.
%
%The function just chooses a random combination
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

totalCombo=binomial(n,m);
%The min is for the (presumably zero probability) case that the random
%variable is 1.
rank=min(fix(rand(1)*totalCombo),totalCombo-1);

combo=unrankCombination(rank,n,m);
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
