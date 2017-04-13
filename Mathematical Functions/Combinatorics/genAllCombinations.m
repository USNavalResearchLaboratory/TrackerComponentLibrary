function theCombos=genAllCombinations(n,k)
%%GENALLCOMBOS Generate all combinations of values from 0 to (n-1) in
%              lexicographic order. there are binomial(n,k) possible
%              combinations.
%
%INPUTS: n The number of items from which k items are chosen for the ranked
%          combinations.
%        k The number of items chosen.
%
%OUTPUTS: theCombos A kXnumCombos matrix containing all possible
%                   combinations of values.
%
%This function just calls the getNextCombo function in a loop.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numCombos=binomial(n,k);

theCombos=zeros(k,numCombos);

theCombos(:,1)=0:(k-1);
for curCombo=2:numCombos
    theCombos(:,curCombo)=getNextCombo(theCombos(:,curCombo-1),n);
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
