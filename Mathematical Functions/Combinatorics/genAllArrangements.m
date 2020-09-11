function arrs=genAllArrangements(n,m)
%%GENALLARRANGEMENTS Return all possible arrangement of n items put into m
%                    spaces. The arrangements are not in lexicographic
%                    order.
%
%INPUTS: n The total number of items from which one draws. n>=0
%        m The number of slots into which these items can be placed. m>=0
%
%OUTPUTS: arrs An mX(binomial(n,m)*factorial(m)) matrix of all possible
%              arrangements of n items into m slots. If m>n or m=0, then an
%              empty matrix is returned.
%
%Arragements are generating by going through all possible combinations of
%which of the n elements are in the m slots and then permuting the ordering
%of the elements into the slots. The getNextCombo function is used to
%generate the combinations and the perms function is used to generate the
%permutations.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(m>n||m==0)
   arrs=[];
   return;
end

numElements=binomial(n,m)*factorial(m);
%Allocate space
arrs=zeros(m,numElements);

numPerm=factorial(m);

curEl=1;
%Go through all combinations of which m elements are chosen.
curCombo=0:(m-1);
while(~isempty(curCombo))
    %For the given combo, go through all permutations of the elements.
    arrs(:,curEl:(curEl+numPerm-1))=perms(curCombo)';
    curEl=curEl+numPerm;
    curCombo=getNextCombo(curCombo,n);
end

arrs=arrs+1;
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
