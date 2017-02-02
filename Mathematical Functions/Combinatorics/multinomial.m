function val=multinomial(repList)
%%MULTINOMIAL Compute a multinomial coefficient. The multinomial
%           coefficient is the number of permutation of a multiset. A
%           multiset is a set with repeated elements. Multinomial
%           coefficients are usually written like binomial coefficients but
%           with multiple terms below.
%
%INPUTS: repList A numUniqueX1 list of how many times each unique item in
%                the multiset is repeated. All elements >=1 and
%                n=sum(repList), where n is the total number of elements in
%                the multiset.
%
%OUTPUTS: val The value of the multinomial coefficient.
%
%The multinomial coefficient is factorial(n)/prod(factorial(repList)).
%However, that formulation is subject to overflow issues. Thus, it is
%implemented here as a product of boinomial coefficients.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numUniqueEls=length(repList);
    val=1;
    cumSum=0;
    for curEl=1:numUniqueEls
       cumSum=cumSum+repList(curEl);
       val=val*binomial(cumSum,repList(curEl));
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
