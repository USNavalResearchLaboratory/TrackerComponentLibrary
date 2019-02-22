function val=subfactorial(n)
%%SUBFACTORIAL Compute the subfactorial of positive integers. Whereas the
%              factorial of a number n is the number of permutations of n
%              items, the subfactorial of n is the number of deragements of
%              n items. That is, the number of ways that n items can be
%              rearranged such that no items remain in their original
%              positions.
%
%INPUTS: n A scalar or matrix of positive, real integers whose
%          subfactorials are desired; n>=1.
%
%OUTPUTS: val The subfactorial of the values in n. This is often written as
%             !n as opposed to the notation n! for the factorial.
%
%As noted in [1], the subfactorial is just round(factorial(n)/exp(1)).
%here, we use round(exp(gammaln(n+1)-1)) instead to lessen the effects of
%overflow errors.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Subfactorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Subfactorial.html
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=round(exp(gammaln(n+1)-1));

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