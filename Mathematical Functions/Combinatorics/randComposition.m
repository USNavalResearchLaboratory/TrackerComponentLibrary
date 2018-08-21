function q=randComposition(t,n)
%%RANDCOMPOSITION Generate a random composition of n unlableled items and t
%                 slots. That is, a random assignment of putting n
%                 unlabeled balls into t labeled slots. This will not work
%                 if the total number of possible compositions
%                 (binomial(n-1,d-1)) is so large as to overflow.
%
%INPUTS: t The number of slots that can hold items.
%        n The number of items that are composed into slots.
%
%OUTPUTS: q An mX1 vector holding the random composition, whose elements
%           sum to n. Each element is the number of "balls" in that slot.
%
%This function just generates a random rank and unranks that composition
%using the unrankComposition function.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

totalCompositions=binomial(n-1,t-1);

if(totalCompositions==0)
    q=[];
    return;
end

%The min is for the (presumably zero probability) case that the random
%variable is 1.
rank=min(fix(rand(1)*totalCompositions),totalCompositions-1);

q=unrankComposition(rank,t,n);
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
