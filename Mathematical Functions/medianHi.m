function medHi=medianHi(x)
%%MEDIANHI Find the high (upper) median of a set of real values. That is if
%          the n values were sorted, then return the value in index
%          fix(n/2)+1.
%
%INPUTS: x An nX1 or 1Xn vector. The elements do not need to be in any
%          particular order.
%
%OUTPUTS: medHi The value in x that holds the high median.
%
%The function kthOrderStat is used to obtain the value without sorting.
%
%EXAMPLE:
%We look at what one obtains with an odd and an even number of values.
% medHi10=medianHi(1:10)
% medHi11=medianHi(1:11)
%One gets medHi10=6 and medHi11=6.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);

k=fix(n/2)+1;
medHi=kthOrderStat(x,k);

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
