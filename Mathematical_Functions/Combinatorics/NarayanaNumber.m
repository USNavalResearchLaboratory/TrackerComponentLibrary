function val=NarayanaNumber(n,k)
%%NARAYANANUMBER Generate the Narayana number specified by n and k. This is
%           the number of Dyck paths of length n with k peaks (local
%           maxima). Narayana numbers are defined, for example, in [1] and
%           their relation to Dyck paths is given in [2].
%           
%INPUTS: n The path length.
%        k The number of peaks.
%
%OUTPUTS: val The specified Narayana number.
%
%EXAMPLE:
%Here, we confirm that the sum from 1 to k of the Narayana numbers equals
%the nth Catalan number.
% n=20;
% sumVal=0;
% for k=1:n
%     sumVal=sumVal+NarayanaNumber(n,k);
% end
% RelErr=(sumVal-CatalanNumber(n))/CatalanNumber(n)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Narayana Number." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/NarayanaNumber.html
%[2] D. Merlino, R. Sprugnoli, and M. C. Verri, "Some statistics on Dyck
%    paths," Journal of Statistical Planning and Inference, vol. 101, no.
%    1-2, pp. 211-227, 15 Feb. 2002.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=binomial(n,k)*binomial(n,k-1)/n;

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
