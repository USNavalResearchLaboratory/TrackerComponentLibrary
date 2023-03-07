function y=hav(x)
%%HAV Compute the haversine of x. The haversine arises in navigation
%     problems. As in [1], it is just sin(x/2).^2. Haversines tend to arise
%     in navigation problems.
%
%INPUTS: x A numeric matrix.
%
%OUTPUTS: y The haversine of x. y has the same dimensionality as x.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Haversine." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/Haversine.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

y=sin(x/2).^2;

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
