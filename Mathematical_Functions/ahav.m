function y=ahav(x)
%%AHAV Compute the inverse haversine (the archaversine) of x . The
%      haversine of x is sin(x/2).^2, so the inverse haversine is
%      2*asin(sqrt(x)), as in [1]. Haversines tend to arise in navigation
%      problems.
%
%INPUTS: x A numeric matrix.
%
%OUTPUTS: y The archaversine of x. y has the same dimensionality as x.
%
%EXAMPLE:
%This shows that the inverse works (within the limits of the uniqueness of
%inputs to hav). The error is 0.
% y=1.1;
% Err=ahav(hav(y))-y
%
%REFERENCES:
%[1] Weisstein, Eric W. "Inverse Haversine." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/InverseHaversine.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

y=2*asin(sqrt(x));

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
