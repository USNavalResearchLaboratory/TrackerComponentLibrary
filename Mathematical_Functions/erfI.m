function val=erfI(z)
%%ERFI Evaluate the imaginary error function (often referred to as erfi).
%      As given in [1], this is simply -1j*erf(1j*z).
%
%INPUTS: z A scalar or matrix of real or complex values.
%
%OUTPUTS: val The value of the imaginary error function evaluated at the
%             points in z. This has the same dimensions as z.
%
%We must use erfComplex to evaluate the error function, because Matlab's
%built-in error function (as of R2022b) cannot handle complex terms.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Erfi." From MathWorld--A Wolfram Web Resource.
%    https://mathworld.wolfram.com/Erfi.html
%
%August 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=-1j*erfComplex(1j*z);

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
