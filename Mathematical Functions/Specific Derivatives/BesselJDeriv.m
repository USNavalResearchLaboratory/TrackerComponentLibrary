function val=BesselJDeriv(m,z)
%%BESSELJDERIV Evaluate the derivative of a Bessel function of the first
%              kind. That is, the derivative of besselj(m,z) with respect
%              to z.
%
%INPUTS: m The order of the function. This can be any real number.
%        z The real or complex argument of the function. Passing a matrix
%          will evaluate the derivative at all of the points in the matrix.
%
%OUTPUTS: val The value fo the derivative of the Bessel function of the
%             first kind. This has the same size as z.
%
%The derivative of a Bessel function can be expressed as the difference of
%two Bessel functions of differing orders. The formula used below can be
%verified in common symbolic differentiation programs such as Mathematica.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=(1/2)*(besselj(m-1,z)-besselj(m+1,z));

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
