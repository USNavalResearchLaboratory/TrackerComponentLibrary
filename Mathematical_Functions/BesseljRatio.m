function val=BesseljRatio(nu,z)
%%BESSELJRATIO Evaluate the ratio of Bessel functions of the first
%              kind of the form x=J_{nu-1}(z)/J_{nu}(z). 
%
%INPUTS: nu The positive (upper) order of the Bessel function of the first
%           kind in the ratio. This does not need to be an integer.
%         z The complex argument of the Bessel function of the first kind;
%           z~=0.
%
%OUTPUTS: val The value of the ratio J_{nu-1}(z)/J_{nu}(z).
%
%The continued fraction representation of [1] is used with the function
%evalContFrac.
%
%REFERENCES:
%[1] W. J. Lentz, "Generating Bessel functions in Mie scattering
%    calculations using continued fractions," Applied optics, vol. 15, no.
%    3, pp. 668-671, Mar. 1976.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

a=@(m)1;
b=@(m)(-1)^(m+1)*2*(nu+m-1)/z;

val=evalContFrac(b,a);

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
