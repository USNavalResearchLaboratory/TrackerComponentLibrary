function val=erfcxComplex(z)
%%ERFCXCOMPLEX Evaluate the scaled complementary error function in a manner
%             valid for complex arguments. The complementary error function
%             is one minus the error function (erf). The scaled
%             complementary error function is exp(z^2)*erfc(z). The
%             functions erfcx and erfc that are built into Matlab cannot
%             handle complex arguments.
%
%INPUTS: z A scalar, vector, or matrix of values at which one wishes to
%          evaluate the scaled complementary error function.
%
%OUTPUTS: val A matrix having the same dimensions as z in which the values
%            of the scaled complementary error function are computed for
%            the values in z.
%
%The scaled complementary error function is just Faddeeva(1i*z) based on
%the definition of the Faddeeva function.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=Faddeeva(1i*z);

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
