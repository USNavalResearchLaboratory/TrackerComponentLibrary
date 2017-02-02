function realRoot=realRootCubEq(a,b,c)
%%REALROOTCUBEQ  All cubic equations with real coefficients have at least
%                one real root. This function finds the real root of a
%                cubic equation of the form x^3+a*x^2+b*x+c=0.
%
%INPUTS: a,b,c  The real values of the coefficients in a cubic equation of
%               the form x^3+a*x^2+b*x+c=0.
%
%OUTPUTS: realRoot The real root of the cubic equation. If complex roots
%                  exist, this is unique. If multiple real roots exist,
%                  this is just one of the values.
%
%The solution comes from the general cubic formula given in [1] and noting
%that if the arguments of the square roots are negative, the complex parts
%will cancel in the subsequent equation.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Cubic Formula." From MathWorld --A Wolfram Web
%    Resource. October 2014. http://mathworld.wolfram.com/CubicFormula.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 22 on Mathworld.
Q=(3*b-a^2)/9;
%Equation 23 on Mathworld.
R=(9*a*b-27*c-2*a^3)/54;

%When summing S and T, the complex parts will cancel. This gets rid of any
%complex parts to begin with.
S=nthroot(R+real(sqrt(Q^3+R^2)),3);
T=nthroot(R-real(sqrt(Q^3+R^2)),3);

%Equation 52 on Mathworld.
realRoot=S+T-a/3;
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
