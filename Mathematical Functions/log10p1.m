function y=log10p1(x)
%%LOG10P1 Evaluate the function log10(x+1) without a loss of precision for
%         values of x near 0.
%
%INPUTS: x A vector of matrix of real values at which the function should
%          be evaluated.
%
%OUTPUTS: y The value(s) of the function log10(x+1) at the point(s) in x.
%
%For values of x having magnitude over 1e-2 and for complex values, the
%function is evaluated as written. For smaller values, a Taylor series
%expansion is used to avoid the singularity near the origin.
%
%EXAMPLE:
%Compare
% log10(1e-16+1)
% log10p1(1e-16)
%For the first one, one will get 0, for the second,
%4.342944819032518e-17. Compare this to a value computed using
%higher-precision arithmetic of 4.3429448190325180593640482375402*10^-17
%and it is clear that this function is more accurate than explicitly
%evaluating the expression.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate space
y=zeros(size(x));

%For large enough values, use the standard function.
sel=abs(x)>1e-2|~isreal(x);
x1=x(sel);
y(sel)=log10(x1+1);

%For smaller magnitude values, use an ninth-order Taylor series expansion.
%The relative error at 10^-2 is on the order of 1e-20, which is less than
%the precision of a double floating point number at that value. The series
%is given in Horner form.
log10Val=log(10);
x1=x(~sel);
y(~sel)=x1.*(x1.*(x1.*(x1.*(x1.*(x1.*(x1.*(x1.*(-(1/(8*log10Val))+x/(9*log10Val))+1/(7*log10Val))-1/(6*log10Val))+1/(5*log10Val))-1/(4*log10Val))+1/(3*log10Val))-1/(2*log10Val))+1/log10Val);
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
