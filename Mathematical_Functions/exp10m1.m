function y=exp10m1(x)
%%EXP10M1 Evaluate the function 10^x-1 without a loss of precision for 
%         values of x near 0.
%
%INPUTS: x A vector of matrix of real values at which the function should
%          be evaluated.
%
%OUTPUTS: y The value(s) of the function 10^x-1 at the point(s) in x.
%
%For values of x having magnitude over 1e-2, the function is evaluated as
%written. For smaller values, a Taylor series expansion is used to avoid
%the singularity near the origin.
%
%EXAMPLE:
%Compare
% 10^(1e-17)-1
% exp10m1(1e-17)
%For the first one, one will get 0, for the second,
%2.302585092994046e-17. Compare this to a value computed using
%higher-precision arithmetic of 2.3025850929940457105274820070764*10^-17
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
y(sel)=10.^x1-1;

%For smaller magnitude values, use an eigth-order Taylor series expansion.
%The relative error at 10^-2 is on the order of 1e-19, which is less than
%the precision of a double floating point number at that value. The series
%is given in Horner form.
log10=log(10);
x1=x(~sel);
y(~sel)=log10*x1.*(1+log10*x1.*(1/2+log10*x1.*(1/6+log10*x1.*(1/24+log10*x1.*(1/120+log10*x1.*(1/720+log10*x1.*(1/5040+(x1.*log10)/40320)))))));

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
