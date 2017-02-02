function y=sinSlope3(x)
%%SINSLOPE3  Evaluate the function -6*(sin(x)-x)/x^3 avoiding
%            problems at the point x=0, where the limit equals 1. This
%            function arises in range enclosure algorithms for interval
%            arithmetic.
%
%INPUTS: x A vector of matrix of real values at which the function should
%          be evaluated.
%
%OUTPUTS: y The value(s) of the function -6*(sin(x)-x)/x^3 at the
%           point(s) in x.
%
%For values of x having magnitude over 1e-2 or for complex values, the
%function is evaluated as written. For smaller real values, a Taylor series
%expansion is used to avoid the singularity near the origin.
%
%This is two times one of the slope function in Table 10.6 of Section
%10.6.2 of [1]. As the standard says the value at 0 is 1, but the value of
%the function in the standard at 0 is 1/2, it is assumed that there is a
%typo in the standard. Hence this function is twice that listed in [1].
%This file implements the function for non-intervals. See the Interval
%class for an implementation for Interval arithmetic.
%
%REFERENCES:
%[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
%    pp.1-97, June 30 2015.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate space
y=zeros(size(x));

%For large enough values, use the standard function.
sel=abs(x)>1e-2|~isreal(x);
x1=x(sel);
y(sel)=-6*(sin(x1)-x1)./x1.^3;

%For smaller magnitude values, use an sixth-order Taylor series expansion.
%The relative error at 10^-2 is on the order of 1e-23, which is less than
%the precision of a double floating point number at that value. The series
%is given in Horner form.
x1=x(~sel);
y(~sel)=1+x1.^2.*(-(1/20)+x1.^2.*(1/840-x1.^2/60480));

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
