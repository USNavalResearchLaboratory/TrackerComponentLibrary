function y=expSlope1(x)
%%EXPSLOPE1  Evaluate the function (exp(x)-1)/x avoiding problems at the
%            point x=0, where the limit equals 1. This function arises in
%            range enclosure algorithms for interval arithmetic.
%
%INPUTS: x A vector of matrix of real values at which the function should
%          be evaluated.
%
%OUTPUTS: y The value(s) of the function (exp(x)-1)/x at the point(s) in x.
%
%The function is implemented as in Chapter 1.14.1 of [2], which only 
%
%This is one of the slope function in Table 10.6 of Section 10.6.2 of [1].
%This file implements the function for non-intervals. See the Interval
%class for an implementation for Interval arithmetic.
%
%REFERENCES:
%[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
%    pp.1-97, June 30 2015.
%[2] N. J. Higham, Accuracy and Stability of Numerical Algorithms.
%    Philadelphia: SIAM, 1996.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate space
y=exp(x);
sel=(y==1);

y(sel)=1;
y(~sel)=(y-1)/log(y);

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
