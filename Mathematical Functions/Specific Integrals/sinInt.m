function val=sinInt(x)
%%SININT Evaluate the sine integral of a real value x. The sine integral of
%        x is defined to be the integral from 0 to x of sin(t)/t dt. The
%        derivative of the sine integral function is the sinc function (as
%        defined without pi in the numerator and denominator); that is
%        sin(x)/x.
%
%INPUTS: x A matrix of real values, at which one wishes to evaluate the
%          sine integral. x does not have to be finite.
%
%OUTPUTS: val A matrix of the real values of the sine integral evaluated at
%             the points in x.
%
%The sine integral can be related to the exponential integral as discussed
%in [1]. This function uses Matlab's built-in expint function to evaluate
%the sine integral using the identity relating the two when x is finite and
%negative. Due to the odd symmetry of the sine integral, the same identity
%can be used when x is positive by flipping the sign of x and the sign of
%the result.
%
%If the result obtained using the identity and the expint function is not
%finite, then x was too large, too small or too close to zero and
%asymptotic values of the sine integral are substituted. If x was too large
%in magnitude the -/+ pi/2 is substituted for small/ large x. If x was too
%close to 0, then 0 is substituted.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Sine Integral." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/SineIntegral.html
%
%June 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The real function is used to deal with any precision problems that
    %might cause the value of the sine integral to be imaginary.
    signX=sign(x);
    x=abs(x);
    val=signX.*(real(1/(2*1i)*(expint(1i*x)-expint(-1i*x)))+pi/2);

    %For values that are not finite, substitute the asymptotic values of
    %the sine integral.
    %If x was too close to zero. 1 is an arbitrary choice.
    sel=~isfinite(val)&x<1;
    val(sel)=0;

    %If the magnitude of x was too large.
    sel=~isfinite(val)&x>1;
    val(sel)=signX(sel)*pi/2;
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
