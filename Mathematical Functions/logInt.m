function val=logInt(x)
%%LOGINT Evaluate the logarithmic integral of a real, non-negative x. The
%        logarthmic integral of x is equal to the integral from 0 to x of
%        1/log(t) dt.
%        
%INPUTS:    x A matrix of real values, at which one wishes to evaluate the
%             logarithmic integral. x does not have to be finite.
%
%OUTPUTS: val A matrix of the real values of the logarithmic integral
%             evaluated at the points in x. Note that val=-Inf if x=1.
%
%As noted in [1], the logarithmic integral of x can be written in terms of
%the exponential integral and natural logarithm functions. The exponential
%integral is a built-in function in Matlab, so this simply uses the expint
%function.
%
%%If the result obtained using the identity and the expint function is not
%finite, then x was too large, or too close to zero and asymptotic values
%of the logarithmic integral are substituted. If x was too large
%in magnitude then Inf is substituted. If x was too close to 0, then 0 is
%substituted.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Logarithmic Integral." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/LogarithmicIntegral.html
%
%June 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The real command deals with precision issues.
    val=real(-expint(-log(x))-1j*pi);

    %For values that are not finite, substitute the asymptotic values of
    %the logarithmic integral.
    %If x was too close to zero. 1 is an arbitrary choice.
    sel=~isfinite(val)&x<1;
    val(sel)=0;

    %If the magnitude of x was too large.
    sel=~isfinite(val)&x>1;
    val(sel)=Inf;
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
