function val=cosInt(x)
%%COSINT Evaluate the cosine integral of a real non-negative value x. The
%        cosine integral of x is defined to be the negative of the integral
%        from x to infinity of cos(t)/t dt. The derivative of the cosine
%        integral function is cos(x)/x.
%
%INPUTS:    x A matrix of non-negative real values, at which one wishes to
%             evaluate the cosine integral. x does not have to be finite.
%
%OUTPUTS: val A matrix of the real values of the cosine integral evaluated 
%             at the points in x.
%
%The cosine integral can be related to the exponential integral as
%discussed in [1]. This function uses Matlab's built-in expint function to
%evaluate the cosine integral using the expint function.
%
%The cosine integral asymptotically approaches zero. If the results from
%the exponential integral functions are not finite, then it is assumed that
%the corresponding values of x were too large and the asymptotic result (0)
%is substituted.
%
%x should be non-negative. If negative values of x are (incorrectly)
%provided, then the value of cosInt(abs(x)) is returned. Strict definitions
%of the cosine integral are usually such that negative values return
%complex numbers.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Cosine Integral." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/CosineIntegral.html
%
%June 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The real function is in case of finite precision problems.
    val=-0.5*real(expint(1j*x)+expint(-1j*x));
    val(~isfinite(val))=0;
    
    %Make sure that NaN inputs remain NaN on the output.
    val(isnan(x))=NaN;
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
