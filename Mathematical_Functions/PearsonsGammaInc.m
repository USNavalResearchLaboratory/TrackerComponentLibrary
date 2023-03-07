function val=PearsonsGammaInc(u,p)
%%PEARSONSGAMMAINC Evaluate Pearson's incomplete gamma function. This is
%          equal to (1/gamma(p+1))*integral_0^(u*sqrt(p+1)) exp(-t)*t^p dt.
%
%INPUTS: u The positive u value (or matrix for multiple values to be
%          computed in parallel) in the upper bound of the integral. u>=0.
%          It is allowable to pass a scalar for this and a matrix for p.
%        p The exponent of the t term in the integral, or a matrix for
%          multiple values to be done in parallel. p>-1. It is allowable to
%          pass a scalar for this and a matrix for u. 
%
%OUTPUTS: val The value or values of Person's gamma function at the desired
%             points.
%
%Pearson's incomplete gamma function is defined in terms of an integral in
%Section 6.5.6 of [1]. In Equation 1.7 of [2], Pearson's incomplete gamma
%function is expressed in terms of Tricomi's incomplete gamma function and
%Equation 1.3 expressed tricomi's gamma function in terms of the type of
%agmma function used in gammainc. Thus, this function uses the
%transformations and the gammainc function to express the result.
%
%REFERENCE
%[1] Abramowitz, M. and Stegun, I. A. (Eds.). "Gamma function and related
%    functions." in Ch. 6 in Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables, 9th printing. New York:
%    Dover, 1972.
%[2] W. Gautschi, "A computational procedure for incomplete gamma
%    functions," ACM Transactions on Mathematical Software, vol. 5, no. 4,
%    pp. 466-481, Dec. 1979.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isscalar(u)&&~isscalar(p))
    u=u*ones(size(p));
elseif(~isscalar(u)&&isscalar(p))
    p=p*ones(size(u));
end

val=gammainc(sqrt(1+p).*u,1+p);

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
