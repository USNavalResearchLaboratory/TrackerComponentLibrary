function val=sincFun(x)
%%SINCFUN Evaluate the sinc function. This is defined to be
%         sin(pi*x)/(pi*x). Note that some authors use the definition
%         sin(x)/x. This function accounts for the singularity at zero.
%
%INPUTS: x A matrix of the values at which to evaluate the sinc function.
%
%OUTPUTS: vals A matrix the same size as x where the sinc function has been
%              evaluted at each of the elements.
%
%The sinc function is described in [1]. It arises in many areas of signal
%processing.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Sinc Function." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/SincFunction.html
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    y=pi*x;
    sel=(y==0);
    val=sin(y)./y;
    val(sel)=1;
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
