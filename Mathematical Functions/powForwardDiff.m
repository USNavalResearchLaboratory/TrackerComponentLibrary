function val=powForwardDiff(k,w,epsVal)
%%POWFORWARDDIFF This function evaluate the forward difference for two
%         power functions: (w^(k+epsVal)-w^k)/epsVal . This in done in
%         such a manner that the result is numerically stable for
%         epsVal>=0, whereas simply evaluating (w^(k+epsVal)-w^k)/epsVal
%         would fail as epsVal approaches zero whereas this approaches the
%         asymptotic limit of w^k*log(w).
%
%INPUTS:     k The real exponent without the offset. -Inf<k<Inf;
%            w The real base that is being exponentiated, w>0.
%       epsVal The offset amount. epsVal>=0. If epsVal<0 is passed, it is
%              substituted with abs(epsVal);
%
%OUTPUTS: val The value (w^(k+epsVal)-w^k)/epsVal.
%
%This function does a check. If abs(epsVal*log(w)) is greater than 1/8,
%then the function is evaluated directly, Otherwise, a Taylor series
%expansion is used.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Arbitrary: If it doesn't converge in 500 terms, it probably is not going
%to.
maxNumTerms=5000;

test=epsVal*log(w);

if(abs(test)>=1/8)
    temp=(exp(test)-1)/epsVal;
    val=w^k*temp;
else
    %This is a Taylor series expansion of (exp(test)-1)/epsVal.
    temp=-log(w);
    i=1;
    didSucceed=false;
    for curTerm=1:maxNumTerms
        nextTerm=((-epsVal)^i*(-log(w))^(i+1))/gamma(i+2);
        if(abs(nextTerm)<eps())
            didSucceed=true;
            break;
        end
        if(~isfinite(temp))
            break;
        end
        temp=temp+nextTerm;
        i=i+1;
    end
    if(didSucceed==false)
        warning('Convergence not obtained.')
    end
    
    val=-w^k*temp;
end
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
