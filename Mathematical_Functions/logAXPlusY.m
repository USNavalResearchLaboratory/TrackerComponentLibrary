function val=logAXPlusY(logA,logX,logY)
%%LOGAXPLUSY  Find the element-by-element natural logarithm of a.*x+y, for
%             non-negative, vector or scalar a,x, and y, given the
%             logarithms of a, x, and y, trying to avoid overflow and
%             underflow errors in the product and sum. This lets one
%             evaluate log(a*x+y) even if e^(log(A), e^(log(X)) and
%             e^(log(Y)) all would overflow.
%
%INPUTS: logA, logX, logY The real natural logarithms of the non-zero
%                         a, x, and y that one want in the evaluation of
%                         log(a.*x+y). They can be scalars or vectors/
%                         matrices. The logarithms are element-by element,
%                         not matrix.
%
%OUTPUTS: val The value log(a.*x+y) for each corresponding element in a,x,
%             and y.
%
%The algorithm just uses some identities of logarithms and exponentials to
%minimize overflow and underflow problems. Large enough logarithmic values
%can still eventually have underflow/ overflow problems, but at values
%significantly larger than when a, x and y individually would overflow. The
%algorithm is still subject to the same loss of precision one would expect
%to see when evaluating log(x+1) and x is close to zero (e.g x=eps(1)).
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

logm=min(logA+logX,logY);
logM=max(logA+logX,logY);

%Underflow problems can occur cause 1+exp(logm-logM)=1. However, this seems
%inevitable under any method chosen.
val=logM+log(1+exp(logm-logM));
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
