function F=logMat(A)
%%LOGMAT Evaluate the matrix logarithm base e of a particular matrix nXn A.
%       This is most easily defined in terms of the Taylor series
%       expansion:
%       log(A)=sum_{k=1}^Inf (-1)^(k+1) (A-eye(n))^k/factorial(k)
%       The opposite of this function is expMat.
%
%INPUTS: A An nXn real or complex matrix.
%
%OUTPUTS: F The matrix logarithm (base e) of A.
%
%The scalar log function is applied to the matrix using the function
%SchurMatFunEval.
%
%EXAMPLE:
% X=randn(7,7);
% diff=logMat(expMat(X))-X
%One will see that all the values in diff are on the order of 1e-14
%indicating the there is agreement to finite precision limits.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    F=SchurMatFunEval(@(x)log(x),A);
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
