function [Q,P]=polarDecomp(A)
%%POLARDECOMP Compute the polar decomposition of the square matrix A. This
%             means find Q and P such that A=Q*P where Q is orthogonal and
%             P is symmetric and positive semidefinite.
%
%INPUTS: A An nXn matrix.
%
%OUTPUTS: Q An nXn orthogonal matrix.
%         P An nXn positive (semi)definite matrix.
%
%The technique for computing a polar decomposition in term of a singular
%value decomposition is given in Chapters 6.4.2 and 9.4.3 of [1].
%
%EXAMPLE:
% A=randn(6,6);
% [Q,P]=polarDecomp(A);
% max(abs(vec(A-Q*P)))
%One will see that the maximum error is on the order of finite precision
%limits.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[U,S,V]=svd(A);

Q=U*V';
P=V*S*V';

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
