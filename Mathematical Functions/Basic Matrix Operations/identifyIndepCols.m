function colIdx=identifyIndepCols(A,tol)
%%IDENTIFYINDEPCOLS  Given a matrix A, which can be square or rectangular,
%                    find a linearly independent subset of the columns.
%                    The total number of linearly independent columns is
%                    determined by the rank function.
%
%INPUTS: A A matrix whereby a linearly independent subset of columns is
%          desired.
%      tol An optional parameter specifying the tolerance for determinaing
%          the rank of the matrix. If omitted, the value
%          tol=max(size(A))*eps(norm(A)) is used.
%
%OUTPUTS: colIdx A vector of the indices of linearly independent columns of
%                A. The indices are ordered in terms of increasing
%                magnitude of the diagonal of R is a QR decomposition.
%
%The use of a QR decomposition with column pivoting for finding linearly
%independent columns is described in Chapter 5.5 of [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    tol=max(size(A))*eps(norm(A));
end

numCol=rank(A,tol);

[~,~,permVec] = qr(A,'vector');
colIdx=permVec(1:numCol);
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
