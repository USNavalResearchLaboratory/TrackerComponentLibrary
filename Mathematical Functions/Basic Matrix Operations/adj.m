function AAdj=adj(A)
%%ADJ Compute the adjugate of a matrix (also known as the "classical
%     adjoint of a matrix, which is not to be confused with the use of
%     "adjoint" to mean the conjugate transpose of a matrix). The
%     adjugate of a matrix is defined in Ch. 4 of [1] as the transpose of
%     the cofactor matrix C. The cofactor matrix of a matrix A is one such
%     that C(i,j)=(-1)^(i+j)det(A_{i,j}), where A_{i,j} is the matrix A
%     with row i and column j removed.
%
%INPUTS: A An nXn square matrix. The elements can be real or complex. The
%          matrix can be singular.
%
%OUTPUTS: AAdj The nXn adjugate of the matrix A.
%
%In order to efficiently compute the matrix adjugate, a singular value
%decomposition (SVD) method is used. As noted in Chapter 4.3 of [1], 
%adj(A*B)=adj(B)*adj(A) (The adjugate operator is an antihomomorphism). On
%the other hand, the SVD of a matrix A is A=U*S*V', where S is diagonal,
%and U and V are unitary matrices (meaning that U is orthonormal and
%U'=inv(U)) (See Ch. 2.4.4 of [2]). Thus, we can write
%adj(A)=adj(U*S*V')=adj(V')*adj(S)*adj(U). For univary matrices, the
%adjugate is just the conjugate transpose of the matrix times the
%determinant of the matrix (which will be zero or 1). This can be proven
%from the identity that inv(U)=adj(U)/det(U) (in Chapter 4.3 of [1]) for
%any invertible matrix U, and that the determinant of a unitary
%matrix is +/-1. Thus, two of the three terms in adj(V')*adj(S)*adj(U) are
%easy to compute. Finally, the adjugate of a diagonal matrix is such that
%each non-diagonal entry is zero and the diagonal entries in row/ column i
%is the product of all of the diagonal elements excluding row/column i.
%
%The matrix adjugate is also equal to adj(A)=det(A)*inv(A) if the matrix A
%is invertible. This identity is not used here.
%
%REFERENCES:
%[1] K. M. Abadir and J. R. Magnus, "Matrix algebra," New York, 2005.
%[2] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[U,S,V]=svd(A);
s=diag(S);%Just keep the diagonal elements (the rest are zeros).

%To hold the diagonal elements of the adjoint of S.
nS=length(s);
sAdj=zeros(nS,1);
for k=1:nS
   sAdj(k)=prod(s([1:(k-1),(k+1):nS])); 
end

AAdj=(det(V')*V)*diag(sAdj)*(det(U)*U');

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
