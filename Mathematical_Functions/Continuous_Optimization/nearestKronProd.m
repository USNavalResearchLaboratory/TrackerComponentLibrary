function [B,C]=nearestKronProd(A,m1,m2,n1,n2)
%%NEARESTKRONPRODSVD Given an (m1*m2)X(n1*n2) matrix, determine the set of
%           matrices B and C that minimize norm(A-kron(B,C),'fro'). B is
%           m1Xn1 and C is m2Xn2. B and C are not unique. One can always
%           replace them with x*B and (1/x)*C for some arbitrary scalar x.
%
%INPUTS: A An (m1*m2)X(n1*n2) real matrix which is an m1Xn1 collection of
%          blocks of size m2Xn2. It has the form
%          A=[A_{1,1} ... A_{1,n1};
%             ...     ...      ...;
%             A_{m1,1} ... A_{m1,n1}];
% m1,m2,n1,n2 The dimensions specifying A and its submatrices as described
%          above. These specify the dimensions of B and C.
%
%OUTPUTS: B The m1Xn1 left matrix in the Kronecker product.
%         C The m2Xn2 right matrix in the Kronecker product.
%
%This function implements the algorithm described in Section 12.3.6 of [1].
%
%EXAMPLE:
%Here, we compute A explicitely from a Kronecker product of B and C
%matrices. The nearest Kronecker product obtained by this function should
%thus have zero error, which, within finite precision limits, is
%demonstrated.
% m1=2;
% m2=3;
% n1=5;
% n2=4;
% B=randn(m1,n1);
% C=randn(m2,n2);
% A=kron(B,C);
% 
% [BK,CK]=nearestKronProd(A,m1,m2,n1,n2);
% norm(A-kron(BK,CK),'fro')%The error will be close to 0.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[sigma,U,V]=KroneckerProdSVD(A,m1,n1,m2,n2);

sqrtSigma1=sqrt(sigma(1));

B=sqrtSigma1*U(:,:,1);
C=sqrtSigma1*V(:,:,1);

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
