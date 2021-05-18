function K=commutationMatrix(m,n)
%%COMMUTATIONMATRIX Return a commutation matrix K for an mXn matrix. This
%          means that given an mXn matrix A that K*A(:)=B(:) where
%          B=transpose(A). As shown below, this matrix can also be used for
%          reversing the order of Kronecker products. As in [2], this
%          matrix is also known as a vec-permutation matrix.
%
%INPUTS: m,n The number of rows and the number of columns of the matrix
%            that is to be commutated by this commutation matrix.
%
%OUTPUTS: K An (m*n)X(m*n) commutation matrix.
%
%The matrix K is just a type of permutation matrix with 1's in the correct
%places.
%
%Properties of the commutation matrix are given in [1] and in [2]. For
%example, commutationMatrix(m,n)=commutationMatrix(n,m)' and also there is
%an inverse relation that
%commutationMatrix(m,n)*commutationMatrix(n,m)==eye(n*m) and
%Note that the Vec permutation matrix I_{m,n} as in [2] is actually
%commutationMatrix(n,m).
%
%EXAMPLE 1:
%Here, we demonstrate the property that the vec-permutation matrix
%(commutation matrix) allows one to get vec(A') from vec(A):
% m=12;
% n=17;
% A=randn(m,n);
% Imn=commutationMatrix(n,m);%The vec-permutation matrix as in [2].
% all(vec(A)==Imn*vec(A'))
%The result will be true, showing that
%vec(A)==commutationMatrix(n,m)*vec(A')
%which agrees with the identity given in the abstract to [2]. 
%
%EXAMPLE 2:
%Here, we give a demonstration of how to reverse the ordering of a
%Kronecker product of two matrices:
% m=12;
% n=17;
% A=randn(m,n);
% p=4;
% q=6;
% B=randn(p,q);
% Imp=commutationMatrix(p,m);
% Iqn=commutationMatrix(n,q);
% %Reversing the order of the matrix Kronecker product as in Equation 25 in
% %[2].
% all(vec(kron(B,A)==Imp*kron(A,B)*Iqn))
%The result will be true, showing that one can use two commutation matrices
%to reverse the ordering of Kronecker products of matrices.
%
%EXAMPLE 3:
%Now, we demonstrate that one can use a single commutation matrix (vec-
%permutation matrix) to reverse the order of the Kronecker product of two
%vectors.
% m=13;
% n=7;
% a=randn(m,1);
% b=randn(n,1);
% Imat=commutationMatrix(m,n);
% all(vec(kron(a,b)==Imat*kron(b,a)))
%The result will be true.
%
%REFERENCES:
%[1] J. R. Magnus and H. Neudecker, "The elimination matrix: Some lemmas
%    and applications," SIAM Journal on Algebraic Discrete Methods, vol. 1,
%    no. 4, pp. 422-449, 1980.
%[2] H. V. Henderson and S. R. Searle, "The vec-permutation matrix, the vec
%    operator and Kronecker products: A review," Linear and Multilinear
%    Algebra, vol. 9, no. 4, pp. 271-288, 1981.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

K=zeros(m*n,m*n);

idx=0;
col=0;

for i=1:n
    idx=idx+1;
    row=idx;
    
    for j=1:m
        col=col+1;
        K(row,col)=1;
        row=row+n;
    end
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
