function K=commutationMatrix(m,n)
%%COMMUTATIONMATRIX Return a commutation matrix K for an mXn matrix. This
%          means that given an mXn matrix A that K*A(:)=B(:) where
%          B=transpose(A).
%
%INPUTS: m,n The number of rows and the number of columns of the matrix
%            that is to be commutated by this commutation matrix.
%
%OUTPUTS: K An (m*n)X(m*n) commutation matrix.
%
%The matrix K is just a type of permutation matrix with 1's in the correct
%places.
%
%Properties of the commutation matrix are given in [1].
%
%REFERENCES:
%[1] J. R. Magnus and H. Neudecker, "The elimination matrix: Some lemmas
%    and applications," SIAM Journal on Algebraic Discrete Methods, vol. 1,
%    no. 4, pp. 422-449, 1980.
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
