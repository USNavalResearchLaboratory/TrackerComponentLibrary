function L=eliminationMatrix(n)
%%ELIMINATIONMATRIX Find an elimination matrix of a particular
%          dimensionality. An elimination matrix is such that
%          L*A(:)=vech(A) for an nXn matrix A. The elimination matrix
%          selects just the diagonal and lower-triangular elements of A
%          stacked by column.
%
%INPUTS: n The dimensionality of the square nXn matrix that is to be
%          transformed.
%
%OUTPUTS: L An (n*(n-1)/2)Xn^2 elimination matrix for an nXn matrix.
%
%Basic properties of the elimination matrix are given in [1]. Some notable
%ones are that L*L'=eye((n*(n-1)/2),(n*(n-1)/2)).
%
%REFERENCES:
%[1] J. R. Magnus and H. Neudecker, "The elimination matrix: Some lemmas
%    and applications," SIAM Journal on Algebraic Discrete Methods, vol. 1,
%    no. 4, pp. 422-449, 1980.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

L=zeros(n*(n+1)/2,n^2);

rowOut=0;
for curCol=1:n
    for curRow=curCol:n
        rowOut=rowOut+1;
        
        idx=n*(curCol-1)+curRow;
        L(rowOut,idx)=1;
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
