function idx=sub2VechInd(n,row,col)
%%SUB2VECHIND Given a pair of indices for a value on the diagonal or lower
%            half triangular part of a matrix, obtain the corresponding
%            index of that value in vech(A) (the vech operation returns the
%            diagonal and lower-half of A stacked column-wise). The
%            function can also accept indices above the diagonal, on the
%            assumption that the matrix A is symmetric.
%
%INPUTS: n The dimension of the nXn matrix A in which one wants to map row,
%          col pairs into indices in vech(A) either with row, col all on or
%          below the main diagonal, or assuming that A is symmetric if
%          above the main diagonal.
%  row,col Matrices of pairs of indices to convert.
%
%OUTPUTS: idx The linear indices in vech(A) correspoinding to the elements
%             at the given rows and columns of a matrix A.
%
%The total number of elements taken by appending the diagonal and lower
%elements of the first C columns is a matrix is (C/2)*(1+2*n-C). This
%function first finds that value and derives the number of rows based on
%the column.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

row1=max(row,col);
col1=min(row,col);

idxBase=(1/2)*(1+2*n-(col1-1))*(col1-1);
idx=idxBase+row1-col1+1;

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
