function [row,col]=vechInd2Sub(n,idx)
%%VECHIND2SUB Given a linear index for a point in vech(A) (the vech
%             operation returns the diagonal and lower-half of A stacked
%             column-wise), return the row and column of A in which the
%             given element in vech(A) resides.
%
%INPUTS: n The dimension of the nXn matrix A in which one wants to map
%          linear indices of vech(A).
%      idx A matrix of the linear indices that are to be converted to row
%          and column indices.
%
%OUTPUTS: row, col Matrices of the row and column indices in A
%                  corresponding to the linear indices in vech(A) given by
%                  idx.
%
%The total number of elements taken by appending the diagonal and lower
%elements of the first C columns is a matrix is (C/2)*(1+2*n-C). This
%function first finds C and then the left over columns.
%
%Linear indices with values requiring the number of bits in a double
%floating point mantissa will fail due to finite precision limits. However,
%it is unlikely that matrices will be that large in most applications.
%
%Invalid idx values will lead to invalid (possibly complex) row,col values.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

col=ceil((1/2)+n-(1/2)*sqrt(1-8*idx+4*n+4*n^2));

%colDelta should be an integer.
colDelta=(1/2)*(1+2*n-col).*col;
row=n+(idx-colDelta);

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
