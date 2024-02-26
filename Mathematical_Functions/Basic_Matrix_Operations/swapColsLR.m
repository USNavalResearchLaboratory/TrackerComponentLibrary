function XSwapped=swapColsLR(X)
%%SWAPCOLSLR Given a 2D matrix, move columns from the left half of the
% matrix to the right half and from the right half to the left half. The
% ordering of the columns in the individual halves is not reversed, so this
% is not the same as the fliplr function. If there is an odd number of
% columns, then the middle column is considered part of the leftmost half
% of the matrix.
%
%INPUTS: X An MXN matrix.
%
%OUTPUTS: XSWapped The NXM matrix given by taking the first half of the
%                  columns of X, moving them to the end and moving the
%                  second half of the columns of X to the beginning.
%
%EXAMPLE:
%Here, we have a matrix with an even number of elements. We show the matrix
%before swapping and after; the columns have been swapped, left and right.
% X=magic(4)
% XSwapped=swapColsLR(X)
%
%August 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[~,N]=size(X);
n1=1:ceil(N/2);
n2=(ceil(N/2)+1):N;

XSwapped=[X(:,n2) X(:,n1)];

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
