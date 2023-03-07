function d = extendedDiag(A,k)
%EXTENDEDDIAG Returns the extended diagonal of a matrix determined by a
%             shift of indices 1:max(size(A)). In the event that A is not
%             square, the matrix is transposed so that the shortest
%             dimension is the column dimension and then the full diagonal
%             is computed on that matrix.
%
%INPUTS:
% A: An m-by-n matrix.
% k: An integer. Positive integers start at the top of column mod(k+1,n).
%    Negative integers start at the beginning of row mod(k+1,m). Zero
%    chooses the main diagonal. Note that no check is performed to ensure k
%    is in a particular interval.
%
%OUTPUTS:
% d: A max([m,n])-by-1 vector containing the extended diagonal of A.
%
%NOTE: If only the standard diagonal is desired, use MATLAB's built-in diag
%      function instead. For an m-by-n matrix, if d is to be of length
%      min([m,n]), then the diag function will provide the correct result.
%
%EXAMPLE 1: Computes some extended diagonals for a magic square.
% A = magic(5);
% d0 = extendedDiag(A,0);
% d3 = extendedDiag(A,3);
% dm3 = extendedDiag(A,-3);
% d8 = extendedDiag(A,8);
% dm8 = extendedDiag(A,-8);
% assert(all(diag(A)==d0))
% assert(all(d3==d8))
% assert(all(dm3==dm8))
%
%EXAMPLE 2: Computes extended diagonals for a non-square matrix.
% A = [1,2,3;4,5,6;7,8,9;10,11,12];
% d0 = extendedDiag(A,0);
% d3 = extendedDiag(A,2);
% dm3 = extendedDiag(A,-3);
% d8 = extendedDiag(A,5);
% dm8 = extendedDiag(A,-7);
% assert(all(diag(A)==d0(1:3)))
% assert(all(d3==d8))
% assert(all(dm3==dm8))
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if isempty(A)
    d = [];
    return
end

sz = size(A);
maxDim = max(sz);

if k<0
    swap = true;
    k = -k;
else
    swap = false;
end

d = zeros(maxDim,1);
for i = 1:maxDim
    if ~swap
        d(i) = A(mod(i-1,sz(1))+1,mod(i+k-1,sz(2))+1);
    else
        d(i) = A(mod(i+k-1,sz(1))+1,mod(i-1,sz(2))+1);
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
