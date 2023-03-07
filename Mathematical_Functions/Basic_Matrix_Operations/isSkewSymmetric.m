function boolVal=isSkewSymmetric(A)
%%ISSKEWSYMMETRIC This function returns true if a matrix is skew symmetric
%        (if real) or skew-Hermitian (if complex). A real matrix is skew
%        symmetric if A.'=-A. For complex matrices, the same concept
%        applies, one is looking for A'=-A and it is called skew Hermetian.
%        The only difference is the transpose turned into a Hermitian
%        operation. Here, we use the Hermitian operation to test for both.
% 
%INPUTS A An nXn matrix.
%
%OUTPUTS: boolVal This is true if A'=-A.
%
%EXAMPLE:
%Matrix A is skew symmetric; matrix B is not.
% A=[0, 2, 3,-9;
%   -2, 0,-1, 0;
%   -3, 1, 0,-3;
%    9, 0, 3, 0];
% B=abs(A);
% boolVal=isSkewSymmetric(A)
% boolVal=isSkewSymmetric(B)
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

boolVal=all(all(A==-A'));

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
