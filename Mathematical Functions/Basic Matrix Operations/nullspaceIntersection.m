function Y=nullspaceIntersection(A,B)
%%NULLSPACEINTERSECTION Given an mXn matrix A and a pXn matrix B, find a
%               basis for the intersection of the nullspaces of the
%               matrices. That is, the output Y is a vectors whose columns
%               are such that A*Y(:,i)=0 and B*Y(:,i) =0 for all i columns.
%               If there is no nullspace intersection, then return an empty
%               matrix.
%
%INPUTS: A An mXn matrix.
%        B A pXn matrix.
%
%OUTPUTS: Y A basis for the common nullspace of A and B. If there is no
%           common nullspace, then an empty matrix is returned.
%
%The common nullspace can be found by finidng the nullspace of an augmented
%matrix as described in Chapter 6.4.2 of [1]. The lower computational-
%complexity algorithm listed in that chapter is not used, because it was
%noticed that determining the rank of the intermediate matrix C was
%difficult to do accurately. For example, when trying to find the
%intersection of [A;A], comparing singular values using a tolerance
%similar to that used in Matlab's rank function would incorrectly lead to
%the conclusion that the matrix C was not rank 0 (even if A was not full
%rank), thus causing Y to be an empty matrix.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

Y=null([A;B]);
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
