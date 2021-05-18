function [R,L]=nullspace(A,rankAlg)
%%NULLSPACE Determine the right and left nullspaces of a matrix using a
%           singular value decomposition. The technique for determining the
%           rank (hence determining the number of vectors in the basis) can
%           be specified. A vector v is in the right nullspace if A*v=0
%           (within finite precision bounds). Similarly , a vector is in
%           the left nullspace if v*A=0.
%
%INPUTS: A The mXn matrix.
%  rankAlg An optional parameter specifying the algorithm used to determine
%          the rank of the matrix. This corresponds to the algorithm input
%          in the matrixRank function. If omitted or an empty matrix is
%          passed, then the default algorithm in matrixRank is used.
%
%OUTPUTS: R The right nullspace of A or an empty matrix if there is no
%           right nullspace.
%         L The left nullspace of A or an empty matrix if there is no
%           right nullspace.
%
%The relation of the SVD and the matrix is discussed in Chapter 5.4.1 of
%[1]. For a matrix A of rank r, the right singular vectors corresponding to
%the (n-r) smallest singular values form the nullspace. The same goes for
%the left singular vectors and the left nullspace.
%
%EXAMPLE 1:
% A=[9,  7,  6, 12;
%    4, 14, 15,  1];
% [R,L]=nullspace(A)
% A*R
%One will see that A*R is close to zero, within finite precision limits and
%that L is an empty matrix, indicating that there is no left nullspace.
%
%EXAMPLE 2:
%In this example, there are both left and right nullspaces.
%  A=[4, 3;
%     8, 6];
% [R,L]=nullspace(A)
% L*A
% A*R
%Both L*A and A*R will be close to zero, withing finite precision limits.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    rankAlg=[];
end

[rankVal,V,U]=matrixRank(A,rankAlg);
R=V(:,(rankVal+1):end);
L=U(:,(rankVal+1):end)';

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
