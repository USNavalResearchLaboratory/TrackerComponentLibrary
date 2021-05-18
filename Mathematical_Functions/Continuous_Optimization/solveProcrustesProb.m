function Q=solveProcrustesProb(A,B)
%%SOLVEPROCRUSTESPROB Solve the procrustes problem. This means find the
%           orthogonal matrix Q with determinant +/-1 such that
%           norm(A-B*Q,'fro) is minimized.
%
%INPUTS: A An mXp matrix.
%        B An mXp matrix.
%
%OUTPUTS: Q The pXp orthogonal matrix minimizing norm(A-B*Q,'fro) and
%           having determinant magnitude 1.
%
%Algorithm 6.4.1 in [1] is used.
%
%EXAMPLE:
%In this case, we create A directly via an orthogonal transformation of B.
%We then show that this function finds the orthgonal transformation matrix.
% m=60;
% p=30;
% QTrue=randOrthoMat(p);
% B=randn(m,p);
% A=B*QTrue;
% Q=solveProcrustesProb(A,B);
% max(abs(vec(A-B*Q)))
%The error will be on the order of what one would expect with finite
%precision limits.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[U,~,V]=svd(B'*A);
Q=U*V';
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
