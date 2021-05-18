function X=solveQuadraticMatrixEq(A,B,C)
%%SOLVEQUADRATICMATRIXEQ Solve the quadratic matrix equation A*X*X+B*X+C=0.
%           for X. Only one solution is obtained.
%
%INPUTS: A, B, C The nXn matrices in the equation A*X*X+B*X+C=0. A B and C
%                can be complex.
%
%OUTPUTS: X The solution to the Equation A*X*X+B*X+C=0. If there is no
%           solution, X will generally be full of NaNs and Inf terms and a
%           warning about matrix singularity should occur. X can be complex
%           even if A, B, and C are real. Note that finite precision errors
%           might force X to be slightly complex even when it should be
%           real.
%
%This function implements the algorithm of [1], which solves the problem
%using a qz (generalized Schur) decompozition.
%
%REFERENCES:
%[1] N. J. Higham and H.-M. Kim, "Numerical analysis of a quadratic matrix
%    equation," IMA Journal of Numerical Analysis, vol. 20, no. 4, pp.
%    499-519, Oct. 2000.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

F=[zeros(n,n),eye(n,n);
   -C,        -B];
G=[eye(n,n),    zeros(n,n);
   zeros(n,n),  A];

%The QZ decomposition is the generalized Schur decomposition.
[~,~,~,Z]=qz(F,G,'complex');

Z11=Z(1:n,1:n);
Z21=Z((n+1):(2*n),1:n);

X=Z21/Z11;
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
