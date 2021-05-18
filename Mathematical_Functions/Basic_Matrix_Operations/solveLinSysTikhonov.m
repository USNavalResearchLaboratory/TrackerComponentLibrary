function xEst=solveLinSysTikhonov(A,b,alpha,B)
%%SOLVELINSYSTIKHONOV When solving a possibly overdetermined system A*x=b,
%               one typically finds x to minimize norm(A*x-b)^2. However,
%               that can be subject to errors if A is badly conditioned.
%               This function avoids the finite precision issues by instead
%               solving the problem norm(A*x-b)^2+alpha^2*norm(x)^2 for a
%               fixed real alpha. As alpha->0, the problem approaches the
%               traditional least-squares solution. This method is known as
%               ridge regression or Tikhonov regularization. A generalized
%               form of the technique makes use of the matrix B, resulting
%               in the cost function norm(A*x-b)^2+alpha^2*norm(B*x)^2.
%
%INPUTS: A An mXn matrix with m>=n. This can be complex.
%        b An mX1 vector.
%    alpha A real parameter.
%        B An optional nXn matrix. If provided, a different algorithm is
%          used to solve the problem. If omitted or an empty matrix is
%          passed, the problem is equivalent to passing the identity
%          matrix.
%
%OUTPUTS: xEst The value of x minimizing norm(A*x)^2+alpha^2+norm(x)^2.
%
%The method of avoiding numerical instability in solving A*x=b is taken
%from Chapter 1.4.7 of [1]. This is the same as the singular value
%decomposition approach that is described in Chapter 6.1.4 of [2].
%
%If B is provided,  the augmented system is directly solved. A generalized
%singular value decomposition approach described in Section 6.1.6 of [2]
%might be better. Note that the gsvd function in Matlab is not the same as
%the generalized svd in [2], because the X produced by gsvd needs to be
%transformed as inv(X') to be the X that is in [2].
%
%Chapters 6.1.4 and 6.1.5 of [2], make a difference between ridge
%regression and Tikhonov resularization. In [2], ridge regression omits B
%and  Tikhonov regularization includes B.
%
%REFERENCES:
%[1] W. Wasylkiwskyj, Signal and Transforms in Linear Systems Analysis.
%    Heidelberg: Springer, 2013.
%[2] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>3&&~isempty(B))
    n=size(A,2);
    xEst=[A;alpha*B]\[b;zeros(n,1)];
else%Basic ridge regression.
    [UR,SigmaR,VR]=svd(A,0);
    sigmaR=diag(SigmaR);

    %Equation 1.200
    xEst=VR*diag(sigmaR./(sigmaR.^2+alpha^2))*UR'*b;
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
