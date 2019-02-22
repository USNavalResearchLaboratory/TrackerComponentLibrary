function xEst=solveLinSysTikhonov(A,b,alpha)
%%SOLVELINSYSTIKHONOV When solving a possibly overdetermined system A*x=b,
%               one typically finds x to minimize norm(A*x-b)^2. However,
%               that can be subject to errors if A is badly conditioned.
%               This function avoids the finite precision issues by instead
%               solving the problem norm(A*x-b)^2+alpha^2*norm(x)^2 for a
%               fixed real alpha. As alpha->0, the problem approaches the
%               traditional least-squares solution. This method is known as
%               Tikhoniv regularization.
%
%INPUTS: A An mXn matrix with m>=n. This can be complex.
%        b AnmX1 vector.
%    alpha A real parameter.
%
%OUTPUTS: xEst The value of x minimizing norm(A*x)^2+alpha^2+norm(x)^2.
%
%This method of avoiding numerical instability in solving A*X=b is taken
%from Chapter 1.4.7 of [1].
%
%REFERENCES:
%[1] W. Wasylkiwskyj, Signal and Transforms in Linear Systems Analysis.
%    Heidelberg: Springer, 2013.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[UR,SigmaR,VR]=svd(A,0);

sigmaR=diag(SigmaR);

%Equation 1.200
xEst=VR*diag(sigmaR./(sigmaR.^2+alpha^2))*UR'*b;
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
