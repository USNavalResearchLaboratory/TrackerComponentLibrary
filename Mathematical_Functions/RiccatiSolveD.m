function X=RiccatiSolveD(A,B,Q,R,S,E)
%%RICCATISOLVED Solve the general discrete-time Riccati equation.
%
%INPUTS: A An nXn matrix.
%        B An nXm matrix, where m<=n.
%        Q An nXn matrix such that Q=Q' and all eigenvalues are non-
%          negative.
%        R An mXm matrix such that R=R' and all eigenvalues are non-
%          negative.
%        S An optional nXm matrix. If omitted, zeros(n,m) is used.
%        E An optional nXn matrix. If omitted, eye(n) is used.
%
%OUTPUTS: X The nXn nonnegative definite solution to the discrete-time
%           algebraic Ricatti equation.
%
%This function finds the nonnegative definite solution to the discrete-time
%Ricatti equation having the form
%E'*X*E=A'*X*A-(A'*X*B+S)*inv(B'*X*B+R)*(A'*X*B+S)'+Q
%The algorithm of [1], which uses a qz decomposition, is used. Note that
%another explicit solution that uses an eigenvalue decomposition is given
%in [3].
%
%The discrete-time Ricatti equation arises when solving for the  steady-
%state covariance of a continuous-time linear Kalman filter, as described
%in Chapter 5.2.5 of [2].
%
%REFERENCES:
%[1] W. F. Arnold III and A. J. Laub, "Generalized eigenproblem algorithms
%    and software for algebraic Riccati equations," Proceedings of the
%    IEEE, vol. 72, no. 12, pp. 1746-1754, Dec. 1984.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[3] D. R. Vaughan, "A nonrecursive algebraic solution for the discrete
%    Riccati equation," IEEE Transactions on Automatic Control, vol. 15,
%    no. 5, pp. 597-599, Oct. 1970.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);
m=size(B,2);

if(nargin<5)
    S=zeros(n,m);
end

if(nargin<6)
   E=eye(n); 
end

L=[E,            zeros(n,n),zeros(n,m);
   zeros(n,n),   A',        zeros(n,m);
   zeros(m,n),  -B',        zeros(m,m)];
M=[A,  zeros(n,n),  B;
   -Q, E',         -S;
   S', zeros(m,n), R];

%We use complex for stability (otherwise ordqz can sometimes fail), but we
%have to use a "real" on the final solution, because it should be real.
[MHat,LHat,V,U]=qz(M,L,'complex');
[~,~,~,U]=ordqz(MHat,LHat,V,U,'udi');

W=[E,           zeros(n,n),zeros(n,m);
   zeros(n,n),  eye(n),zeros(n,m)]*U;

W11=W(1:n,1:n);
W21=W((n+1):(2*n),1:n);

%W21/W11; Use a pseudoinverse for stability.
X=real(W21*pinv(W11));

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
