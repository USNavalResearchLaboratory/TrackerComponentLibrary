function X=RiccatiSolveC(A,B,Q,R,S,E)
%%RICCATISOLVEC Solve the general continuous-time Riccati equation.
%
%INPUTS: A An nXn matrix.
%        B An nXm matrix, where m<=n.
%        Q An nXn matrix such that Q=Q' and all eigenvalues are non-
%          negative.
%        R An optional mXm matrix such that R=R' and all eigenvalues are
%          non-negative. If omitted, eye(m,m) is used.
%        S An optional nXm matrix. If omitted, zeros(n,n) is used.
%        E An optional nXn matrix. If omitted, eye(n) is used.
%
%OUTPUTS: X The nXn nonnegative definite solution to the discrete-time
%           algebraic Ricatti equation.
%
%This function finds the nonnegative definite solution to the 
%continuous-time Ricatti equation having the form
%A'*X*E+E'*X*A-(E'*X*B+S)*inv(R)*(B'*X*E+S')+Q=0
%or, if R, S, and E are omitted, the equation under consideration becomes
%A'*X+X*A-X*B*B'*X+Q=0
%The algorithm of  [1] is used. Note that this function does not work with
%problems of the form -X*B*B'*X+Q=0 due to numerical issues.
%
%The continuous-time Ricatti equation arises when solving for the  steady-
%state covariance of a continuous-time linear Kalman filter, as described
%in Chapter 9.2.3 of [2].
%
%REFERENCES:
%[1] W. F. Arnold III and A. J. Laub, "Generalized eigenproblem algorithms
%    and software for algebraic Riccati equations," Proceedings of the
%    IEEE, vol. 72, no. 12, pp. 1746-1754, Dec. 1984.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);
m=size(B,2);

if(nargin<4)
    R=eye(m,m);
end

if(nargin<5)
    S=zeros(n,m);
end

if(nargin<6)
   E=eye(n); 
end

L=[E,            zeros(n,n),zeros(n,m);
   zeros(n,n),   E',        zeros(n,m);
   zeros(m,n),   zeros(m,n),zeros(m,m)];
M=[A,   zeros(n,n),  B;
  -Q,  -A',         -S;
   S',  B',          R];

[MHat,LHat,V,U]=qz(M,L,'real');
[~,~,~,U]=ordqz(MHat,LHat,V,U,'lhp');

W=[E,           zeros(n,n),zeros(n,m);
   zeros(n,n),  eye(n),zeros(n,m)]*U;

W11=W(1:n,1:n);
W21=W((n+1):(2*n),1:n);

X=W21/W11;
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

