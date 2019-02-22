function [X,E0,R0]=totalLeastSquares(A,B,d,t)
%%TOTALLEASTSQUARES Solve the total least squares problem. That is, given
%          an mXn matrix A and an mXk matrix or vector B, find the an nXk
%          matrix or vector X such that (A+E0)*X=(B+R0), where the error
%          terms E0 and R0 are chosen to minimize the Frobenius norm of
%          D*[E0,R0]*T where D=diag(d) and T=diag(t) are nonsingular
%          diagonal matrices. This problem is essentially just solving
%          A*X=B for X where both A and B can be noisy, whereas the normal
%          least squares problem solves for the case where only B is noisy.
%
%INPUTS: A A real mXn matrix.
%        B A real nXk matrix.
%        d Optional  mX1 vector of real diagonal elements for the matrix D,
%          which affects how the errors are weighted. If omitted, a vector
%          of ones is used (D is just an identity matrix).
%        t Optional (n+k)X1 vector of real diagonal elements for the matrix
%          T, which affects show errors are weighted. If omitted, a vector
%          of ones is used (T is just an identity matrix).
%
%OUTPUTS: X The nXk real solution of the total least squares problem.
%     E0,R0 The error matrices associated with the solution to the problem.
%           They are such that (A+E0)*X=(B+R0) and the Frobenius norm of
%           D*[E0,R0]*T, which is computed via norm(D*[E0,R0]*T,'fro'), is
%           minimized.
%
%The algorithm for solving the total least squares problem is taken from
%Chapter 6.3.1 of [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
n=size(A,2);
k=size(B,2);

if(nargin<4)
    t=ones(n+k,1);
end

if(nargin<3)
    d=ones(m,1); 
end

D=diag(d);
T=diag(t);

[U,Sigma,V]=svd(D*[A,B]*T);

V12=V(1:n,(n+1:end));
V22=V((n+1):end,(n+1):end);

T1=diag(t(1:n));
T2=diag(t((n+1):end));

X=-T1*V12/V22/T2;

if(nargout>1)
    U2=U(:,(n+1):end);
    Sigma2=Sigma((n+1):end,(n+1):end);

    ERMat=-D\(U2*Sigma2*[V12',V22'])/T;

    E0=ERMat(:,1:n);
    R0=ERMat(:,(n+1):end);
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
