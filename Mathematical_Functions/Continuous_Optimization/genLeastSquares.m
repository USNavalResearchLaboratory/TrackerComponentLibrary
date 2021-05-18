function x=genLeastSquares(A,b,B,algorithm,matType)
%%GENLEASTSQUARES Solve the generalized least squares problem, also known
%            as the weighted least squares problem. This finds the vector x
%            that minimizes (A*x-b)'*inv(B*B')*(A*x-b), where B is a lower-
%            triangular matrix. This is the same as norm(inv(B)*(A*x-b)),
%            which is also equivalent to minimizing v'*v under the
%            constraint that b=A*x+B*v.
%
%INPUTS: A An mXn matrix with m>=n. It is assumed that A has a rank of at
%          least n.
%        b An mX1 vector.
%        B An mXm lower-triangular, symmetric matrix. If matType=1, then
%          this is actually B*B' given above.
% algorithm An optional parameter selecting the algorithm to use. Possible
%         values are:
%         0 (The default if omitted or an empty matrix is passed). Use the
%           algorithm of Chapter 6.1.2 of [1], which makes use of the QR
%           decomposition and assumes that A and B are not rank deficient.
%           This is a variant of what is decribed in [2] and [4].
%         1 This is kind of the algorithm of [2], but using the augmented
%           matrix expression in [3]. This is not efficient.
%         2 The simplest direct approach. An equivlaent system replaces A
%           with B\A and b with B\b. Solve that system and transform back.
%         3 In [4], it is shown that the problem is equivalent to an
%           equality-constrained least squares problem. This solution feeds
%           that formulation to the constrainedLSEq function. This is
%           typicaly the slowest approach.
% matType A parameter specifying whether B is the B matrix above (0) or if
%         it is actually B*B' in the formulation above. The default if
%         omitted or an empty matrix is passed is 0.
%
%OUTPUTS: x The nX1 estimate.
%
%Note that [4] also gives an approach to handling the rank deficient case,
%which is not implemented here.
%
%%EXAMPLE:
%Here, we show that the different methods produce results on a simple
%random problem that are pretty close (within finite precision limits), but
%that the execution time can differ a lot.
% m=1000;
% n=900;
% A=randn(m,n);
% A(end,:)=A(1,:);
% xTrue=randn(n,1);
% b=A*xTrue+0*randn(m,1);
% 
% R=randOrthoMat(m);
% W=R*diag(10*rand(m,1))*R';
% B=chol(W,'lower');
% 
% tic
% x0=genLeastSquares(A,b,B,0);
% toc
% tic
% x1=genLeastSquares(A,b,B,1);
% toc
% tic
% x2=genLeastSquares(A,b,B,2);
% toc
% tic
% x3=genLeastSquares(A,b,B,3);
% toc
% 
% %All have comparable error:
% norm(B\(A*x0-b))
% norm(B\(A*x1-b))
% norm(B\(A*x2-b))
% norm(B\(A*x3-b))
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%[2] C. C. Paige, "Fast numerically stable computations for generalized
%    linear least squares problems," SIAM Journal on Numerical Analysis,
%    vol. 16, no. 1, pp. 165-171, Feb. 1979.
%[3] J. Y. Yuan, "Numerical methods for generalized least squares
%    problems," Journal of Computational and Applied Mathematics, vol. 66,
%    no. 1-2, pp. 571-584, 31 Jan. 1996.
%[4] C. C. Paige, "Computer solution and perturbation analysis of
%    generalized linear least squares problems," Mathematics of
%    Computation, vol. 33, no. 145, pp. 171-183, Jan. 1979.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[m,n]=size(A);

if(nargin<4||isempty(algorithm))
    algorithm=0;
end

if(nargin<5||isempty(matType))
    matType=0;
end
     
if(m<n)
    error('This algorithm requires m>=n.');
end

if(matType==1)
    B=chol(B,'lower');
end

switch(algorithm)
    case 0%The algorithm of [1], which is a variant of [2].

        %Take the QR decomposition of Q and find the submatrices specified
        %in Chapter 6.1.2 of [1].
        [Q,R]=qr(A);
        Q1=Q(:,1:n);
        Q2=Q(:,(n+1):m);
        R1=R(1:n,1:n);
        
        %Get the orthogonal matrix Z and the (m-n column) matrix S.
        tempMat=B'*Q2;
        tempMat=tempMat(:,(m-n):-1:1);
        [Z,S]=qr(tempMat);
        S=S((m-n):-1:1,(m-n):-1:1)';
        Z2=Z(:,(m-n):-1:1);

        opts.UT=true;
        opts.LT=false;
        opts.RECT=false;
        
        %Equation 6.1.9 in [1].
        %u=S\(Q2'*b);
        u=linsolve(S,Q2'*b,opts);
        
        %Equation 6.1.10 in [1].
        %x=R1\(Q1'*(b-B*(Z2*u)));
        x=linsolve(R1,Q1'*(b-B*(Z2*u)),opts);
    case 1%The algorithm of [2] as expressed in [3].
        c=[zeros(m,m), B, A;
           B', -eye(m), zeros(m,n);
           A', zeros(n,m), zeros(n,n)]\[b;zeros(m+n,1)];

        x=c((2*m+1):(2*m+n));
    case 2%The "standard" direct way.
        A=B\A;
        b=B\b;
        x=A\b;  
    case 3%The reformulation as an equality constrained problem, solved
           %using constrainedLSEq.
        AAlt=[zeros(m,n),eye(m,m)];
        bAlt=zeros(m,1);
        BAlt=[A,B];
        
        xv=constrainedLSEq(AAlt,bAlt,BAlt,b);
        x=xv(1:n);
    otherwise
        error('Unknown algorithm specified.')
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
