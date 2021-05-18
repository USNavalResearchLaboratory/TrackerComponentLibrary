function [x,z]=constrainedLSEq(A,b,B,d,algorithm)
%%CONSTRAINEDLSEQ Find x to minimize norm(A*x-b)^2 under the constraint
%                 that B*x=d. It is assumed that the (m+p)Xn stacked matrix
%                 [A;B] has linearly independent columns (left invertible)
%                 and B has linearly independent rows (right invertible).
%                 This means that p<=n<=m+p. This is equality constrained
%                 least squares.
%
%INPUTS: A A real mXn matrix
%        b A real mX1 vector.
%        B A real pXn matrix.
%        d A real pX1 vector.
% algorithm An optional prameter specfying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Algorithm  6.2.2 in Chapter 6.2.3 of [1].
%          1 In Chapter 14 of [2], it is shown that x solves the
%            constrained least squares problem if and only if there exists
%            a z such that [A'*A, B'; B, zeros(p,p)]*[x;z]=[A'*b;d]
%            This is a set of n+p equations and unknowns. This solution
%            just solves that system.
%          2 Use the algorithm of Chapter 6.2.5 of [1], which makes use of
%            the generalized singular value decomposition (GSVD).
%
%OUTPUTS: x The solution to the equality constrained least squares problem.
%           An mX1 vector.
%         z If algorithm=0, this is an (n-p)X1 unconstrained solution that
%           is an intermediate result. If algorithm=1, this is the pX1 set
%           of Lagrangian multipliers used in the solution. If algorithm=2,
%           this is an empty matrix.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%[2] S. Boyd and L. Vandenberghe, Vectors, Matrices, and Least Squares,
%    Sep. 2015, draft version of book. [Online].
%    Available: http://www.seas.ucla.edu/~vandenbe/133A/133A-textbook.pdf
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%Algorithm 6.2.2 in Chapter 6.2.3 of [1].
        [~,n1]=size(A);
        [m2,~]=size(B);

        [Q,R]=qr(B');
        opts.UT=false;
        opts.LT=true;
        opts.RECT=false;
        y=linsolve(R(1:m2,1:m2)',d,opts);

        A=A*Q;
        z=linearLeastSquares(A(:,(m2+1):n1),b-A(:,1:m2)*y);
        x=Q(:,1:m2)*y+Q(:,(m2+1):n1)*z;
    case 1%From Chapter 14 of [2].
        m=size(A,2);
        p=size(B,1);
        H=A'*A;
        c=A'*b;
        xz=[H,B';
            B,zeros(p,p)]\[c;d];

        x=xz(1:m);
        z=xz((m+1):end); %The Lagrangian multipliers.
    case 2%The GSVD approach of Chapter 6.2.5 of [1].
        [~,n1]=size(A);
        [m2,~]=size(B);

        [U1,U2,X,DA,DB]=gsvd(A,B);
        bTilde=U1'*b;
        dTilde=U2'*d;

        %The definition of the GSVD used in Chapter 6.1.6 of [1] is not the
        %same as that implemented by Matlab's gsvd function. Specifically,
        %we have to modify X as follows:
        X=inv(X');

        alphaVals=diag(DA);
        betaVals=diag(DB);

        x=zeros(n1,1);
        for i=1:m2
            x=x+(dTilde(i)/betaVals(i))*X(:,i);
        end

        for i=(m2+1):n1
            x=x+(bTilde(i)/alphaVals(i))*X(:,i);
        end
        z=[];
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
