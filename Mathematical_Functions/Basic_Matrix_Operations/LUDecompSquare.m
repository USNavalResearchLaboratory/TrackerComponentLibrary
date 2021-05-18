function [L,U,piv1,piv2,sgn1,sgn2]=LUDecompSquare(A,algorithm)
%%LUDECOMPSQUARE Perform an LU decomposition of the matrix A. That is, find
%                L and U such that A=L*U or A=P'*L*U or A=P'*L*U*Q, where L
%                is a lower-triangular matrix, U is an upper-triangular
%                matrix, and P and Q are permutation matrices, assuming
%                that an algorithm is selected that produces piv1 and piv2
%                from which P and Q can be derived.
%
%INPUTS: A An nXn matrix. If algorithm=0 or algorithm=2, the matrix must be
%          positive definite.
% algorithm A parameter specifying the algorithm to use, which also affects
%          what is returned. Possible values are:
%          0 Use Algorithm 3.2.2 in [1]. This is for A=L*U.
%          1 Use Algorithm 3.4.2 in [1]. This is for A=P'*L*U. P can be
%            obtained from piv1 as described below.
%          2 Use Algorithm 3.2.1 in [1]. This is for A=L*U.
%          3 Use Algorithm 3.4.1 in [1]. This is for A=P'*L*U. 
%          4 Use Algorithm 3.4.3 in [1]. This is for A=P'*L*U*Q. P and Q
%            can be obtained from piv1 and piv2 as described below. This
%            algorithm can be much slower than algorithm 5.
%          5 Use the Algorithm in Section 3.4.7 in [1], which modifies
%            algorithm 3.4.3 of [1] and is faster. This is for A=P'*L*U*Q.
%
%OUTPUTS: L An nXn lower-triangular matrix.
%         U An nXn upper triangular matrix.
% piv1,piv2 If an algorithm producing these is used, these are permutation
%           vectors that are turned into the correspinding permutation
%           matrices as described below.
% sgn1,sgn2 If an even number of swaps went into forming piv1 (or piv1 is
%           empty), then sgn1 is 1; otherwise it is -1. The same goes for
%           sgn2 and piv2. These sign values play a role in the computation
%           of determinants.
%
%For algorithms with partial pivoting, piv1 is filled in and A=P'*L*U. For
%algorithms with full or rook pivoting, piv1 and piv2 are filled in and
%A=P'*L*U*Q. P and Q are found as
% P=eye(n,n);
% Q=eye(n,n);
% for k=1:(n-1)
%    P([k,piv1(k)],:)=P([piv1(k),k],:);
%    Q([k,piv2(k)],:)=Q([piv2(k),k],:);
% end
%
%EXAMPLE:
%We generate a large random matrix and perform an lu decomposition using
%all of the different algorithms and then display the maximum absolute
%error from each algorithm.
% n=500;
% A=rand(n,n);
% 
% [L0,U0]=LUDecompSquare(A,0);
% P0=eye(n,n);
% [L1,U1,piv1]=LUDecompSquare(A,1);
% P1=eye(n,n);
% for k=1:(n-1)
%    P1([k,piv1(k)],:)=P1([piv1(k),k],:); 
% end
% [L2,U2]=LUDecompSquare(A,2);
% P2=eye(n,n);
% [L3,U3,piv1]=LUDecompSquare(A,3);
% P3=eye(n,n);
% for k=1:(n-1)
%    P3([k,piv1(k)],:)=P3([piv1(k),k],:); 
% end
% [L4,U4,piv1,piv2]=LUDecompSquare(A,4);
% P4=eye(n,n);
% Q4=eye(n,n);
% for k=1:(n-1)
%    P4([k,piv1(k)],:)=P4([piv1(k),k],:);
%    Q4([k,piv2(k)],:)=Q4([piv2(k),k],:); 
% end
% [L5,U5,piv1,piv2]=LUDecompSquare(A,5);
% P5=eye(n,n);
% Q5=eye(n,n);
% for k=1:(n-1)
%    P5([k,piv1(k)],:)=P5([piv1(k),k],:);
%    Q5([k,piv2(k)],:)=Q5([piv2(k),k],:);
% end
% Alg1MaxErr0=max(abs(vec(A-P0'*L0*U0)))
% Alg1MaxErr1=max(abs(vec(A-P1'*L1*U1)))
% Alg1MaxErr2=max(abs(vec(A-P2'*L2*U2)))
% Alg1MaxErr3=max(abs(vec(A-P3'*L3*U3)))
% Alg1MaxErr4=max(abs(vec(A-P4'*L4*U4*Q4)))
% Alg1MaxErr5=max(abs(vec(A-P5'*L5*U5*Q5)))
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);
switch(algorithm)
    case 0%Algorithm 3.2.2 in [1] (gaxpy, no pivoting)
        L=eye(n,n);
        U=zeros(n,n);
        for j=1:n
            if(j==1)
                v=A(:,1);
            else
                a=A(:,j);
                %Solve linear system. L(1:j-1,1:(j-1))\a(1:(j-1));
                z=forwardSubstitution(L(1:j-1,1:(j-1)),a(1:(j-1)));

                U(1:(j-1),j)=z;
                v(j:n)=a(j:n)-L(j:n,1:(j-1))*z;
            end
            U(j,j)=v(j);
            L((j+1):n,j)=v((j+1):n)/v(j);
        end
        piv1=[];
        piv2=[];
        sgn1=1;
        sgn2=1;
    case 1%Algorithm 3.4.2 in [1] (gaxpy, partial pivoting)
        sgn1=1;
        L=eye(n,n);
        U=zeros(n,n);
        piv1=1:(n-1);
        for j=1:n
            if(j==1)
                v=A(:,1);
            else
                a=A(:,j);
                for k=1:(j-1)
                    a([piv1(k),k])=a([k,piv1(k)]);
                end
                %Solve the linear system. L(1:(j-1),1:(j-1))\a(1:(j-1));
                z=forwardSubstitution(L(1:(j-1),1:(j-1)),a(1:(j-1)),0);
                
                U(1:(j-1),j)=z;
                v(j:n)=a(j:n)-L(j:n,1:(j-1))*z;
            end
            if(j<n)
                [~,mu]=max(abs(v(j:n)));
                mu=mu+j-1;
                piv1(j)=mu;
                if(mu~=j)
                    sgn1=-sgn1;
                end

                v([j,piv1(j)])=v([piv1(j),j]);
                L([j,piv1(j)],1:(j-1))=L([piv1(j),j],1:(j-1));
            end
            U(j,j)=v(j);
            if(v(j)~=0)
                L((j+1):n,j)=v((j+1):n)/v(j);
            end
        end
        piv2=[];
        sgn2=1;
    case 2%Algorithm 3.2.1 in [1] (outer, no pivoting).
        for k=1:(n-1)
            rho=(k+1):n;
            A(rho,k)=A(rho,k)/A(k,k);
            A(rho,rho)=A(rho,rho)-A(rho,k)*A(k,rho);
        end
        L=eye(n,n)+tril(A,-1);
        U=triu(A);
        piv1=[];
        piv2=[];
        sgn1=1;
        sgn2=1;
    case 3%Algorithm 3.4.1 in [1] (outer, partial pivoting).
        sgn1=1;
        piv1=1:(n-1);
        for k=1:(n-1)
            [~,mu]=max(abs(A(k:n,k))); 
            mu=k-1+mu; 
            piv1(k)=mu; 
            if(mu~=k)
                sgn1=-sgn1;
            end
            
            A([k,mu],:)=A([mu,k],:);
            
            if(A(k,k)~=0)
                rho=(k+1):n;
                A(rho,k)=A(rho,k)/A(k,k);
                A(rho,rho) = A(rho,rho) - A(rho,k)*A(k,rho);
            end
        end
        L=eye(n,n)+tril(A,-1);
        U=triu(A);
        piv2=[];
        sgn2=1;
    case 4%Algorithm 3.4.3 in [1] (outer, full pivoting).
        sgn1=1;
        sgn2=1;
        piv1=1:(n-1);
        piv2=1:(n-1);
        for k=1:(n-1)
            [mCol,muCol]=max(abs(A(k:n,k:n))); 
            [~,tau]=max(mCol);
            piv1(k)=muCol(tau)+k-1;
            piv2(k)=tau+k-1;
            
            if(k~=piv1(k))
                sgn1=-sgn1;
            end

            if(k~=piv2(k))
                sgn2=-sgn2;
            end
            
            A([piv1(k),k],:)=A([k,piv1(k)],:);
            A(:,[piv2(k),k])=A(:,[k,piv2(k)]);
            if(A(k,k)~=0)
                rho=(k+1):n;
                A(rho,k)=A(rho,k)/A(k,k);
                A(rho,rho)=A(rho,rho)-A(rho,k)*A(k,rho);
            end
        end
        L=eye(n,n)+tril(A,-1);
        U=triu(A);
    case 5%The Algorithm in Section 3.4.7 in [1] (outer, rook pivoting).
        %This modifies part of Algorithm 3.4.3 (outer, full pivoting).
        sgn1=1;
        sgn2=1;
        piv1=1:(n-1);
        piv2=1:(n-1);
        for k=1:(n-1)
            mu=k;
            lambda=k;
            rho=abs(A(k,k));
            s=0;
            while(rho<max(abs(A(k:n,lambda)))||rho<max(abs(A(mu,k:n))))
                if(mod(s,2)==0)
                    [~,tau]=max(abs(A(k:n,lambda)));
                    mu=tau+k-1;
                else
                    [~,tau]=max(abs(A(mu,k:n)));
                    lambda=tau+k-1;
                end
                rho=abs(A(mu,lambda));
                s=s+1;
            end
            piv1(k)=mu;
            piv2(k)=lambda;
            
            if(piv1(k)~=k)
                sgn1=-sgn1;
            end
            
            if(piv2(k)~=k)
                sgn2=-sgn2;
            end
            
            A([piv1(k),k],:)=A([k,piv1(k)],:);
            A(:,[piv2(k),k])=A(:,[k,piv2(k)]);
            if(A(k,k)~=0)
                rho=(k+1):n;
                A(rho,k)=A(rho,k)/A(k,k);
                A(rho,rho)=A(rho,rho)-A(rho,k)*A(k,rho);
            end
        end
        L=eye(n,n)+tril(A,-1);
        U=triu(A);
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
