function [B,Lambda,CLS,exitCode]=jointMatDiagLS(A,w,maxIter,CLSTol,RelTol,AbsTol)
%%JOINTMATDIAGLS Given K real or complex Hermitian or symmetric matrices in
%              A, try to find a common diagonalizing matrix B and sets of
%              values Lambda such that A(:,:,k)=B*diag(Lambda(:,k))*B' for
%              all k. An exact solution is only guaranteed for positive
%              (semi-)definite matrices with K=2. However, for general
%              problems the B and Lambda values are chosen trying to
%              minimize the least squares cost function
%              CLS=\sum_{k=1}^K w(k)*norm(A_k-B*diag(lambda_k)*B','fro')^2 
%              when the matrices A are Hermitian, lambda_k is an NX1 vector
%              and 
%              CLS=\sum_{k=1}^K w(k)*norm(A_k-B*diag(lambda_k)*B.','fro')^2
%              when the matrices A re symmetric. Real matrices count as
%              Hermitian; one cannot mix Hermitian and Symmetric complex
%              matrices in a single A.
%
%INPUTS: A An NXNXK set of K NXN Hermitian or symmetric matrices. This
%          cannot contain Inf or NaN values.
%        w An optional set of positive that go into the cost function.
%          These do not need to sum to any particular value. The default if
%          omitted or an empty matrix is passed is ones(K,1).
%  maxIter The maximum number of iterations to perform. The default if
%          omitted or an empty matrix is passed is 10000.
%   CLSTol If the least squares cost is <= this threshold, then the
%          algorithm will terminate early. This value is checked every 10
%          iterations. The default if omitted or an empty matrix is passed
%          is 1e-6.
% RelTol, AbsTol Absolute and relative tolerances on the parameters being
%          estimated. let theta=[B(:);Lambda(:)] and thetaPrev the estimate
%          from the previous step. Say that diff=abs(theta-thetaPrev).
%          Convergence is declared if
%          all((diff<=AbsTol)|diff<=RelTol*abs(theta))
%          The default values if omitted or empty matrices are passed are
%          1e-9 and 1e-12.
%
%OUTPUTS: B The NXN matrix that approximately turns the Lambda values into
%           A. These are the diagonalizer values. B is generally not
%           an orthonormal matrix.
%    Lambda The NXK set of the diagonals of the K diagonal matrices.
%       CLS The value of the cost function on termination.
%  exitCode A parameter indicating how the function terminated. Possible
%           values are:
%           -1 The RelTol and/or AbsTol criteria were satisfied.
%            0 The CLSTol criterion was satisfied.
%            1 The maximum number of iterations was reached.
%
%This function implements the algorithm of [1]. Though exact solutions are
%possible for pairs of positive semidefinite matrices, note that
%convergence is often faster for pairs of positive definite matrices than
%for pairs of positive semidefinite matrices, particularly when the
%matrices are complex.
%
%EXAMPLE 1:
%Here, we show that we can obtain an exact solution for a pair of positive
%definite real, symmetric (Hermitian) matrices where one is positive
%semi-definite.
% C1=[57,  7, 12, 17;
%      7, 47, 17, 22;
%     12, 17, 37, 27;
%     17, 22, 27, 27];
% [V,D,U] = svd(C1);
% D(2,2)=0;%Make uninformative.
% D(3,3)=0;%Make uninformative.
% C1=V*D*U';
% C2=[87, 25, 18, 31;
%     25, 63, 20, 17;
%     18, 20, 65, 29;
%     31, 17, 29, 65];
% A=zeros(4,4,2);
% A(:,:,1)=C1;
% A(:,:,2)=C2;
% [B,Lambda,CLS,exitCode]=jointMatDiagLS(A)
% %The CLS is verified as:
% CLSCheck=norm(A(:,:,1)-B*diag(Lambda(:,1))*B','fro')^2+norm(A(:,:,2)-B*diag(Lambda(:,2))*B','fro')^2
%One will see that convergence is by the CLS bound.
%
%EXAMPLE 2:
%Here, we have two positive definite complex Hermitian matrices. Again, an
%exact solution is found rapidly. If the additive identity matrices to C1
%and C2 are removed, then the matrices are no longer positive definite and
%only a least-squares approximation can be found. In such an instance, the
%algorithm is much slower.
% C1=[  9+  0*1i,   -65+  0*1i,  -11-153*1i,  -91-173*1i;
%     -65+  0*1i,    83+  0*1i,   54- 38*1i,   31+ 28*1i;
%     -11+153*1i,    54+ 38*1i,  130+  0*1i,   16- 47*1i;
%     -91+173*1i,    31- 28*1i,   16+ 47*1i,   22+  0*1i]+215*eye(4);
% 
% C2=[-16+  0*1i,  -32- 56*1i,  -12-128*1i,   16+114*1i;
%     -32+ 56*1i,   79+  0*1i,  -87- 67*1i,  -48+ 51*1i;
%     -12+128*1i,  -87+ 67*1i,   76+  0*1i,   -7- 96*1i;
%      16-114*1i,  -48- 51*1i,   -7+ 96*1i, -147+  0*1i]+259*eye(4);
% A=zeros(4,4,2);
% A(:,:,1)=C1;
% A(:,:,2)=C2;
% [B,Lambda,CLS,exitCode]=jointMatDiagLS(A)
% %The CLS is verified as:
% CLSCheck=norm(A(:,:,1)-B*diag(Lambda(:,1))*B','fro')^2+norm(A(:,:,2)-B*diag(Lambda(:,2))*B','fro')^2
%Again, convergence is due to satisfying the cost function criterion.
%
%EXAMPLE 3:
%Though only approximations can be obtained for more than two matrices in
%general, here we show that if three matrices originate from a common set
%of basis vectors, then this function will find the correct
%diagonalization.
% [V,D]=eig(magic(4));
% C1=V*D*V';
% C2=V*diag([1;2;3;5])*V';
% C3=V*diag([12;1;0;8])*V';
% A=zeros(4,4,3);
% A(:,:,1)=C1;
% A(:,:,2)=C2;
% A(:,:,3)=C3;
% [B,Lambda,CLS,exitCode]=jointMatDiagLS(A)
%Convergence of the cost function is quickly achieved for the three
%matrices. Note, however, that B is not the same as V.
%
%REFERENCES:
%[1] A. Yeredor, "Non-orthogonal joint diagonalization in the least-squares
%    sense with application in blind source separation," IEEE Transactions
%    on Signal Processing, vol. 50, no. 7, pp. 1545-1553, Jul. 2002.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    K=size(A,3);
    N=size(A,1);

    if(nargin<6||isempty(AbsTol))
        AbsTol=1e-12;
    end
    
    if(nargin<5||isempty(RelTol))
        RelTol=1e-9;
    end

    if(nargin<4||isempty(CLSTol))
        CLSTol=1e-6;
    end
    
    if(nargin<3||isempty(maxIter))
        maxIter=10000;
    end

    if(nargin<2||isempty(w))
        w=ones(K,1);
    end

    %Check whether all (assume symmetric) matrices are Hermitian. If so,
    %then a different algorithm can be used.
    isHermitian=true;
    if(~isreal(A))
        for k=1:K
            if(any(any(A(:,:,k)~=A(:,:,k)')))
                isHermitian=false;
                break;
            end
        end
    end

    w=w(:);

    %Initial values. This doesn't necessarily help speed up the
    %convergence, but it will work even if the matrices in A are singular.
    B=eye(N,N);
    Lambda=ones(N,K);

    exitCode=1;
    if(isHermitian)
        for curIter=1:maxIter
            thetaPrev=[B(:);Lambda(:)];
            
            for iter=1:4
                for l=1:N
                    B=ACPhase(A,Lambda,w,B,l);
                end
            end
            Lambda=DCPhase(A,B);

            %Check for convergence.
            if(mod(curIter,10))
                %Cost criterion.
                CLS=0;
                for k=1:K
                    CLS=CLS+w(k)*norm(A(:,:,k)-B*diag(Lambda(:,k))*B','fro')^2;
                end
                if(CLS<=CLSTol)
                    exitCode=0;
                    return;
                end
                
                %Parameter change criterion.
                theta=[B(:);Lambda(:)];
                
                diff=abs(theta-thetaPrev);
                if(all((diff<=AbsTol)|diff<=RelTol*abs(theta)))
                    exitCode=-1;
                    break;
                end
            end
        end
    else%Assume it is symmetric, not Hermitian.
        for curIter=1:maxIter
            thetaPrev=[B(:);Lambda(:)];
            
            for iter=1:4
                for l=1:N
                    B=ACPhaseSym(A,Lambda,w,B,l);
                end
            end
            Lambda=DCPhaseSym(A,B);

            %Check for convergence.
            if(mod(curIter,10))
                %Cost criterion.
                CLS=0;
                for k=1:K
                    CLS=CLS+w(k)*norm(A(:,:,k)-B*diag(Lambda(:,k))*B.','fro')^2;
                end
                if(CLS<=CLSTol)
                    exitCode=0;
                    return;
                end
                
                %Parameter change criterion.
                theta=[B(:);Lambda(:)];
                
                diff=abs(theta-thetaPrev);
                if(all((diff<=AbsTol)|diff<=RelTol*abs(theta)))
                    exigtCode=-1;
                    break;
                end
            end
        end
    end

    %If the algorithm did not convergence, compute the current CLS if it should
    %be returned.
    if(nargout>2)
        CLS=0;
        if(isHermitian)
            for k=1:K
                CLS=CLS+w(k)*norm(A(:,:,k)-B*diag(Lambda(:,k))*B','fro')^2;
            end
        else
            for k=1:K
                CLS=CLS+w(k)*norm(A(:,:,k)-B*diag(Lambda(:,k))*B.','fro')^2;
            end
        end
    end
end

function B=ACPhase(A,Lambda,w,B,l)

    N=size(A,1);
    K=size(A,3);

    %Step 1
    P=zeros(N,N);
    for k=1:K
        innerSum=B*diag(Lambda(:,k))*B'-Lambda(l,k)*(B(:,l)*B(:,l)');

        P=P+w(k)*Lambda(l,k)*(A(:,:,k)-innerSum);
    end

    %Step 2
    [V,D]=eig(P);
    d=diag(D);
    [~,idx]=max(d);
    mu=d(idx);
    betaVec=V(:,idx);

    %The paper suggests using the first nonzero index of betaVec. We just
    %use the index having the largest magnitude value. in it.
    [~,nonZeroIdx]=max(abs(betaVec));

    %Force the first nonzero element of the unit eigenvector to be real and
    %positive.
    if(isreal(betaVec))
        if(betaVec(nonZeroIdx)<0)
            betaVec=-betaVec;
        end
    else
        betaVec=exp(-angle(betaVec(nonZeroIdx))*1i)*betaVec;
        %betaVec(nonZeroIdx) should be real and positive, but this line
        %deals avoids any possible finite precision issues.
        betaVec(nonZeroIdx)=abs(betaVec(nonZeroIdx));
    end

    if(mu<0)
        B(:,l)=zeros(N,1);
    else
        B(:,l)=betaVec*sqrt(mu)/sqrt(sum(w.*(Lambda(l,:).^2).'));
    end
end

function Lambda=DCPhase(A,B)

    N=size(A,1);
    K=size(A,3);

    %Step 1
    prodVal=B'*B;
    prodVal=conj(prodVal).*prodVal;
    %Use a pseudoinverse for poorly conditioned problems.
    if(rcond(prodVal)<=1e-14)
        G=pinv(prodVal);
    else
        G=inv(prodVal);
    end

    %Step 2
    Lambda=zeros(N,K);
    for k=1:K
        Lambda(:,k)=G*diag(B'*A(:,:,k)*B);
    end
end

function B=ACPhaseSym(A,Lambda,w,B,l)

    N=size(A,1);
    K=size(A,3);

    %Step 1
    P=zeros(N,N);
    for k=1:K
        innerSum=A(:,:,k);
        for n=1:N
            if(n~=l)
                innerSum=innerSum-Lambda(n,k)*(B(:,n)*B(:,n).');
            end
        end

        P=P+w(k)*Lambda(l,k)*conj(innerSum);
    end

    %Step 2
    [V,D]=eig([real(P), -imag(P);
              -imag(P), -real(P)]);

    %D should be real. This line is just in case some finite precision
    %issue in matlab makes it not real.
    D=real(D);
    
    d=diag(D);
    [~,idx]=max(d);
    mu=d(idx);
    xi=V(:,idx);

    %Step 3
    if(mu<0)
        B(:,l)=zeros(N,1);
    else
        gammaVal=xi(1:N);
        delta=xi((N+1):(2*N));

        B(:,l)=(gammaVal+1j*delta)*sqrt(mu)/sqrt(sum(w.*(abs(Lambda(l,:)).^2).'));
    end
end

function Lambda=DCPhaseSym(A,B)

    N=size(A,1);
    K=size(A,3);

    %Step 1
    prodVal=B'*B;
    prodVal=prodVal.*prodVal;
    %Use a pseudoinverse for poorly conditioned problems.
    if(rcond(prodVal)<=1e-14)
        G=pinv(prodVal);
    else
        G=inv(prodVal);
    end

    %Step 2
    Lambda=zeros(N,K);
    for k=1:K
        Lambda(:,k)=G*diag(B'*A(:,:,k)*conj(B));
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
