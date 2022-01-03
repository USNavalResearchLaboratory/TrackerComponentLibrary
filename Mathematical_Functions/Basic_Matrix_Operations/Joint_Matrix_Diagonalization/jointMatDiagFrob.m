function [V,C,FCost,exitCode]=jointMatDiagFrob(A,w,maxIter,RelTol,AbsTol)
%%JOINTMATDIAGFROB Given K real NXN symmetric matrices in A, try to find a
%               common diagonalization matrix V such that V*A(:,:,k)*V' is
%               diagonal for all k. V is generally not orthogonal. When
%               K>2, it is possible that no common matrix exists, in which
%               case, the sum of the squared Frobenius norms of the
%               matrices E_k for k=1:K is minimized. E_k is the matrix 
%               V*A(:,:,k)*V' with the main diagonal elements set to zero.
%               Alternatively, positive weights can be provided in w to
%               allow for the Frobenius norms to be weighred going into the
%               cost function. In the unweighted case, the fast algorithm
%               of [1] is used. In the weighted case, the the slower
%               algorithm of [2] is used.
%
%INPUTS: A An NXNXK set of K NXN real, symmetric matrices. This cannot
%          contain Inf or NaN values. The matrices do not have to be
%          positive definite.
%        w If all of the Frobenius norms in the cost function are meant to
%          be equally weighted, then an empty matrix can be passed for w
%          and the fast algorithm of [1] will be used. Otherwise, the
%          slower algorithm of [2] will be used.
%  maxIter The maximum number of iterations to perform. The default if
%          omitted or an empty matrix is passed is 1000.
% RelTol, AbsTol Tolerances on the cost function for declaring convergence.
%          Let E be the sum of the squared Frobenius norms of the matrices
%          E_k, where E_k is V*A(:,:,k)*V' with the diagonal elements
%          zeroed. Let D be the sum of the squared Frobenis norms of the
%          matrices D_k=diag(diag(V*A(:,:,k)*V')). Convergence is declared
%          if E/(N*(N-1))<=AbsTol or if E/(N*(N-1))<D/N*RelTol. The
%          inclusion of N deals with how many matrix elements are present.
%          The defaults if omitted or empty matrices are passed are
%          RelTol=eps() and AbsTol=eps().
%
%OUTPUTS: V The NXN matrix such that V*A(:,:,k)*V' is approximately
%           diagonal for all k.
%         C The set of matrices C=V*A(:,:,k)*V'.
%     FCost The sum of the squared Frobenius norms of the (approximately)
%           in A matrices after removing the diagonal elements divided by
%           (N*(N-1)). This is a termination cost criterion.
%  exitCode A parameter indicating how the algorithm terminated. Possible
%           values are:
%           0 The algorithm converged.
%           1 the maximum number of iterations elapsed.
%
%This function is an implementation of the algorithm described in [1].
%
%EXAMPLE 1:
%This is a demonstration given two positive semidefinite matrices. An exact
%solution is found.
% %C1 is positive semidefinite with one zero eigenvalue.
% C1=(magic(4)+magic(4)')/2;
% [V,D]=eig(C1);
% %C2 is positive semidefinite with two zero eigenvalues.
% C2=V*diag([1;0;0;8])*V';
% %Force numeric symmetry.
% C2=(C2+C2')/2;
% 
% A=zeros(4,4,2);
% A(:,:,1)=C1;
% A(:,:,2)=C2;
% [V,C,FCost,exitCode]=jointMatDiagFrob(A)
% %The algorithm quickly converges.
%
%EXAMPLE 2:
%Though only approximations can be obtained for more than two matrices in
%general, here we show that if three matrices originate from a common set
%of basis vectors, then this function will find the correct
%diagonalization.
% [V1,D]=eig(magic(4));
% C1=V1*D*V1';
% C2=V1*diag([1;2;3;5])*V1';
% C3=V1*diag([12;1;0;8])*V1';
% A=zeros(4,4,3);
% A(:,:,1)=C1;
% A(:,:,2)=C2;
% A(:,:,3)=C3;
% [V,C,FCost,exitCode]=jointMatDiagFrob(A)
%Convergence of the cost function is quickly achieved for the three
%matrices. Note, however, that the V obtained is not the same as V1.
%
%REFERENCES:
%[1] A. Ziehe, P. Laskov, G. Nolte, and K.-R. Müller, "A fast algorithm for
%    joint diagonalization with non-orthogonal transformations and its
%    application to blind source separation," Journal of Machine Learning,
%    vol. 5, pp. 777-800, Jul. 2004.
%[2] R. Vollgraf and K. Obermayer, "Quadratic optimization for simultaneous
%    matrix diagonalization," IEEE Transactions on Signal Processing, vol.
%    54, no. 9, pp. 3270-3278, Sep. 2006.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isreal(A))
    error('This algorithm only supports real matrices A.')
end

if(nargin<3||isempty(maxIter))
    maxIter=1000;
end

if(nargin<4||isempty(RelTol))
    RelTol=eps();
end

if(nargin<5||isempty(AbsTol))
    AbsTol=eps();
end

if(nargin>1&&~isempty(w)&&~all(w(1)==w))
    w=w/sum(w);
    [V,C,FCost,exitCode]=VollAlg(A,w,maxIter,RelTol,AbsTol);
else
    [V,C,FCost,exitCode]=ZieheAlg(A,maxIter,RelTol,AbsTol);
end
end

function [V,C,FCost,exitCode]=ZieheAlg(A,maxIter,RelTol,AbsTol)
%%ZIEHEALG This function implements the algorithm of [1] using Equation 17
%          instead of 18.
%
%REFERENCES:
%[1] A. Ziehe, P. Laskov, G. Nolte, and K.-R. Müller, "A fast algorithm for
%    joint diagonalization with non-orthogonal transformations and its
%    application to blind source separation," Journal of Machine Learning,
%    vol. 5, pp. 777-800, Jul. 2004.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

N=size(A,1);
K=size(A,3);

%The bound for the Frobenius norm of W. This is an arbitrary value <1.
theta=0.9;

I=eye(N,N);

%Initial values
W=zeros(N,N);
V=I;
C=A;

diagIdx=diagElIdx(N);
%With V initially being the identity matrix, we do not need to update C.
exitCode=1;
for curIter=1:maxIter
    %Implementation of y and z for Equation 17.
    z=zeros(N,N);
    y=zeros(N,N);
    for i=1:N
        for j=1:N
            for k=1:K
                z(i,j)=z(i,j)+C(i,i,k)*C(j,j,k);
                y(i,j)=y(i,j)+C(j,j,k)*(C(i,j,k)+C(j,i,k))/2;
            end
        end
    end
    
    %Equation 17
    for i=1:N
        for j=(i+1):N
            denom=z(j,j)*z(i,i)-z(i,j)^2;
            
            W(i,j)=(z(i,j)*y(j,i)-z(i,i)*y(i,j))/denom;
            W(j,i)=(z(i,j)*y(i,j)-z(j,j)*y(j,i))/denom;
        end
    end  

    %Limit the Frobenius norm of W.
    WNorm=norm(W,'fro');
    if(WNorm>theta)
         W=(theta/WNorm)*W;
    end
    
    %Equation 7
    V=(I+W)*V;

    for k=1:K
        C(:,:,k)=(I+W)*C(:,:,k)*(I+W)';
    end
    
    %Test for convergence.
    FCost=0;
    DCost=0;
    for k=1:K
        E=C(:,:,k);
        d=E(diagIdx);%diagonal elements.
        E(diagIdx)=0;
        
        %Cumulative squared Frobenius norm of the off-diagonal elements.
        FCost=FCost+norm(E,'fro')^2;
        
        %Cumulative squared Frobenius norm fo the diagonal elements.
        DCost=DCost+sum(d.^2);
    end
    DCost=DCost/N;
    FCost=FCost/(N*(N-1));

    if(FCost<=AbsTol||FCost<=RelTol*DCost)
        exitCode=0;
        break;%It converged.
    end
end
end

function [W,CD,FCost,exitCode]=VollAlg(A,alphaVals,maxIter,RelTol,AbsTol)
%%VOLLALG This function implements the algorithm of [1].
%
%REFERENCES:
%[1] R. Vollgraf and K. Obermayer, "Quadratic optimization for simultaneous
%    matrix diagonalization," IEEE Transactions on Signal Processing, vol.
%    54, no. 9, pp. 3270-3278, Sep. 2006.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

N=size(A,1);
K=size(A,3);

%We just set C0=eye(N,N) and thus P=eye(N,N) so P'*C0*P=eye(N,N);
%There does not seem to be a reason to randomize the W matrix, so we just
%set it to eye(N,N) and thus diag(W*C0*W')=1 for all i.
W=eye(N,N);
M=zeros(N,N);
for k=1:K
    %W is the identity matrix, so it can be omitted from M1 and M2.
    M1=A(:,:,k);
    %M2=A(:,:,k)';
   
    M1M1p=M1*M1';%M2*M2'=(M1*M1')'
    
    M=M+alphaVals(k)*(M1M1p+M1M1p');
end

diagIdx=diagElIdx(N);
exitCode=1;
CD=zeros(N,N,K);
for curIter=1:maxIter
    for i=1:N
        wi=W(:,i);
        for k=1:K
            m1=A(:,:,k)*wi;
            m2=A(:,:,k)'*wi;
            M=M-alphaVals(k)*(m1*m1'+m2*m2');
        end
        
        [V,D]=eig(M);
        d=diag(D);
        [~,idx]=min(abs(d));
        
        W(:,i)=V(:,idx);
        wi=W(:,i);
        
        for k=1:K
            m1=A(:,:,k)*wi;
            m2=A(:,:,k)'*wi;
            M=M+alphaVals(k)*(m1*m1'+m2*m2');
        end
    end
    
    %Compute the LW cost in Equation 4. We only do this every few
    %iterations so that it is faster.
    if(mod(curIter,25)==0)
        FCost=0;
        DCost=0;
        for k=1:K
            CD(:,:,k)=W'*A(:,:,k)*W;
            E=CD(:,:,k);
            d=E(diagIdx);%diagonal elements.
            E(diagIdx)=0;

            %Cumulative squared Frobenius norm of the off-diagonal elements.
            FCost=FCost+alphaVals(k)*norm(E,'fro')^2;

            %Cumulative squared Frobenius norm fo the diagonal elements.
            DCost=DCost+alphaVals(k)*sum(d.^2);
        end
        DCost=DCost/N;
        FCost=FCost/(N*(N-1));
        
        if(FCost<=AbsTol||FCost<=RelTol*DCost)
            exitCode=0;
            W=W';
            return;%It converged.
        end
    end
end

%Compute an up to date FCost
FCost=0;
for k=1:K
    CD(:,:,k)=W'*A(:,:,k)*W;
    E=CD(:,:,k);
    E(diagIdx)=0;

    %Cumulative squared Frobenius norm of the off-diagonal elements.
    FCost=FCost+alphaVals(k)*norm(E,'fro')^2;
end
FCost=FCost/(N*(N-1));

W=W';

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
