function [A,CB,exitCode]=jointEig(M,maxIter)
%%JOINTEIG Given K NXN real matrices that are not defective, try to find a
%          set of joint eigenvectors. That is find A such that
%          inv(A)*M(:,:,k)*A is a diagonal matrix for all k. Because such
%          solutions cannot be found exactly for general matrices, a
%          minimization of the sum of the squared Frobenius norms of the
%          matrices E_k=inv(A)*M(:,:,k)*A after zeroing the diagonal
%          elements is performed.
%
%INPUTS: M An NXNXK set of K NXN real, non-defective matrices.
%  maxIter The maximum number of iterations to try to perform. The default
%          if omitted or an empty matrix is passed is 500.
%
%OUTPUTS: A The approximate NXN matrix of joint eigenvectors. 
%        CB The value of the pentalty term at termination.
%  exitCode A value indicating how the algorithm exited. Possible values
%           are:
%           0 The algorithm converged.
%           1 The maximum number of iterations was reached.
%           
%This implements the algorithm of [1]. Note that a matrix is considered
%defective if it does not have a complete basis of eigenvectors and is thus
%not diagonalizable. This type of joint eigenvalue decomposition arises in
%some independent component analysis (ICA) techniques.
%
%If CB is the cost (sum of Frobenius norm squared of off-diagonal elements
%of the approximately diaognalized matrices), and CBPrev is the cost at the
%previous iteration, then convergence is determined when
%abs(CB-prevCost)<=prevCost*eps()
%
%EXAMPLE 1:
%Here we get an exact solution for two matrices when there exists an exact
%common set of eigenvectors.
% C1=[87, 25, 18, 31;
%     25, 63, 20, 17;
%     18, 20, 65, 29;
%     31, 17, 29, 65];
% [V,D]=eig(C1);
% C2=V*diag([1;2;3;-10])*inv(V);
% M=zeros(4,4,2);
% M(:,:,1)=C1;
% M(:,:,2)=C2;
% [A,CB,exitCode]=jointEig(M)
% %Verify that the solutions are diagonal.
% inv(A)*C1*A
% inv(A)*C2*A
%
%EXAMPLE 2:
%Here, we perturb the matrices with some noise and see that reasonable
%diagonalization is still possible.
% C1=[87, 25, 18, 31;
%     25, 63, 20, 17;
%     18, 20, 65, 29;
%     31, 17, 29, 65];
% [V,D]=eig(C1);
% C2=V*diag([1;2;3;-10])*inv(V);
% C1=C1+1e-3*randn(4,4);
% C2=C2+1e-3*randn(4,4);
% M=zeros(4,4,2);
% M(:,:,1)=C1;
% M(:,:,2)=C2;
% [A,CB,exitCode]=jointEig(M)
% %Verify that the solutions are approximately diagonal desite the noise.
% inv(A)*C1*A
% inv(A)*C2*A
%
%REFERENCES:
%[1] R. André, T. Trainini, X. Luciani, and E. Moreau, "A fast algorithm
%    for joint eigenvalue decomposition of real matrices," in Proceedings
%    of the 23rd European Signal Processing Conference, Nice, France, 31
%    Aug.-4Sep.2015,pp.1316-1320.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maxIter))
    maxIter=1000;
end

N=size(M,1);
K=size(M,3);

%Allocate space.
Z=zeros(N,N);%Diagonals are initialized to zero.

prevCost=0;

diagIdx=diagElIdx(N);

exitCode=1;

T=M;%Initial value suggested after Equation 3.
A=eye(N,N);%Initialize with the identity matrix.
for curIter=1:maxIter
    for m=1:N
        for n=1:N
            if(m~=n)
                %Equation 15, using the definitions of O and Lambda from
                %Equation 9, where the kth Lambda matrix holds just the
                %diagonal elements of T(:,:,k) (The others are zeroed) and
                %O holds the non-diagonal elements (the others are zeroed).

                num=0;
                denom=0;
                for k=1:K
                    diffVal=T(m,m,k)-T(n,n,k);
                    num=num+T(m,n,k)*diffVal;
                    denom=denom+diffVal^2;
                end
                Z(m,n)=-num/denom;
                if(num==0)
                    Z(m,n)=0;
                end
            end
        end
    end
    
    B=eye(N,N)+Z;
    if(any(~isfinite(B(:))))
       error('Finite precision errors prevent continuation. One or more matrices in M might be defective.') 
    end
    for k=1:K
        T(:,:,k)=B\T(:,:,k)*B; 
    end
    
    A=A*B;

    %Compute the cost.
    CB=0;
    for k=1:K
        FMat=A\M(:,:,k)*A;
        FMat(diagIdx)=0;
        CB=CB+norm(FMat,'fro')^2;
    end

    %Check for convergence.
    if(abs(CB-prevCost)<=prevCost*eps())
        exitCode=0;
        break;
    end
    prevCost=CB;
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
