function [L,V,D,exitCode]=jointEigLU(M,maxIter)
%%JOINTEIGLU Given K NXN real or complex matrices that are not defective,
%          try to find a set of joint eigenvectors. That is find A such
%          that inv(A)*M(:,:,k)*A is a diagonal matrix for all k. Rather
%          than directly finding A, this function actually finds a lower
%          triangular matrix L, an upper triangular matrix V, and K
%          diagonal matrices such that inv(L)*M(:,:,k)*L=V*D(:,:,k)*inv(V),
%          or such that the diagonalization is as close as possible
%          for all M according to a cost function. The diagonalizing matrix
%          A is this A=L*V (which completes a generalized eigenvalue
%          decomposition). The values inv(L)*M(:,:,k)*L are also upper
%          triangular, and the cost function is based on the sum of the
%          squares of the lower triangular elements of inv(L)*M(:,:,k)*L
%          (See more below).
%
%INPUTS: M An NXNXK set of K NXN real, non-defective matrices.
%  maxIter The maximum number of iterations to try to perform. The default
%          if omitted or an empty matrix is passed is 500.
%
%OUTPUTS: L An NXN lower-triangular matrix such inv(L)*M(:,:,k)*L is
%           approximately upper triangular for all k.
%         V An NXN matrix such that V*D(:,:,k)*inv(V) is upper triangular
%           and is has ben optimized to be close to inv(L)*M(:,:,k)*L.
%         D An NXNXK set of K NXN diagonal matrices.
%  exitCode A value indicating how the algorithm exited. POssible values
%           are:
%           0 The algorithm converged.
%           1 The maximum number of iterations was reached.
%
%This implements the algorithm of [1] using the zeta_U cost function,
%which means that sequential updates are performed trying to drive each
%element of inv(L)*M(:,:,k)*L to a lower-triangular form. This is not a
%global cost function (unlike zeta_O in [1]), but was demonstrated to be
%better in simulation in [1]. However, a step of the optimization involves
%finding a real root of a polynomial. In Equation 19, one would tend to
%think that the root that should be taken would be the one totally
%minimizing the cost function. However, it was observed that this behavior
%can cause iterative cycles of worsening and improving global cost (global
%cost is the sum of the squares of all lower-triangular elements of
%inv(L)*M(:,:,k)*L) when the function is given certain structures matrices,
%such as magic(4). It was observed that this behaviour could be avoided by
%always taking the smallest zero of the derivative of zeta_U, which is not
%necessarily a global minimizer. Though heuristic, there is no explicit
%global convergence proof in [1], so there is no indication that this
%heuristic is theoretically worse.
%
%Note that a matrix is considered defective if it does not have a complete
%basis of eigenvectors and is thus not diagonalizable. This type of joint
%eigenvalue decomposition arises in some independent component analysis
%(ICA) techniques.
%
%If CB is the cost (sum of the squared magnitudes of the lower-triangular
%elements of inv(L)*M(:,:,k)*L summed across all K) and CBPrev is the cost
%at the previous iteration, then convergence is determined when
%abs(CB-prevCost)<=prevCost*eps()
%
%EXAMPLE:
%Here, we have three matrix, two of which are complex, that share a common
%basis. We demonstrate that this function is able to find a unique common
%basis (though not identical to the normalized basis used at the start).
% C1=[87, 25, 18, 31;
%     25, 63, 20, 17;
%     18, 20, 65, 29;
%     31, 17, 29, 65];
% [V,~]=eig(C1);
% C2=V*diag([1;2-10*1j;3;-10+1j])*inv(V);
% C3=V*diag([9;-18;-4+2*1j;-10])*inv(V);
% 
% M=zeros(4,4,3);
% M(:,:,1)=C1;
% M(:,:,2)=C2;
% M(:,:,3)=C3;
% M=real(M);
% [L,V,D,exitCode]=jointEigLU(M)
% A=L*V;
% 
% %Verify that the solutions are diagonal.
% inv(A)*C1*A
% inv(A)*C2*A
% inv(A)*C3*A
%
%REFERENCES:
%[1] X. Luciani and L. Albera, "Joint eigenvalue decomposition of non-
%   defective matrices based on the LU factorization with application to
%   ICA," IEEE transactions on Signal Processing, vol. 63, no. 17, pp.
%   4594-4608, 1 Sep. 2015.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maxIter))
    maxIter=500;
end

N=size(M,1);
K=size(M,3);

%This will be used to select all of the elements of M below the main
%diagonal for each matrix. This plays a role in determining convergence.
MSubEls=zeros(N,N,K);
for i=1:N
    for j=1:(i-1)
        MSubEls(i,j,:)=1;
    end
end

exitCode=1;
L=eye(N,N);%Initialize.
if(isreal(M))
    prevCost=Inf;
    for curIter=1:maxIter
        %Compute the polynomial coefficients of x and y from the appendix and
        %use them to directly evaluate Eij and update L and M.
        for j=1:(N-1)
            for i=(j+1):N
                %Get the coefficients from Appendix A. for zeta_U.
                xi=sum(M(j,i,:).^2);
                lambda=0;
                mu=0;
                nu=0;
                for k=1:K
                    MDiff=M(i,i,k)-M(j,j,k);

                    lambda=lambda+M(i,j,k)*MDiff;
                    mu=mu+MDiff^2-2*M(i,j,k)*M(j,i,k);
                    nu=nu+M(j,i,k)*MDiff;
                end
                lambda=2*lambda;
                nu=-2*nu;

                %Solve for the zeros of the derivative of Equation 47.
                zVals=roots([4*xi;3*nu;2*mu;lambda]);
                if(~isreal(zVals))
                    %Take the one real solution. That will be the one with
                    %the smallest absolute imaginary component.     
                    [~,idx]=min(abs(imag(zVals)));
                    z=zVals(idx);
                else%If there are three real solutions, take the minimum
                    %one (not necessarily the one minimizing the cost
                    %function).

                    z=min(zVals);
                end

                Eij=eye(N,N);
                Eij(i,j)=z;
                Eijn=eye(N,N);
                Eijn(i,j)=-z;
                L=L*Eij;%Update L
                %Update M
                for k=1:K
                    M(:,:,k)=Eijn*M(:,:,k)*Eij;
                end
            end
        end

        %Compute the cost, the sum of the squared magnitudes of the sub-
        %diagonal terms, then check for convergence.
        curCost=sum(MSubEls(:).*M(:).^2);
        diff=abs(curCost-prevCost);
        if(diff<prevCost*eps())
            exitCode=0;%It converged.
            break
        end
        prevCost=curCost;
    end
else
    prevCost=Inf;
    for curIter=1:maxIter
        for j=1:(N-1)
            for i=(j+1):N
                %Compute the polynomial coefficients of x to evaluate EIj
                %and update L and M. Then compute those of y to update L
                %and M separately.
                U=real(M);
                V=imag(M);
                
                %Get the coefficients from Appendix B for x
                epsilon=sum(U(j,i,:).^2+V(j,i,:).^2);
                beta=0;
                gamma=0;
                delta=0;
                for k=1:K
                    UDiff=U(i,i,k)-U(j,j,k);
                    VDiff=V(i,i,k)-V(j,j,k);

                    %Coefficients for x.
                    beta=beta+U(i,j,k)*UDiff+V(i,j,k)*VDiff;
                    gamma=gamma+UDiff^2-2*U(i,j,k)*U(j,i,k)+VDiff^2-2*V(i,j,k)*V(j,i,k);
                    delta=delta+U(j,i,k)*UDiff+V(j,i,k)*VDiff;
                end
                beta=2*beta;
                delta=-2*delta;

                %Solve for the zeros of the derivative of Equation 48.
                xVals=roots([4*epsilon;3*delta;2*gamma;beta]);
                if(~isreal(xVals))     
                    %Take the one real solution.      
                    [~,idx]=min(abs(imag(xVals)));
                    x=xVals(idx);
                else%If there are three real solutions, take the minimum
                    %one (not necessarily the one minimizing the cost
                    %function).

                    x=min(xVals);
                end

                %Update L and M
                EijR=eye(N,N);
                EijR(i,j)=x;
                EijRn=eye(N,N);
                EijRn(i,j)=-x;

                L=L*EijR;%Update L
                %Update M
                for k=1:K
                    M(:,:,k)=EijRn*M(:,:,k)*EijR;
                end
                
                %The x update is done. Proceed to the y update.
                U=real(M);
                V=imag(M);
                
                %Get the coefficients from Appendix B for y
                xi=sum(U(j,i,:).^2+V(j,i,:).^2);
                lambda=0;
                mu=0;
                nu=0;
                for k=1:K
                    UDiff=U(i,i,k)-U(j,j,k);
                    VDiff=V(i,i,k)-V(j,j,k);

                    %Coefficients for y
                    lambda=lambda+U(i,j,k)*VDiff-V(i,j,k)*UDiff;
                    mu=mu+UDiff^2+2*U(i,j,k)*U(j,i,k)+VDiff^2+2*V(i,j,k)*V(j,i,k);
                    nu=nu+U(j,i,k)*VDiff-V(j,i,k)*UDiff;
                end
                lambda=-2*lambda;
                nu=-2*nu;

                %Solve for the zeros of the derivative of Equation 49.
                yVals=roots([4*xi;3*nu;2*mu;lambda]);
                if(~isreal(yVals))
                    %Take the one real solution.           
                    [~,idx]=min(abs(imag(yVals)));
                    y=yVals(idx);
                else%If there are three real solutions, take the minimum
                    %one (not necessarily the one minimizing the cost
                    %function).

                    y=min(yVals);
                end
                
                EijI=eye(N,N);
                EijI(i,j)=1i*y;
                EijIn=eye(N,N);
                EijIn(i,j)=-1i*y;

                L=L*EijI;%Update L
                %Update M
                for k=1:K
                    M(:,:,k)=EijIn*M(:,:,k)*EijI;
                end
            end
        end
        
        %Compute the cost, the sum of the squared magnitudes of the sub-
        %diagonal terms, then check for convergence.
        curCost=sum(MSubEls(:).*abs(M(:)).^2);
        diff=abs(curCost-prevCost);
        if(diff<prevCost*eps())
            exitCode=0;%It converged.
            break
        end
        prevCost=curCost;
    end 
end

%Compute V as in Section IV-B. The M matrices already hold what is called R
%in the text.
V=eye(N,N);%The diagonal is known to be 1.
aij=zeros(K,1);
bij=zeros(K,1);
for j=2:N
    for i=(j-1):-1:1 
        %Equation 25
        for k=1:K
            aij(k)=M(j,j,k)-M(i,i,k);

            bij(k)=0;
            for p=(i+1):j
                bij(k)=bij(k)+M(i,p,k)*V(p,j);
            end
        end

        V(i,j)=aij'*bij/norm(aij)^2;
    end
end

if(nargout>1)
    D=zeros(N,N,K);
    for k=1:K
       D(:,:,k)=diag(diag(M(:,:,k))); 
    end
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
