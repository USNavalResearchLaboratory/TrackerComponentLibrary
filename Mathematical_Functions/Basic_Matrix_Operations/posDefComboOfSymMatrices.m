function [B,alpha,exitCode]=posDefComboOfSymMatrices(A,alphaInit,maxIter)
%%POSITIVEDEFCOMBOOFSYMMATRICES Given K real, symmetric matrices, the kth
%           one being A(:,:,k), find weights alpha such that
%           B=sum_{k=0}^K alpha(k)*A(:,:,k) is positive definite. Note that
%           it is possible to provide a set of A matrices such that no
%           solution exists. If a solution exists, other solutions exist
%           within a scale factor. However, there are often multiple
%           solutions even when not considering just a scale factor. This
%           function simply returns one feasible solution without regard
%           for an optimality criterium.
%
%INPUTS: A An NXNXnumMats set of real, symmetric matrices.
% alphaInit This is an optional parameter. it is an initial estimate for
%           alpha. If this is omitted or an empty matrix is passed, then a
%           alpha=randn(numMats,1) is used.
%   maxIter The maximum number of iterations to perform. If no solution
%           exists, then the function will iterate this many times. The
%           default if omitted or an empty matrix is passed is 200.
%
%OUTPUTS: B The NXN real, symmetric matrix formed as a linear combination
%           of the A matrices. If a solution was found, this matrix is
%           positive definite.
%     alpha The numMatsX1 weights used to form B. These can be positive and
%           negative.
%  exitCode This indicates how the algorithm terminated. Possible values
%           are:
%           0 The algorithm converged within maxIter iterations and B
%             should be positive definite.
%           1 The algorithm did not convergence. One can assume that B is
%             not positive definite.
%
%This function implements the algorithm of [1].
%
%EXAMPLE 1:
%Here is an example where a solution exists, but it is not unique. All of
%the A matrices have negative eigenvalues, but this function finds a
%weighting to create a positive definite solution. Here, we find two
%solutions, display their eigenvalues (which will all be positive) and show
%that the different weight vectors are not just multiples of each other (if
%they were, the ratio of the weights for the two solutions would be the
%same for all elements of alpha).
% A=zeros(3,3,3);
% A(:,:,1)=[-0.0904,    0.2294,    0.1721;
%            0.2294,   -0.5820,   -0.4367;
%            0.1721,   -0.4367,   -0.3276];
% A(:,:,2)=[2.5574,    1.7122,    0.8704;
%           1.7122,    1.1463,    0.5828;
%           0.8704,    0.5828,    0.2963];
% A(:,:,3)=[0.4337,    0.8187,    0.6504;
%           0.8187,   -0.3554,    0.9081;
%           0.6504,    0.9081,    0.9217];
% [B,alpha1]=posDefComboOfSymMatrices(A);
% eigVals1=eig(B)%All positive
% [B,alpha2]=posDefComboOfSymMatrices(A);
% eigVals2=eig(B)%All positive.
% alpha1./alpha2
%
%EXAMPLE 2:
%This is an example where no possible solution exists. The function doesn't
%identify infeasibility, but it will fail to solve the problem, thus
%returning exitCode=1.
% A=zeros(2,2,2);
% A(:,:,1)=[-1,2;
%            2,-1];
% A(:,:,2)=[1,0;
%           0,-1];
% [B,alpha,exitCode]=posDefComboOfSymMatrices(A);
% eig(B)%One can see that B has a negative eigenvalue.
%
%REFERENCES:
%[1] A. Zaidi, "Positive definite combination of symmetric matrices," IEEE
%    Transactions on Signal Processing, vol. 53, no. 11, pp. 4412-4416,
%    Nov. 2005.
%
%November 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(maxIter))
    maxIter=200;
end

p=size(A,1);
K=size(A,3);

%Create the A projection matrix described in Section IV-C.
AProj=pinv(reshape(A,[p^2,K]));

%Initialization, Step 1.
sigma=1;

%Initialization, Step 2.
if(nargin<3||isempty(alphaInit))
    alpha=randn(K,1);
else
    alpha=alphaInit;
end

%Initialization, step 3.
B=alpha(1)*A(:,:,1);
for k=2:K
    B=B+alpha(k)*A(:,:,k);
end

%Initialization, step 4.
lambda=eig(B);
mu=min(abs(lambda(lambda~=0)));
if(length(lambda(lambda>=0))>=length(lambda(lambda<0)))
    epsilon=1;
else
    epsilon=-1;
end
B=(epsilon/mu)*B;

exitCode=1;
for curIter=1:maxIter
    [U,D]=eig(B);
    lambda=diag(D);
    lambdap=min(lambda);
    if(lambdap>0)
        %The algorithm converged.
        exitCode=0;
        return;
    end

    delta=lambda-sigma;
    delta(lambda>sigma)=0;
    Delta=diag(delta);
    %The paper has D+Delta.
    alpha=AProj*vec((U*(D-Delta)*U'));
    B=alpha(1)*A(:,:,1);
    for k=2:K
        B=B+alpha(k)*A(:,:,k);
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
