function [x,M]=linearLeastSquaresRankDef(A,b,algorithm,rankAlg)
%%LINEARLEASTSQUARESRANKDEF Solve the problem arg min_x norm(A*x-b) where A
%           can be square or rectangular (overdetermined or
%           underdetermined), and can be rank deficient.
%
%INPUTS: A A real mXn matrix. This can (but does not ahve to be) be rank
%          deficient.
%        b An mXnumCol vector.
% algorithm An optional parameter specifying which solution should be
%         taken (rank-deficient matrices have an infinite number of
%         solutions). Possible values are:
%         0 (The default if omitted or an empty matrix is passed) Get the
%           basic solution using Algorithm 5.1.1 in Chapter 5.5.7 of [1].
%           If the matrix A is rank r, then the solution will have at least
%           n-r zeros.
%         1 Obtain the minimum-norm solution. This is given by the
%           algorithm implcit in theorem 5.5.1 in Chapter 5.5.1 of [1].
%
%OUTPUTS: x The chosen nXnumCol minimum l2 norm solutions, one for each
%           column of b.
%         M If x is the minimum norm solution, then one can obtain any
%           other solution using x+M*y, where y is any vector.
%
%To obtain M, one considers that norm(A*x-b)^2 is the same cost function as
%an unconstrained quadratic programming problem. Consequently, the relation
%between the minimum norm solution and any other solution in Equation 3.2
%of [2] can be used to get M.
%
%EXAMPLE: 
%In this example, we have a rank deficient matrix A and we choose x to have
%one zero. Solving the problem for the basic solution produces the original
%x (because it is "sparse"). Solving for the minimum norm solution produces
%a different solution. One sees that both solutions solve the equality
%within finite precision limits and that the norm of the second solution is
%indeed less than the norm of the first. If one were to change the zero in
%x to something else, then one can expect that neither solutions equal the
%original x.
% A=magic(4);
% xTrue=[1;0;2;3];
% b=A*xTrue;
% x0=linearLeastSquaresRankDef(A,b,0)%Basic solution
% x1=linearLeastSquaresRankDef(A,b,1)%Minimum norm solution.
% 
% %Both solutions solve the equality.
% A*x0-b
% A*x1-b
% 
% %x1 has a smaller norm than x0.
% norm(x0)
% norm(x1)
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%[2] D. L. Nelson, "Quadratic programming techniques using matrix
%    pseudoinverses," Ph.D. dissertation, Texas Tech University, Lubbock,
%    TX, Dec. 1969.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

if(nargin<4)
    rankAlg=[];
end

n=size(A,2);
numCol=size(b,2);
[rankVal,V,U,S]=matrixRank(A,rankAlg);

if(rankVal==n)
    %It is a full rank problem.
    x=linearLeastSquares(A,b);
    return
end

switch(algorithm)
    case 0%The basic solution; Algorithm 5.5.1 in [1].
        [~,~,P]=qr(V(:,1:rankVal)');

        AP=A*P;
        B1=AP(:,1:rankVal);

        z=linearLeastSquares(B1,b);

        x=P*[z;zeros(n-rankVal,numCol)];
    case 1%The minimum norm solution.
        s=diag(S);
        
        x=zeros(n,numCol);
        for curCol=1:numCol
            for k=1:rankVal
                x(:,curCol)=x(:,curCol)+(U(:,k)'*b(:,curCol)/s(k))*V(:,k);
            end
        end
        
        %The above is equivalent to:
        % T11=S(1:rankVal,1:rankVal);
        % 
        % cd=U'*b;
        % c=cd(1:rankVal);
        % 
        % x=V*[T11\c;zeros(n-rankVal,1)];
    otherwise
        error('Unknown algorithm specified.')
end

if(nargout>1)
    B=A'*A;
    pinvB=pseudoInverse(B);
    %As in Equation 3.2 in [2].
    M=eye(n,n)-pinvB*B;
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
