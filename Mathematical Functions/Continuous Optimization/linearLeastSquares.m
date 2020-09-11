function x=linearLeastSquares(A,b)
%%LINEARLEASTSQUARES Solve the problem arg min_x norm(A*x-b) where A can
%           be square or rectangular (overdetermined or underdetermined),
%           but is assumed to be full rank. If A is rank deficient, then
%           one might want to use linearLeastSquaresRankDef.
%
%INPUTS: A A real mXn matrix.
%        b An mXnumCol vector.
%
%OUTPUTS: x The nXnumCol minimum l2 norm solution for each column of b.
%
%For the overdetermined case, Algorithm 5.3.2 of Chapter 5.3 of [1] is
%used. For the underdetermined case, Algorithm 5.6.1 in Chapter 5.6
%of [1] is used. The the square case, a QR decomposition with
%backsubstitution is used.
%
%EXAMPLE:
%The cost value can be driven to zero in the underdetermined case even
%if the measurements are corrupted with noise.
% n=200;
% m=190;
% x=randn(n,1);
% A=randn(m,n);
% b=A*x+randn(m,1);
% xLS=linearLeastSquares(A,b);
% costFun=norm(A*xLS-b)
%The cost function will be zero, within finite precision limits.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[m,n]=size(A);

if(n>m)
    %The underdetermined case. Algorithm 5.6.2. in [1].
    [Q,R]=qr(A');
    z=R(1:m,1:m)'\b;
    x=Q(:,1:m)*z;
elseif(n==m)
    [Q,R]=qr(A);
    y=Q'*b;
    %Solve taking advantage of the upper-triangular structure.
    %(Back substitution)
    opts.UT=true;
    opts.LT=false;
    opts.RECT=false;
    x=linsolve(R,y,opts);
else%Algorithm 5.3.2 in [1] for the overdetermined case.
    [A,R]=HouseholderQR(A);
    
    numCol=size(b,2);
    for j=1:n
        v=[1;A((j+1):m,j)];
        betaVal=2/(v'*v);
        for curCol=1:numCol
            b(j:m,curCol)=b(j:m,curCol)-betaVal*(v'*b(j:m,curCol))*v;
        end
    end

    %Solve taking advantage of the upper-triangular structure.
    %(Back substitution)
    opts.UT=true;
    opts.LT=false;
    opts.RECT=false;
    x=linsolve(R(1:n,1:n),b(1:n,:),opts);
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
