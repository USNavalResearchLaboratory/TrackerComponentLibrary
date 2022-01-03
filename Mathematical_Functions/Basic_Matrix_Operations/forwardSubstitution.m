function b=forwardSubstitution(L,b,algorithm)
%%FORWARDSUBSTITUTION Given a lower-triangular matrix, use forward
%       substitution to solve the system L*x=b for x. The mex version of
%       this is only implemented for real matrices.
%
%INPUTS: L An nXn positive-definite lower-triangular matrix.
%        b An nX1 vector.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            row-oriented algorithm, Algorithm 3.1.1 in [1].
%          1 Use the column-oriented algorithm, Algorithm 3.1.3 in [1].
%
%OUTPUTS: b The x in the solution to L*x=b.
%
%EXAMPLE:
%In this example, we generate a random lower-triangular system and solve
%it using the two algorithms and with whatever Matlab uses as its default
%with the \ command. Often, algorithm 0 is the best and outperforms
%Matlab's default.
% numMC=1e5;
% n=20;
% err0=0;
% err1=0;
% errMat=0;
% for k=1:numMC
%     L=tril(randn(n,n));
%     x=(1:n).';
%     b=L*x;
% 
%     x0=forwardSubstitution(L,b,0);
%     x1=forwardSubstitution(L,b,1);
%     xM=L\b;
% 
%     err0=err0+max(abs(x0(:)-x(:)));
%     err1=err1+max(abs(x1(:)-x(:)));
%     errMat=errMat+max(abs(xM(:)-x(:)));
% end
% err0=err0/numMC
% err1=err1/numMC
% errMat=errMat/numMC
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
	algorithm=0;
end

n=length(b);
switch(algorithm)
    case 0%Algorithm 3.1.1 in [1], row-oriented.
        n=length(b);
        b(1)=b(1)/L(1,1);
        for k=2:n
            b(k)=(b(k)-L(k,1:(k-1))*b(1:(k-1)))/L(k,k);
        end
    case 1%Algorithm 3.1.3 in [1], column-oriented.
        for k=1:(n-1)
            b(k)=b(k)/L(k,k);
            b((k+1):n)=b((k+1):n)-L((k+1):n,k)*b(k);
        end
        b(n)=b(n)/L(n,n);
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
