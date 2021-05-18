function F=expMat(A,algorithm)
%%EXPMAT Evaluate the matrix exponential of a particular matrix nXn A. This
%       is most easily defined in terms of the Taylor series expansion:
%       exp(A)=eye(n)+A+A*A/factorial(2)+A*A*A/factorial(3)+...
%       The inverse of the function is logMat.
%
%INPUTS: A An nXn real or complex matrix.
% algorithm An optional parameter specifying how to compute the matrix
%          exponential. Possible values are:
%          0 Use Algorithm 9.3.1 or Chapter 9.3.1 of [1]. This is the Padé
%            approximation method.
%          1 (The default if omitted or an empty matrix is passed) Use the
%            Schur-Parlett method from Chapter 9.3.1 of [1] method via the
%            function SchurMatFunEval.
%
%OUTPUTS: F The matrix exponential of A.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(A);

if(nargin<2||isempty(algorithm))
    algorithm=1; 
end

if(size(A,1)~=size(A,2))
   error('A must be a square matrix.'); 
end

switch(algorithm)
    case 0
        delta=eps(max(abs(A(isfinite(A)))));

        j=max(0,1+floor(log2(norm(A,'inf'))));
        A=A/2^j;

        %Determine the order of the expansion.
        q=0;
        %epsilon(q,q) is defined before Algorithm 9.3.1 in Chapoter 9.3.1 of [1]. 8 is
        %the value of epsilon(0,0).
        epsilonQQ=8;
        while(epsilonQQ>delta)
            q=q+1;
            %The book has epsilon(q,q) in the form
            %epsilonQQ=2^(3-2*q)*(factorial(q)/factorial(2*q))^2/(2*q+1)
            %However, this equivalent formulation (which can be verified with
            %Mathematica) is more efficient.
            epsilonQQ=(2^(1-6*q)*pi*(1+2*q))/gamma(3/2+q)^2;
        end
        D=eye(n,n);
        X=eye(n,n);
        N=eye(n,n);

        c=1;
        signVal=-1;
        for k=1:q
            c=c*(q-k+1)/((2*q-k+1)*k);
            X=A*X;
            N=N+c*X;
            D=D+signVal*c*X;
            signVal=-signVal;
        end

        F=D\N;
        for k=1:j
            F=F*F;
        end
    case 1
        F=SchurMatFunEval(@(x)exp(x),A);
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
