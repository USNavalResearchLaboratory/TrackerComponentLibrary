function sqrtA=cholSemiDef(A,upperOrLower,method,epsVal)
%%CHOLSEMIDEF Perform a Cholesky decomposition on a real, symmetric matrix
%             that is positive definite, positive semi-definite, or
%             possibly very slightly non-positive semi-definite due to
%             finite precision errors. Matlab's implementation of the
%             Cholesky decomposition cannot handle positive semi-definite
%             matrices. Multiple approaches to this problem are
%             implemented.
%
%INPUTS: A A symmetric, real positive semi-definite (or nearly positive
%          semi-definite) matrix.
% upperOrLower A string indicating whether an upper-triangular or lower-
%          triangular Cholesky decomposition is desired. If omitted or an
%          empty matrix is passed, the default value of 'upper' is used.
%          This can take the values 'upper' and 'lower'.
%   method This selects the algorithm used to perform the decomposition.
%          The techniques based on an LDL' decomposition are usually the
%          best. Possible values are
%          0 (The default if this parameter is omitted or an empty matrix
%            is passed) Return a lower-triangular matrix such that
%            A=sqrtA*sqrtA'. This method performs an LDL decomposition
%            using Matlab's ldl function followed by a triangularization
%            step to assure that the block-diagonal result is actually
%            diagonal. Chapter 4.2.3 of [1] related the LDL decomposition
%            to the Cholesky decomposition. values in the diagonal matrix D
%            that are less than epsVal*max(diag(D)) are replaced with
%            epsVal*max(diag(D)).
%          1 This is similar to 0, but no triangularization step is
%            performed. Thus, A=sqrtA*sqrtA' but sqrtA might not be
%            completely diagonal. This is faster than option 0.
%          2 Perform a traditional Cholesky decomposition as in Chapter 4.2
%            of [1]. However, if the argument of the square root step in
%            the algorithm is less than epsVal, replace it with epsVal. 
%    epsVal A positive value that affects how the algorithm works with
%          semi-definite and non-positive definite matrices. The default
%          for methos 0 and 1 is 0. The default for method 2 is eps().
%          Using nonzero values in methods 0 and 1 will not guarantee that 
%          A=sqrtA*sqrtA'. Method 2 requires this to be nonzero.   
%
%OUTPUTS: sqrtA The upper or lower-triangular Cholesky decomposition of A.
%               If 'upper' is chosen, then sqrtA'*sqrtA=A. Otherwise
%               sqrtA*sqrtA'=A.
%
%A Cholesky decomposition finds a triangular (or for method 1, almost
%triangular) matrix C such that A=C*C'. A solution only exists if A is
%symmetric and positive semi-definite. Given a non-positive semi-definite
%matrix, this function will return a matrix that is "close" to satisfying
%the equality.
%
%This function uses the fact that Matlab's LDL decomposition algorithm can
%handle positive semi-definite matrices. The algebraic relationship between
%an LDL decomposition and a Cholesky decomposition is then used. More
%information on the LDL decomposition of singular and ill-conditioned
%matrices is given in Chapter 4.2.7 of [1], and Chapter 4.2.3 relates the
%LDL decomposition to the Cholesky decomposition.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(upperOrLower))
    upperOrLower='upper';
end

if(nargin<3||isempty(method))
   method=0; 
end

if(nargin<4||isempty(epsVal))
    if(method==2)
        epsVal=eps();
    else
        epsVal=0;
    end
end

if(method==0||method==1)
    [L,D]=ldl(A);

    %The absolute value is to handle precision limitations that would make an
    %otherwise positive semi-definite matrix not positive semi-definite. The
    %tria function is to deal with the fact that when A is poorly conditioned,
    %Matlab will return an L that is not actually lower-triangular. Thus, while
    %L*sqrt(D) is indeed a matrix square root, it is not a valid Cholesky
    %decomposition. The tria function turns it into a valid Cholesky
    %decomposition.
    dAbs=abs(diag(D));
    diagVal=diag(sqrt(max(diag(D),epsVal*max(dAbs))));

    sqrtA=L*diagVal;
    if(method==0)
        sqrtA=tria(L*diagVal);
    end
elseif(method==2)
    sqrtA=cholTradLoad(A,epsVal);
else
    error('Unknown method specified.')
end

switch(upperOrLower)
    case 'lower'
    case 'upper'
        sqrtA=sqrtA';
    otherwise
        error('An invalid value was entered for the upperOrLower parameter.')
end
end

function L=cholTradLoad(A,epsVal)
n=length(A);
L=zeros(n,n);
for j=1:n
    temp=A(j,j)-sum(L(j,1:(j-1)).*conj(L(j,1:(j-1))));
    
    temp=max(epsVal,temp);
    
    L(j,j)=sqrt(temp);
    for i=(j+1):n
        L(i,j)=(A(i,j)-sum(L(i,1:(j-1)).*conj(L(j,1:(j-1)))))/L(j,j);
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