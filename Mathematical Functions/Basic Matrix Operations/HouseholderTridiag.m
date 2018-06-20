function [T,QFactored,Q]=HouseholderTridiag(A)
%%HOUSEHOLDERTRIDIAG Given a real symmetric nXn matrix A, find a matrix Q
%           such that Q'*A*Q is tridiagonal. That is, it has the form:
%           [a1,b1,0,  0, ... 0;
%            b1,a2,b2, 0, ... 0;
%            0, b2,a3,b3, ... 0;
%            ..................
%            0,  0, 0, 0, ... an]; 
%            Q can be given in a factored form or explicitly.
%
%INPUTS: A A real symmetric nXn matrix.
%
%OUTPUTS: T The real nXn tridiagonalization of the matrix A.
% QFactored The nXn matrix Q in factored form. The function
%           factFormMat2Explicit converts it to explicit form.
%
%This function implements algorithm 8.3.1 of Chapter 8.3.1 of [1].
%Tridiagonalization plays a role in the implementation of symmetric QR
%decompositions.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

if(~issymmetric(A))
   error('This function requires that A be symmetric.'); 
end

QFactored=zeros(n,n);
QFactored(1,1)=1;
for k=1:(n-2)
    [v,beta]=HouseholderVec(A((k+1):n,k));
    p=beta*(A((k+1):n,(k+1):n)*v);
    w=p-(beta*p'*v/2)*v;
    A(k+1,k)=norm(A((k+1):n,k));
    A(k,k+1)=A(k+1,k);

    %Store the Householder vector. This could be stored in the lower half
    %of A, since it reuses one of the elements of the tridiagonal.
    QFactored((k+1):n,k+1)=[A(k+1,k);v(2:(n-k))];
    
    %To properly overwrite A with T, one must zero the values above and
    %below the tridiagonal.
    A((k+2):n,k)=0;
    A(k,(k+2):n)=0;
    
    A((k+1):n,(k+1):n)=A((k+1):n,(k+1):n)-(v*w'+w*v');
end
QFactored(n,n)=A(n-1,n);
    
T=A;
Q=factFormMat2Explicit(QFactored);

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
