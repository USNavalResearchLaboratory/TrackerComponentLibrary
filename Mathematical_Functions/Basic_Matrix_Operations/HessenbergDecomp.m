function [H,U0]=HessenbergDecomp(A)
%%HESSENBERGDECOMP Perform a Hessenberg Decomposition of the square matrix
%                  A. This is U0'*A*U0=H, where H is a Hessenberg matrix
%                  (everything below the main diagonal is zero) and
%                  U0'*U0=eye(n,n).
%
%INPUTS: A An nXn matrix.
%
%OUTPUTS: H An nXn Hessenberg matrix.
%        U0 An nXn matrix such that U0'*U0 is the identity matrix and
%           U0'*A*U0=H.
%
%This function implements Algorithm 7.4.2 in Chapter 7.4.3 of  [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

if(nargout<2)
    for k=1:(n-2)
        [v,beta]=HouseholderVec(A((k+1):n,k));
        A((k+1):n,k:n)=A((k+1):n,k:n)-(beta*v)*(v'*A((k+1):n,k:n));
        A(1:n,(k+1):n)=A(1:n,(k+1):n)-beta*(A(1:n,(k+1):n)*v)*v';
    end
    %The lower-triangular part of A is already zero, except for finite
    %precision limits. This just forces it to 0.
    H=triu(A,-1);
else
    %We also want U0.
    for k=1:(n-2)
        [v,beta]=HouseholderVec(A((k+1):n,k));
        A((k+1):n,k:n)=A((k+1):n,k:n)-(beta*v)*(v'*A((k+1):n,k:n));
        A(1:n,(k+1):n)=A(1:n,(k+1):n)-beta*(A(1:n,(k+1):n)*v)*v';
        %Store the Householder vector in the part of A that will be all zeros.
        A((k+2):n,k)=v(2:(n-k));
    end
    %Extract the Hessenberg matrix.
    H=triu(A,-1);
    
    U0=zeros(n,n);
    U0(1,1)=1;
    U0(2:n,2:n)=factFormMat2Explicit(A(2:n,1:(n-1)));
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
