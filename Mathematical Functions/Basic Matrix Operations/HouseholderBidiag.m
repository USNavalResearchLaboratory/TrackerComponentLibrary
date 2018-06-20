function [B,UFactored,VFactored,U,V]=HouseholderBidiag(A)
%HOUSEHOLDERBIDIAG Given an mXn matrix with m>n, find matrices U and V such
%               that  U'*A*V=B, where B is an mXn bidiagonal matrix of the
%               form
%               [d1,f1,0, 0,...0;
%                 0,d2,f2,0,...0;
%                ................
%                 0, 0, 0, 0,..dn;
%                 0, 0, 0, 0,...0;
%                ................
%                 0, 0, 0, 0,...0];
%                U and V can be given in factored form or explicitly.
%
%INPUTS: A A real or complex mXn matrix with m>=n.
%
%OUTPUTS: B An mXN bidiagonal matrix.
% UFactored, VFactored The mXm matrix U and the nXn matrix V given
%           in factored form. The function factFormMat2Explicit convers
%           them into explicit form.
%      U, V mXm and nXn matrices U and V formed from UFactored and
%           VFactored using factFormMat2Explicit.
%
%The bidiagonalization of a matrix is explained in Chapter 5.4.8 of [1]. It
%plays a role in the implementation fo singular value decomposition
%algorithms.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
n=size(A,2);

if(m<n)
    error('The bidiagonalization algorithm requires that m>=n.')
end

for j=1:n
    [v,beta]=HouseholderVec(A(j:m,j));
    
    A(j:m,j:n)=A(j:m,j:n)-(beta*v)*(v'*A(j:m,j:n));
    A((j+1):m,j)=v(2:(m-j+1));
    
    if(j<=(n-2))
        [v,beta]=HouseholderVec(A(j,(j+1):n)');
        A(j:m,(j+1):n)=A(j:m,(j+1):n)-(A(j:m,(j+1):n)*v)*(beta*v)';
        A(j,(j+2):n) = v(2:(n-j))';
    end
end

UFactored=zeros(m,m);
UFactored(:,1:n)=tril(A,-1);

%B holds the main diagonal of A as well as one diagonal above that. For
%rectangular matrices with numRows>numCols, rows numebred after numCols are
%all zero.
B=tril(triu(A),1);

VFactored=zeros(n,n);
for j=1:(n-2)
   VFactored((j+2):n,j+1)=A(j,(j+2):n)';
end

if(nargout>3)
    U=factFormMat2Explicit(UFactored);
    if(nargout>4)
        V=factFormMat2Explicit(VFactored);
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
