function C=matMultiply(A,B,typename)
%MATMULTIPLY Multiply two matrices. This is the same as A*B in Matlab. The
%            point of this function is so that classes that override basic
%            addition and multiplication operations can use this for
%            multiplication.
%
%INPUTS: A An mXr matrix. This must be of a type where the size function
%          works.
%        B An rXn matrix. This must be of a type where the size function
%          works.
% typename An optional string for how space for the output C should be
%          allocated. The string is used in feval(typename,m,n); For
%          example, if C should be an instance of a class called Interval,
%          then typename='Interval' to call the constructor of the class.
%          If omitted or an empty matrix is passed, then C is allocated as
%          zeros(m,n).
%
%OUTPUTS: C The mXn-sized product of matrices A and B.
%
%Basic matrix multiplication as given in Chapter 1.1.11 of [1] is
%performed. The .* operator is used for multiplying the individual elements
%of A and B, allowing this to be used in a routine for overloading a *
%operator.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
r=size(A,2);
n=size(B,2);

%Allocate space for the result. The optional typename lets this function be
%used with classes that override scalar addition and multiplication
%operations, as might be the case in a class used for interval analysis.
if(nargin>2&&~isempty(typename))
    C=feval(typename,m,n);
else
    C=zeros(m,n);
end

for j=1:n
    for k=1:r
        for i=1:m
          C(i,j)=C(i,j)+A(i,k).*B(k,j); 
        end
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
