function xn=linEqSolveFirstnRows(A,b,n)
%%LINEQSOLVEFIRSTNROWS This function solves the problem A*x=b for the first
%                      n rows of x. The matrix A CAN be singular as long as
%                      a unique solution exists for the first n rows of x.
%                      This function is useful for finding the position
%                      estimate in an information filter state even before
%                      estimates of all of the other target state
%                      components have become observable In such an
%                      instance A is the inverse covariance matrix, x is
%                      the target state and b is the information state.
%
%INPUTS: A The NXN matrix A in the equation A*x=b, where x is unknown. The
%          matrix can be singular, but the first n rows of x should be
%          observable.
%        b The NX1 column vector b in the equation A*x=b.
%        n The number of rows of x, starting from the first row and going
%          down, for which one wishes to solve in the equation A*b. n must
%          be less than or equal to the number of rows in A.
%
%OUTPUTS: xn The first n rows of the column vector of x solved from A*x=b.
%            If any of the components of xn are not finite, then the
%            problem was not completely observable. Because of how the
%            problem is solved, if the kth component is not finite, then
%            all subsequent components will not be finite, regardless of
%            whether they are observable.
%
%The linear equation is solved by performing a modified qr decomposition on
%the matrix A such that A=Q*R, where Q is an rotation matrix and R is a
%LOWER triangular matrix. A standard QR decomposition produces an upper
%triangular matrix. By flipping the rows and columns of a and then
%performing the inverse operations on Q and R, one can get a decomposition
%where R is a lower-triangular marix. One can then write R*x=Q'*b, since
%the transpose of a rotation matrix is its inverse. The first n rows of x
%can then be solved using forward substitution.
%
%The QR decomposition is in Matlab's built-in function qr. It is also
%discussed in Chapter 5.2 of [1].
%
%REFERENCES:
%[1] G. E. Golub and C. F. van Loan, Matrix Computations, 4rd ed.
%    Baltimore, MD: Johns Hopkins University Press, 2013.
%
%June 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Perform a lower-triangular qr decomposition.
[Q,R]=qr(rot90(A,2));
Q=rot90(Q,2);
R=rot90(R,2);
%Now, Q*R=A and R is lower-triangular.

b=Q'*b;

%Perform forward substitution to solve for the first n components of x.
xn=zeros(n,1);
xn(1)=b(1)/R(1,1);
for curRow=2:n
    xn(curRow)=(b(curRow)-sum(xn(1:(curRow-1))'.*R(curRow,1:(curRow-1))))/R(curRow,curRow);
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
