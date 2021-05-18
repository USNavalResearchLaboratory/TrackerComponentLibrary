function C=multiDimPolyMatMultiply(A,B)
%%MULTIDIMPOLYMATMULTIPLY Given 2D matrices of multivariate polynomials,
%          multiply the matrices to get a resulting matrix containing
%          multivariate polynomials. The multivariate polynomials
%          themselves are represented as hypermatrices, so the "matrices"
%          of these polynomials are represented by cell arrays of
%          hypermatrices holding the polynomial coefficients.
%
%INPUTS:   A An mXr cell array where each element contains a hypermatrix
%            representing a multivariate polynomials. The format of these
%            hypermatrices is discussed below.
%          B An rXn cell array where each element contains a hypermatrix
%            representing a multivariate polynomial.
%
%OUTPUTS: C The mXn-sized product of the matrices of multivariate
%           polynomials, A and B. C is a cell array where each entry is a
%           hypermatrix of coefficients of the multivariate polynomial.
%
%The coefficients of a hypermatrix representing a multivariate polynomial
%coeffs are arranged such that coeffs(a1,a2,a3...an) corresponds to the
%coefficient of an x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus,
%the number of indices coeffs takes is equal to the dimensionality of x
%(not counting singleton dimensions at the end of coeffs). Note that this
%ordering is the reverse that used in the 1D polyval function that is built
%into Matlab. The number of elements for each index in coeffs is the
%maximum order of that dimension +1. The indices of the multivariate
%polynomials in the different elements of the cell arrays must represent
%the same variable.
%
%This just implements the rules for basic matrix multiplication as in
%Chapter 1.11 of [1], using the appropriate functions for multiplying
%(convn) and summing (polySumMultiDim) multivariate polynomials.
%
%EXAMPLE:
%COnsider multiplying a 2X2 matrix of polynomials by a 2X1 vector of
%polynomials.
% A=cell(2,2);
% B=cell(2,1);
% %A multivariate linear polynomial with 3 variables.
% %-3+4*x1+5*x1*x2*x3
% A11=zeros(2,2,2);
% A11(0+1,0+1,0+1)=-3;
% A11(1+1,0+1,0+1)=4;
% A11(1+1,1+1,1+1)=5;
% A{1,1}=A11;
% %Both diagonals equal and just a single third-order term.
% %12*x1^3
% A12=zeros(4,1,1);
% A12(3+1,0+1,0+1)=12;
% A{1,2}=A12;
% A{2,1}=A12;
% %A second order polynomial only with x1 and x3.
% %6+x1^2-2*x3^2
% A22=zeros(3,1,3);
% A22(0+1,0+1,0+1)=6;
% A22(2+1,0+1,0+1)=1;
% A22(0+1,0+1,2+1)=-2;
% A{2,2}=A22;
% %Now the B vector.
% %-6+x1*x2+x1^2
% B1=zeros(2,3);
% B1(0+1,0+1)=-6;
% B1(1+1,1+1)=1;
% B1(2+1,0+1)=1;
% B{1}=B1;
% %-12-4*x1^2*x2^2
% B2=zeros(3,3);
% B2(0+1,0+1)=-12;
% B2(2+1,2+1)=-4;
% B{2}=B2;
% C=multiDimPolyMatMultiply(A,B)
% %One finds that C is a 2X1 cell array where the highest exponent in C{1} is
% %3 and the highest in C{2} is 4. We can evaluate the polynomial at a
% %point. Say
% x=[1;2;3];
% f=@(coeffs)polyValMultiDim(coeffs,x);
% val=cellfun(f,C)
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
r=size(A,2);
n=size(B,2);

C=num2cell(zeros(m,n));

for j=1:n
    for k=1:r
        for i=1:m
           polyProd=convn(A{i,k},B{k,j});
           C{i,j}=polySumMultiDim(C{i,j},polyProd);
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
