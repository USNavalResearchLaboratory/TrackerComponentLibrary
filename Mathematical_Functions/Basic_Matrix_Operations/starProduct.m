function C=starProduct(A,B)
%%STARPRODUCT Take the star product of two matrices. The star product of an
%             mXn matrix A and an (m*p)X(n*q) matrix B is sum_{i,j}
%             a_{i,j}B_{i,j} where a_{i,j} is the ijth subelement of A and
%             B_{ij} is the ijth mXn-dimensional submatrix of B.
%
%INPUTS: A An mXn matrix.
%        B An (m*p)X(n*q) matrix B.
%
%OUTPUTS: C The pXq star matrix product of A and B.
%
%Various properties of the star matrix product are given in [1].
%
%EXAMPLE:
%The start product of A and B, if they are both the same size, is just
%trace(A'*B). We demonstrate that here. diffVal will be zero.
% m=3;
% n=2;
% A=eye(m,n);
% B=rand(m,n);
% diffVal=starProduct(A,B)-trace(A'*B)
%
%REFERENCES:
%[1] E. C. MacRae, "Matrix derivatives with an application to an adaptive
%    linear decision problem," The Annals of Statistics, vol. 2, no. 2, pp.
%    337-346, Mar. 1974.
%
%February 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
n=size(A,2);

p=size(B,1)/m;
q=size(B,2)/n;

if(fix(p)~=p||fix(q)~=q)
    error('The dimensions of B are not consistent with the dimensions of A.')
end

C=zeros(p,q);
sel1=1:p;
for i=1:m
    sel2=1:q;
    for j=1:n
        C=C+A(i,j)*B(sel1,sel2);
        sel2=sel2+q;
    end
    sel1=sel1+p;
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
