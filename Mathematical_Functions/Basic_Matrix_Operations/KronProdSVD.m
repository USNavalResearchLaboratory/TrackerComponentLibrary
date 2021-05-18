function [U,sigma,V]=KronProdSVD(A,m1,m2,n1,n2)
%%KRONPRODSVD Computer a Kronecker product singular value decomposition.
%             This finds U, sigma and V such that:
%             A=sum_{k=1}^r sigma(k)*kron(U(:,:,k),V(:,:,k));
%             where r=min(m1*n1,m2*n2) if the recreation is exact.
%
%INPUTS: A An (m1*m2)X(n1*n2) real matrix which is an m1Xn1 collection of
%          block of size m2Xn2. It has the form
%          A=[A_{1,1} ... A_{1,n1};
%             ...     ...      ...;
%             A_{m1,1} ... A_{m1,n1}];
% m1,m2,n1,n2 The dimensions specifying A and its submatrices as described
%          above.
%
%OUTPUTS: U An m1Xn1Xr colection of left multiplication matrices as shown
%           above with r=min(m1*n1,m2*n2).
%     sigma The rX1 set of singular values.
%         V The m2Xn2Xr set of right-multiplication matrices.
%
%This function implements the technique discussed in Theorem 12.3.1 of
%Chapter 12.3.6 of [1].
%
%EXAMPLE:
%Here, we show that A can be reconstructed using U, sigma V and Kronecker
%products.
% m1=2;
% m2=3;
% n1=3;
% n2=4;
% A=randn(m1*m2,n1*n2);
% 
% [U,sigma,V]=KronProdSVD(A,m1,m2,n1,n2);
% r=min(m1*n1,m2*n2);
% AK=zeros(m1*m2,n1*n2);
% for k=1:r
%     AK=AK+sigma(k)*kron(U(:,:,k),V(:,:,k));
% end
% RelErr=max(abs((AK(:)-A(:))./A(:)))
%The relative error will probability be on the order of 1e-13 or better,
%which is on the order of what one would expect in terms of finite
%precision errors.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m1m2=m1*m2;
n1n2=n1*n2;

if(m1m2~=size(A,1)||n1n2~=size(A,2))
    error('The dimensions of A are inconsistent with m1, m2, n1, and n2.')
end

r=min(m1*n1,m2*n2);

if(isempty(A))
   U=[];
   sigma=[];
   V=[];
   return
end

RA=zeros(m1*n1,m2*n2);
for i1=1:n1
    for i2=1:m1
        %Defined in Equation 12.3.16 in [1].
        Aij = A((1+(i2-1)*m2):(i2*m2),(1+(i1-1)*n2):(i1*n2));
        RA(i2+(i1-1)*m1,:)=Aij(:).';
    end
end
%This returns the singular values in descending order.
[UFull,S,VFull]=svd(RA);
sigma=diag(S);

U=zeros(m1,n1,r);
V=zeros(m2,n2,r);
for k=1:r
    U(:,:,k)=reshape(UFull(:,k),[m1,n1]);
    V(:,:,k)=reshape(VFull(:,k),[m2,n2]);
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
