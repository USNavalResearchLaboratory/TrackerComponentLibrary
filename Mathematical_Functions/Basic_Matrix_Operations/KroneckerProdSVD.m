function [sigma,B,C]=KroneckerProdSVD(A,m1,n1,m2,n2)
%%KRONECKERPRODSVD Perform a Kronecker product singular value
%         decomposition. Given the m1*m2Xn1*n2 matrix A that is an m1Xn1
%         matrix composed of blocks of size m2Xn2, this finds sigma, B and
%         C such that A=sum_{k=1}^r(sigma(k)*kron(B(:,:,k),C(:,:,k)),
%         where r is the maximum possible rank of the Kronecker product,
%         r=min(m1*n1,m2*n2). If the first rTilde<r terms in the sum are
%         used, then the result is the closest matrix to A (with regard to
%         the Frobenius norm) that is the sum of rTilde Kronecker products.
%         This function is only implemented for real matrices.
%
%INPUTS: A A real m1*m2Xn1*n2 block matrix consisting of m1Xn1 blocks of
%          size m2Xn2. A must be real.
% m1,n1,m2,n2 The integer sizes making up the dimensions of A.
%
%OUTPUTS: sigma An rX1 vector of the Kronecker singular values, 
%               r=min(m1*n1,m2*n2).
%             B An m1Xn1Xr set of r m1Xn1 matrices.
%             C An m2Xn2Xr set of r matrices.
%
%The Kronecker product singular value decomposition is described in Chapter
%12.3.6 of [1]. The algorithm is taken from Theorem 12.3.1 of that chapter.
%
%EXAMPLE 1:
% B=magic(3);
% m1=size(B,1);
% n1=size(B,2);
% C=[magic(4),ones(4,1)];
% m2=size(C,1);
% n2=size(C,2);
% A=kron(B,C);
% [sigma,BKSVD,CKSVD]=KroneckerProdSVD(A,m1,n1,m2,n2);
% AEst=sigma(1)*kron(BKSVD(:,:,1),CKSVD(:,:,1));
% AEst-A
%One will observe that only sigma(1) is meaningfully above zero. That is
%because A has Kronecker rank-1. Thus, using only the first term, one can
%reconstruct A within expected finite precision limitations. However,
%BKSVD(:,:,1) does not equal B and CKSVD(:,:,1) does not equal C.

%EXAMPLE:
%Here, we show that A can be reconstructed using U, sigma V and Kronecker
%products.
% m1=2;
% m2=3;
% n1=3;
% n2=4;
% A=randn(m1*m2,n1*n2);
% [sigma,U,V]=KroneckerProdSVD(A,m1,m2,n1,n2);
% r=min(m1*n1,m2*n2);
% AK=zeros(m1*m2,n1*n2);
% for k=1:r
%     AK=AK+sigma(k)*kron(U(:,:,k),V(:,:,k));
% end
% RelErr=max(abs((AK(:)-A(:))./A(:)))
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(A))
    sigma=[];
    B=[];
    C=[];
    return;
end

if(size(A,1)~=m1*m2||size(A,2)~=n1*n2)
    error('A has inconsistent dimensions')
end

if(any(~isreal(A(:))))
    error('This function is only implemented for real matrices')
end

RA=makeRA(A,m1,n1,m2,n2);

%The SVD as in Equation 12.3.17
[U,Sigma,V]=svd(RA);
sigma=diag(Sigma);

%The maximum possible Kronecker product rank
r=min(m1*n1,m2*n2);

B=zeros(m1,n1,r);
C=zeros(m2,n2,r);

for n=1:r
    B(:,:,n)=reshape(U(:,n),m1,n1);
    C(:,:,n)=reshape(V(:,n),m2,n2);
end
end

function RA = makeRA(A,m1,n1,m2,n2)
%%MAKERA Let A be an m1Xn1 block matrix where each block has size m2Xn2.
%        This function compute RA below Equation 12.3.16 in Chapter 12.3.6
%        of [1].

RA =zeros(m1*n1,m2*n2);
curRow=1;
for j=1:n1
    for i=1:m1
        Aij=A(1+(i-1)*m2:i*m2,1+(j-1)*n2:j*n2);
        RA(curRow,:)= Aij(:).';
        curRow=curRow+1;
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
