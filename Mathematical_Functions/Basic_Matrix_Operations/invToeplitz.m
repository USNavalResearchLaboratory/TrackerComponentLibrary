function B=invToeplitz(x)
%%INVTOEPLITZ Use Trench's algorithm to find the inverse of a real
%    symmetric positive definite matrix given by toeplitz(x). As noted in
%    [1], this algorithm scales better than traditional matrix inversion
%    for large matrices.
%
%INPUTS: x A real length n vector such that toeplitz(x) is a positive
%          definite real matrix.
%
%OUTPUTS: B The nXn matrix inv(toeplitz(x)) as computed from Trench's
%           algorithm.
%
%This function implements Algorithm 4.7.3 from [1] with a minor scaling
%modification to deal with the first element of x not being 1.
%
%EXAMPLE:
%Here, the inverse of a Topeplits matrix computed with this function is
%compared to one computed directly. The results agree with a relative error
%on the order of finite precision limita.
% x=[5;0.8;-0.7;0.5;-0.1];
% XInv1=invToeplitz(x);
% XInv2=inv(toeplitz(x));
% RelErr=max(abs((XInv1(:)-XInv2(:))./XInv1(:)))
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

scalFact=x(1);
n=length(x);
r=x(2:end)/scalFact;
y=solveYuleWalker([1;r]);
B=zeros(n,n);

gamma=1/(1+r(1:(n-1))'*y(1:(n-1)));
v=gamma*y((n-1):-1:1);%(1:(n-1))
B(1,1)=gamma;
B(1,2:n)=v((n-1):-1:1)';
for i=2:(floor((n-1)/2)+1)
    for j=i:(n-i+1)
        B(i,j)=B(i-1,j-1)+(v(n+1-j)*v(n+1-i)-v(i-1)*v(j-1))/gamma;
    end
end

%Use symmetry and persymmetry to fill in B as described in Chapter 4.7.5 of
%[1].
for i=1:n
    for j=(i+1):n
        B(j,i)=B(i,j);
    end
end
for i=1:n
    for j=1:(n-i+1)
        B(n-i+1,n-j+1)=B(i,j);
    end
end

%Scale the final matrix.
B=B/scalFact;
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
