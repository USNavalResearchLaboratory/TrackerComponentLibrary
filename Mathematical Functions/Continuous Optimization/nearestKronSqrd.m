function X=nearestKronSqrd(A)
%%NEARESTKRONSQRD Given an (m^2)X(m^2) matrix, determine the matrix X that
%           minimize norm(A-kron(X,X),'fro').
%
%INPUTS: A An (m^2)X(m^2) real matrix.
%
%OUTPUTS: X An mXm matrix minimizing the Frobenius norm specified above.
%
%The algorithm is given in Chapter 12.3.8 of [1].
%
%EXAMPLE:
%We shall create A explicitely from the square of a matrix A and show that
%this function will recover the same X.
% m=25;
% X=randn(m,m);
% A=kron(X,X);
% X=nearestKronSqrd(A);
% norm(A-kron(X,X),'fro')%The error will be close to 0.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=sqrt(length(A));

if(m~=fix(m))
    error('The size of the length of A is not a squared number.')
end

RA=zeros(m^2,m^2);
for i1=1:m
    for i2=1:m
        %Defined in Equation 12.3.16 in [1].
        Aij = A((1+(i2-1)*m):(i2*m),(1+(i1-1)*m):(i1*m));
        RA(i2+(i1-1)*m,:)=Aij(:).';
    end
end
[V,D]=eig(RA);

[maxD,idx]=max(diag(D));

X=sqrt(maxD)*reshape(V(:,idx),[m,m]);

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
