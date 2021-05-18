function [Q,Z,QAZ,QBZ]=HessenbergTriRed(A,B)
%%HESSENBERGTRIRED Perform a Hessenberg triangular reduction. Given A and
%             B, two nXn matrices, find orthogonal matrices Q and Z such
%             that Q'*A*Z is an upper Hessenberg matrix and Q'*B*Z is an
%             upper triangular matrix.  An upper Hessenberg matrix has zero
%             entries below the first subdisgonal and is thus almost upper
%             triangular. This implementation is only for real matrices.
%
%INPUTS: A, B Two nXn real matrices.
%
%OUTPUTS Q,Z Orthogonal nXn matrices such that Q'*A*Z is upper Hessenberg
%            and Q'*B*Z is lower triangular.
%        QAZ The matrix Q'*A*Z.
%        QBZThe matrix Q'*B*Z
%
%This implements algorithm 7.7.1 in Chapter 7.7.4 of [1], modified to
%obtain the Q and Z matrices explicitly. This is similar to the hess
%function in Matlab.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

%The initial upper triangularization. We are just using qr and not
%Algorithm 5.2.1.
%Make B upper triangular
[Q,B]=qr(B,0);
A=Q'*A;%Adjust A accordingly

Z=eye(n,n);
for j=1:(n-2)
    for i=n:-1:(j+2)
        %Zero A(i,j)
        [c,s]=GivensCS(A(i-1,j),A(i,j));
        A((i-1):i,j:n)=[c, s;-s,c]'*A(i-1:i,j:n);
        B((i-1):i,(i-1):n)=[c, s;-s,c]'*B(i-1:i,i-1:n);
        %Store the transformation
        Q(:,(i-1):i)=Q(:,(i-1):i)*[c,s;-s,c];
        %Zero B(i,j)
        [c,s]=GivensCS(-B(i,i),B(i,i-1));
        B(1:i,i-1:i) = B(1:i,i-1:i)*[c,s;-s,c];
        A(1:n,i-1:i) = A(1:n,i-1:i)*[c,s;-s,c];
        Z(:,(i-1):i) = Z(:,(i-1):i)*[c,s;-s,c];
    end
end

%Save the results.
QAZ=A;
QBZ=B;
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
