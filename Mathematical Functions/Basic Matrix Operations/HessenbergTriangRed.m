function [uHessMat,uTriMat,Q,Z]=HessenbergTriangRed(A,B)
%%HESSENBERGTRIANGRED Perform a Hessenberg triangular reduction of the nXn
%               matrices A and B. This finds Q and Z such that Q'*A*Z= an
%               upper Hessenberg matrix and Q'*B*Z is an upper triangular
%               matrix.
%
%INPUTS: A, B Two real nXn matrices.
%
%OUTPUTS: uHessMat The upper Hessenberg matrix obtained from Q'*A*Z.
%          uTriMat The upper triangular matrix obtained with Q'*B*Z.
%
%This function implements algorithm 7.7.1 in [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

if(nargout>2)
    Z=eye(n,n);
    %B=Q*R. Replace B with Q'*B and A with Q'*A
    [Q,R]=qr(B);
    B=R;%Q'*B=R.
    A=Q'*A;

    for j=1:(n-2)
        for i=n:-1:(j+2)
            %Zeros A(i,j)
            [c,s]=GivensCS(A(i-1,j),A(i,j));
            rotMat=[c, s;
                   -s, c];
            A((i-1):i,j:n)=rotMat'*A((i-1):i,j:n);
            B((i-1):i,(i-1):n)=rotMat'*B((i-1):i,(i-1):n);

            Q(:,(i-1):i) = Q(:,(i-1):i)*rotMat;

            %Zeros B(i,i-1).
            [c,s]=GivensCS(-B(i,i),B(i,i-1));
            rotMat=[c, s;
                   -s, c];
            B(1:i,(i-1):i) = B(1:i,(i-1):i)*rotMat;
            A(1:n,(i-1):i) = A(1:n,(i-1):i)*rotMat;

            Z(:,(i-1):i) = Z(:,(i-1):i)*rotMat;
        end
    end
    uHessMat=A;
    uTriMat=B;
else
    %If Q and Z are not desired, don't compute them.
    
    %B=Q*R. Replace B with Q'*B and A with Q'*A
    [Q,R]=qr(B);
    B=R;%Q'*B=R.
    A=Q'*A;
    for j=1:(n-2)
        for i=n:-1:(j+2)
            %Zeros A(i,j)
            [c,s]=GivensCS(A(i-1,j),A(i,j));
            rotMat=[c, s;
                   -s, c];
            A((i-1):i,j:n)=rotMat'*A((i-1):i,j:n);
            B((i-1):i,(i-1):n)=rotMat'*B((i-1):i,(i-1):n);

            %Zeros B(i,i-1).
            [c,s]=GivensCS(-B(i,i),B(i,i-1));
            rotMat=[c, s;
                   -s, c];
            B(1:i,(i-1):i) = B(1:i,(i-1):i)*rotMat;
            A(1:n,(i-1):i) = A(1:n,(i-1):i)*rotMat;
        end
    end
    uHessMat=A;
    uTriMat=B;
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
