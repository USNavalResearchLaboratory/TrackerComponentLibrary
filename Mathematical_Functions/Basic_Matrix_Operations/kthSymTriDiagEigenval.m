function lambda=kthSymTriDiagEigenval(T,k)
%%KTHSYMTRIDIAGEIGENVAL This function find the kth largest eigenvalue of a
%   real symmetric tridiagonal matrix with no 0 subdiagonal values (which
%   would lead to repeated eigenvalues). This is done efficiently without
%   finding all of the other eigenvalues.
%
%INPUTS: T An nXn real symmetric tridiagonal matrix where no subdiagonal
%          terms are 0 (which would lead to repeated roots).
%        k The order of the eigenvalue to find 1<=k<=n.
%
%OUTPUTS: lambda The kth eigenvalue of T.
%
%This function implements the algorithm that is described in Section 8.4.2
%of [1].
%
%EXAMPLE:
%The maximum relative error of this function a scompared to usign the eig
%function and choosing the kthsorted eigenvalue is computed and shown to be
%on the order of what one would expect due to finite precision limitations.
% n=10;
% k=3;
% numMCRuns=1000;
% maxRelErr=0;
% for curRun=1:numMCRuns
%     a=randn(n,1);
%     b=randn(n-1,1);
%     T=makeSymTriDiagMat(a,b);
%     lambdak=kthSymTriDiagEigenval(T,k);
%     lambda=sort(eig(T),'descend');
%     RelErr=abs((lambda(3)-lambdak)./lambda(3));
%     maxRelErr=max(maxRelErr,RelErr);
% end
% maxRelErr
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(T,1);
z=norm(T,1);%The maximum absolute column sum.
y=-z;
%The tolerance is arbitrarily set larger than eps.
while(abs(y-z)>eps(abs(y)+abs(z)))
    x=(y+z)/2;
    if(numSymTriDiagEigLT(T,x)>n-k)
        z=x;
    else
        y=x;
    end
end
lambda=(y+z)/2;
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
