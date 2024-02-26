function count=numSymTriDiagEigLT(T,lambda)
%%NUMSYMTRIDIAGEIGLT Given a real symmetric tri-diagonal matrix T with no 0
%       subdiagonal values, find the number of eigenvalues in T less than
%       the real value lambda without performing an eigenvalue
%       decomposition on T.
%
%INPUTS: T An nXn real symmetric tridiagonal matrix where no subdiagonal
%          terms are 0 (which would lead to repeated roots).
%   lambda The real value below which we wish to count the eigenvalues.
%
%OUTPUTS: count The number of eigenvalues of T that are less than lambda.
%
%The algorithm for counting the number of eigenvalues above or below a
%particular value when dealing with a real symmetric tridiagonal matrix
%comes from the number of sign changes in some polynomial defined in
%Chapter 8.4.2 of [1], where it is noted that if there are no zero
%subdiagonal elements, then the eigenvalues are unique.
%
%EXAMPLE:
%This generates a random symmetric tridiagonal matrix. The number of
%eigenvalues below 1 are counted using this function and by performing an
%eigenvalue decomposition and a comparison. The results are equal.
% n=10;
% a=randn(n,1);
% b=randn(n-1,1);
% T=makeSymTriDiagMat(a,b);
% lambda=1;
% numSymTriDiagEigLT(T,lambda)
% sum(eig(T)<lambda)
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(T,1);
prm2=1;
prm1=T(1,1)-lambda;
count=0;
sign=1;
if(prm1<=0)
    count=1;
    sign=-sign;
end  
for r=2:n
    %Equation 8.4.2.
    pr=(T(r,r)-lambda)*prm1-prm2*T(r,r-1)^2;
    if(sign*pr<=0)
        count=count+1;
        sign=-sign;
    end
    prm2=prm1;
    prm1=pr;
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
