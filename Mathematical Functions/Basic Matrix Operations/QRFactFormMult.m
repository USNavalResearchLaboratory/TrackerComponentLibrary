function C=QRFactFormMult(QFactForm,C)
%%QRFACTFORMMULT Given a Q matrix from a QR decomposition in factorial
%           form, multiply the matrix by another matrix C. That is evaluate
%           Q*C where Q is given in factorial form (not explicitly given).
%
%INPUTS: QFactForm The Q matrix from a QR decomposition in factorial form.
%                  This can be obtained from the HouseholderQR function.
%                C The matrix that should right-multiply Q.
%
%OUTPUTS: C The product Q*C.
%
%This algorithm is based on the discussion in Section 5.1.6 of [1]. Using Q
%in factorial form rather than explicitly computing Q is more efficient.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(QFactForm,1);
n=size(QFactForm,2);
v=zeros(m,1);

%Equation 5.1.4
for j=min(n,m-1):-1:1
    v(j:m)=[1;QFactForm((j+1):m,j)];
    beta=2/(1+QFactForm((j+1):m,j)'*QFactForm((j+1):m,j));
    C(j:m,:)=C(j:m,:)-(beta*v(j:m))*(v(j:m)'*C(j:m,:));
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
