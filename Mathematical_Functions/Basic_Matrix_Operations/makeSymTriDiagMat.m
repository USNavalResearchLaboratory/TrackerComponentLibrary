function T=makeSymTriDiagMat(a,b)
%%MAKESYMTRIDIAGMAT Create a symmetric tri-diagonal matrix with the
%     elements of the vector a on the main diagonal and the elements of the
%     vector b on the diagonals just above and below the main diagonal.
%     That is:
%     T=[a(1),b(1),   0, ...,   0;
%        b(1),a(2),b(2), ...,   0;
%           0,b(2),a(3), ...,   0;
%         ...,...,  ..., ..., ...;
%          0, ...,  ..., ...,a(n)];
%
%INPUTS: a A length n vector.
%        b A length n-1 vector.
%
%OUTPUTS: T The nXn symmetric tridiagonal matrix as described above.
%
%EXAMPLE:
% a=[1;2;3;4];
% b=[12;24;36];
% T=makeSymTriDiagMat(a,b)
%The above is:
%T=[1, 12,  0,     0;
%   12,  2, 24,     0;
%    0, 24,  3,    36;
%    0,  0, 36,     4];
%
%September 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(a);

T=zeros(n,n);
%Set the main diagonal.
T(1:(n+1):(n^2))=a;
%Set the lower diagonal.
T(2:(n+1):(n^2-n))=b;
%Set the upper diagonal.
T((n+1):(n+1):(n^2-1))=b;

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
