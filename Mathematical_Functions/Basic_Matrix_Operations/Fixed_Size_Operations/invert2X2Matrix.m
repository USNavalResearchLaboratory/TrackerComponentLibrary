function AInv=invert2X2Matrix(A)
%%INVERT2X2MATRIX Inver the 2X2 matrix A. This function demonstrates the
%       algorithm given in [1], which is more efficient than Gaussian
%       elimination. This could be used as a template for implementation in
%       other languages. In Matlab, this function is not faster than the
%       built-in inv function.
%
%INPUTS: A A 2X2 real or complex invertible matrix.
%
%OUTPUTS: AInv The 2X2 matrix inverse of A. 
%
%EXAMPLE:
%This just shows that this function produces the same result as Matlab
%within finite precision limitations.
% A=randn(2,2);
% C=invert2X2Matrix(A);
% C1=inv(A);%Matlab's way.
% RelErr=max(max(abs(abs((C-C1)./C1))))
%
%REFERENCES:
%[1] V. Strassen, "Gaussian elimination is not optimal," Numerische
%    Mathematik, vol. 13, pp. 354-356, 1965.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

I=1/A(1,1);
II=A(2,1)*I;
III=I*A(1,2);
IV=A(2,1)*III;
V=IV-A(2,2);
VI=1/V;

AInv=zeros(2,2);
AInv(1,2)=III*VI;
AInv(2,1)=VI*II;
VII=III*AInv(2,1);
AInv(1,1)=I-VII;
AInv(2,2)=-VI;

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
