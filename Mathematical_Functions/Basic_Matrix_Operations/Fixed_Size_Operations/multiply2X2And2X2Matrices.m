function C=multiply2X2And2X2Matrices(A,B)
%%MULTIPLY2X2AND2X2MATRICES This just multiplies the two 2X2 matrices A and
%           B to get a 2X2 matrix C. This demonstrates the algorithm
%           developed in [1], which minimizes the number of scalar
%           multiplication operations. This could be used as a template for
%           efficient implementation in other programming languages. In
%           Matlab, this is not faster than the built-in matrix
%           multiplication operation.
%
%INPUTS: A A 2X2 matrix.
%        B A 2X2 matrix.
%
%OUTPUTS: C A 2X2 matrix.
%
%Though only implemented here for 2X2 matrices, the algorithm in [1] is
%more general than just 2X2 matrices and via repeated application to block
%sub-matrices can provide expressions for the multiplication of matrices
%of dimensions other than 2X2.
%
%EXAMPLE:
%This just shows that this function produces the same result as Matlab
%within finite precision limitations.
% A=randn(2,2);
% B=randn(2,2);
% C=multiply2X2And2X2Matrices(A,B);
% C1=A*B;%Matlab's way.
% RelErr=max(max(abs(abs((C-C1)./C1))))
%
%REFERENCES:
%[1] V. Strassen, "Gaussian elimination is not optimal," Numerische
%    Mathematik, vol. 13, pp. 354-356, 1965.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

I=(A(1,1)+A(2,2))*(B(1,1)+B(2,2));
II=(A(2,1)+A(2,2))*B(1,1);
III=A(1,1)*(B(1,2)-B(2,2));
IV=A(2,2)*(-B(1,1)+B(2,1));
V=(A(1,1)+A(1,2))*B(2,2);
VI=(-A(1,1)+A(2,1))*(B(1,1)+B(1,2));
VII=(A(1,2)-A(2,2))*(B(2,1)+B(2,2));

C=zeros(2,2);
C(1,1)=I+IV-V+VII;
C(2,1)=II+IV;
C(1,2)=III+V;
C(2,2)=I+III-II+VI;

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
