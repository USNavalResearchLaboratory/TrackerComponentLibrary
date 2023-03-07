function C=multiply3X3And3X3Matrices(A,B)
%%MULTIPLY3X3And3X3Matrices This just multiplies the two 3X3 matrices A and
%           B to get a 3X3 matrix C. This demonstrates the algorithm
%           developed in [1], which minimizes the number of scalar
%           multiplication operations. This could be used as a template for
%           efficient implementation in other programming languages. In
%           Matlab, this is not faster than the built-in matrix
%           multiplication operation.
%
%INPUTS: A A 3X3 matrix.
%        B A 3X3 matrix.
%
%OUTPUTS: C A 3X3 matrix.
%
%EXAMPLE:
%This just shows that this function produces the same result as Matlab
%within finite precision limitations.
% A=randn(3,3);
% B=randn(3,3);
% C=multiply3X3And3X3Matrices(A,B);
% C1=A*B;%Matlab's way.
% RelErr=max(max(abs((C-C1)./C1)))
%
%REFERENCES:
%[2] J. D. Laderman, "A noncommutative algorithm for multiplying 3 X 3
%matrices using 23 multiplications," Bulletin of the American Mathematical
%Society, vol. 82, no. 1, pp. 126-128, Jan. 1976.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

a11=A(1,1);
a12=A(1,2);
a13=A(1,3);
a21=A(2,1);
a22=A(2,2);
a23=A(2,3);
a31=A(3,1);
a32=A(3,2);
a33=A(3,3);

b11=B(1,1);
b12=B(1,2);
b13=B(1,3);
b21=B(2,1);
b22=B(2,2);
b23=B(2,3);
b31=B(3,1);
b32=B(3,2);
b33=B(3,3);

m1=(a11+a12+a13-a21-a22-a32-a33)*b22;
m2=(a11-a21)*(-b12+b22);
m3=a22*(-b11+b12+b21-b22-b23-b31+b33);
m4=(-a11+a21+a22)*(b11-b12+b22);
m5=(a21+a22)*(-b11+b12);
m6=a11*b11;
m7=(-a11+a31+a32)*(b11-b13+b23);
m8=(-a11+a31)*(b13-b23);
m9=(a31+a32)*(-b11+b13);
m10=(a11+a12+a13-a22-a23-a31-a32)*b23;
m11=a32*(-b11+b13+b21-b22-b23-b31+b32);
m12=(-a13+a32+a33)*(b22+b31-b32);
m13=(a13-a33)*(b22-b32);
m14=a13*b31;
m15=(a32+a33)*(-b31+b32);
m16=(-a13+a22+a23)*(b23+b31-b33);
m17=(a13-a23)*(b23-b33);
m18=(a22+a23)*(-b31+b33);
m19=a12*b21;
m20=a23*b32;
m21=a21*b13;
m22=a31*b12;
m23=a33*b33;

C=zeros(3,3);
C(1,1)=m6+m14+m19;
C(1,2)=m1+m4+m5+m6+m12+m14+m15;
C(1,3)=m6+m7+m9+m10+m14+m16+m18;
C(2,1)=m2+m3+m4+m6+m14+m16+m17;
C(2,2)=m2+m4+m5+m6+m20;
C(2,3)=m14+m16+m17+m18+m21;
C(3,1)=m6+m7+m8+m11+m12+m13+m14;
C(3,2)=m12+m13+m14+m15+m22;
C(3,3)=m6+m7+m8+m9+m23;

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
