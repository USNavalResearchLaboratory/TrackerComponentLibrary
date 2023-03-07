function C=multiply4X5And5X5Matrices(A,B)
%%MULTIPLY4X4And5X5Matrices This just multiplies the 4X4 matrix A by the
%           5X5 matrix B to get a 4X5 matrix C. This demonstrates the
%           algorithm developed in [1], which minimizes the number of
%           scalar multiplication operations. This could be used as a
%           template for efficient implementation in other programming
%           languages. In Matlab, this is not faster than the built-in
%           matrix multiplication operation.
%
%INPUTS: A A 4X5 matrix.
%        B A 5X5 matrix.
%
%OUTPUTS: C A 4X5 matrix.
%
%EXAMPLE:
%This just shows that this function produces the same result as Matlab
%within finite precision limitations.
% A=randn(4,5);
% B=randn(5,5);
% C=multiply4X5And5X5Matrices(A,B);
% C1=A*B;%Matlab's way.
% RelErr=max(max(abs(abs((C-C1)./C1))))
%
%REFERENCES:
%[1] A. Fawzi, M. Balog, A. Huang, T. Hubert, B. Romera-Paredes,
%    M. Barekatain, A. Novikov, F. J. R. Ruiz, J. Schrittwieser, G.
%    Swirszcz, D. Silver, D. Hassabis, and P. Kohli, "Discovering faster
%    matrix multiplication algorithms with reinforcement learning," Nature,
%    vol. 610, pp. 47-53, 2022.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

h1=A(3,2)*(-B(2,1)-B(2,5)-B(3,1));
h2=(A(2,2)+A(2,5)-A(3,5))*(-B(2,5)-B(5,1));
h3=(-A(3,1)-A(4,1)+A(4,2))*(-B(1,1)+B(2,5));
h4=(A(1,2)+A(1,4)+A(3,4))*(-B(2,5)-B(4,1));
h5=(A(1,5)+A(2,2)+A(2,5))*(-B(2,4)+B(5,1));
h6=(-A(2,2)-A(2,5)-A(4,5))*(B(2,3)+B(5,1));
h7=(-A(1,1)+A(4,1)-A(4,2))*(B(1,1)+B(2,4));
h8=(A(3,2)-A(3,3)-A(4,3))*(-B(2,3)+B(3,1));
h9=(-A(1,2)-A(1,4)+A(4,4))*(B(2,3)+B(4,1));
h10=(A(2,2)+A(2,5))*B(5,1);
h11=(-A(2,1)-A(4,1)+A(4,2))*(-B(1,1)+B(2,2));
h12=(A(4,1)-A(4,2))*B(1,1);
h13=(A(1,2)+A(1,4)+A(2,4))*(B(2,2)+B(4,1));
h14=(A(1,3)-A(3,2)+A(3,3))*(B(2,4)+B(3,1));
h15=(-A(1,2)-A(1,4))*B(4,1);
h16=(-A(3,2)+A(3,3))*B(3,1);
h17=(A(1,2)+A(1,4)-A(2,1)+A(2,2)-A(2,3)+A(2,4)-A(3,2)+A(3,3)-A(4,1)+A(4,2))*B(2,2);
h18=A(2,1)*(B(1,1)+B(1,2)+B(5,2));
h19=-A(2,3)*(B(3,1)+B(3,2)+B(5,2));
h20=(-A(1,5)+A(2,1)+A(2,3)-A(2,5))*(-B(1,1)-B(1,2)+B(1,4)-B(5,2));
h21=(A(2,1)+A(2,3)-A(2,5))*B(5,2);
h22=(A(1,3)-A(1,4)-A(2,4))*(B(1,1)+B(1,2)-B(1,4)-B(3,1)-B(3,2)+B(3,4)+B(4,4));
h23=A(1,3)*(-B(3,1)+B(3,4)+B(4,4));
h24=A(1,5)*(-B(4,4)-B(5,1)+B(5,4));
h25=-A(1,1)*(B(1,1)-B(1,4));
h26=(-A(1,3)+A(1,4)+A(1,5))*B(4,4);
h27=(A(1,3)-A(3,1)+A(3,3))*(B(1,1)-B(1,4)+B(1,5)+B(3,5));
h28=-A(3,4)*(-B(3,5)-B(4,1)-B(4,5));
h29=A(3,1)*(B(1,1)+B(1,5)+B(3,5));
h30=(A(3,1)-A(3,3)+A(3,4))*B(3,5);
h31=(-A(1,4)-A(1,5)-A(3,4))*(-B(4,4)-B(5,1)+B(5,4)-B(5,5));
h32=(A(2,1)+A(4,1)+A(4,4))*(B(1,3)-B(4,1)-B(4,2)-B(4,3));
h33=A(4,3)*(-B(3,1)-B(3,3));
h34=A(4,4)*(-B(1,3)+B(4,1)+B(4,3));
h35=-A(4,5)*(B(1,3)+B(5,1)+B(5,3));
h36=(A(2,3)-A(2,5)-A(4,5))*(B(3,1)+B(3,2)+B(3,3)+B(5,2));
h37=(-A(4,1)-A(4,4)+A(4,5))*B(1,3);
h38=(-A(2,3)-A(3,1)+A(3,3)-A(3,4))*(B(3,5)+B(4,1)+B(4,2)+B(4,5));
h39=(-A(3,1)-A(4,1)-A(4,4)+A(4,5))*(B(1,3)+B(5,1)+B(5,3)+B(5,5));
h40=(-A(1,3)+A(1,4)+A(1,5)-A(4,4))*(-B(3,1)-B(3,3)+B(3,4)+B(4,4));
h41=(-A(1,1)+A(4,1)-A(4,5))*(B(1,3)+B(3,1)+B(3,3)-B(3,4)+B(5,1)+B(5,3)-B(5,4));
h42=(-A(2,1)+A(2,5)-A(3,5))*(-B(1,1)-B(1,2)-B(1,5)+B(4,1)+B(4,2)+B(4,5)-B(5,2));
h43=A(2,4)*(B(4,1)+B(4,2));
h44=(A(2,3)+A(3,2)-A(3,3))*(B(2,2)-B(3,1));
h45=(-A(3,3)+A(3,4)-A(4,3))*(B(3,5)+B(4,1)+B(4,3)+B(4,5)+B(5,1)+B(5,3)+B(5,5));
h46=-A(3,5)*(-B(5,1)-B(5,5));
h47=(A(2,1)-A(2,5)-A(3,1)+A(3,5))*(B(1,1)+B(1,2)+B(1,5)-B(4,1)-B(4,2)-B(4,5));
h48=(-A(2,3)+A(3,3))*(B(2,2)+B(3,2)+B(3,5)+B(4,1)+B(4,2)+B(4,5));
h49=(-A(1,1)-A(1,3)+A(1,4)+A(1,5)-A(2,1)-A(2,3)+A(2,4)+A(2,5))*(-B(1,1)-B(1,2)+B(1,4));
h50=(-A(1,4)-A(2,4))*(B(2,2)-B(3,1)-B(3,2)+B(3,4)-B(4,2)+B(4,4));
h51=A(2,2)*(B(2,1)+B(2,2)-B(5,1));
h52=A(4,2)*(B(1,1)+B(2,1)+B(2,3));
h53=-A(1,2)*(-B(2,1)+B(2,4)+B(4,1));
h54=(A(1,2)+A(1,4)-A(2,2)-A(2,5)-A(3,2)+A(3,3)-A(4,2)+A(4,3)-A(4,4)-A(4,5))*B(2,3);
h55=(A(1,4)-A(4,4))*(-B(2,3)+B(3,1)+B(3,3)-B(3,4)+B(4,3)-B(4,4));
h56=(A(1,1)-A(1,5)-A(4,1)+A(4,5))*(B(3,1)+B(3,3)-B(3,4)+B(5,1)+B(5,3)-B(5,4));
h57=(-A(3,1)-A(4,1))*(-B(1,3)-B(1,5)-B(2,5)-B(5,1)-B(5,3)-B(5,5));
h58=(-A(1,4)-A(1,5)-A(3,4)-A(3,5))*(-B(5,1)+B(5,4)-B(5,5));
h59=(-A(3,3)+A(3,4)-A(4,3)+A(4,4))*(B(4,1)+B(4,3)+B(4,5)+B(5,1)+B(5,3)+B(5,5));
h60=(A(2,5)+A(4,5))*(B(2,3)-B(3,1)-B(3,2)-B(3,3)-B(5,2)-B(5,3));
h61=(A(1,4)+A(3,4))*(B(1,1)-B(1,4)+B(1,5)-B(2,5)-B(4,4)+B(4,5)-B(5,1)+B(5,4)-B(5,5));
h62=(A(2,1)+A(4,1))*(B(1,2)+B(1,3)+B(2,2)-B(4,1)-B(4,2)-B(4,3));
h63=(-A(3,3)-A(4,3))*(-B(2,3)-B(3,3)-B(3,5)-B(4,1)-B(4,3)-B(4,5));
h64=(A(1,1)-A(1,3)-A(1,4)+A(3,1)-A(3,3)-A(3,4))*(B(1,1)-B(1,4)+B(1,5));
h65=(-A(1,1)+A(4,1))*(-B(1,3)+B(1,4)+B(2,4)-B(5,1)-B(5,3)+B(5,4));
h66=(A(1,1)-A(1,2)+A(1,3)-A(1,5)-A(2,2)-A(2,5)-A(3,2)+A(3,3)-A(4,1)+A(4,2))*B(2,4);
h67=(A(2,5)-A(3,5))*(B(1,1)+B(1,2)+B(1,5)-B(2,5)-B(4,1)-B(4,2)-B(4,5)+B(5,2)+B(5,5));
h68=(A(1,1)+A(1,3)-A(1,4)-A(1,5)-A(4,1)-A(4,3)+A(4,4)+A(4,5))*(-B(3,1)-B(3,3)+B(3,4));
h69=(-A(1,3)+A(1,4)-A(2,3)+A(2,4))*(-B(2,4)-B(3,1)-B(3,2)+B(3,4)-B(5,2)+B(5,4));
h70=(A(2,3)-A(2,5)+A(4,3)-A(4,5))*(-B(3,1)-B(3,2)-B(3,3));
h71=(-A(3,1)+A(3,3)-A(3,4)+A(3,5)-A(4,1)+A(4,3)-A(4,4)+A(4,5))*(-B(5,1)-B(5,3)-B(5,5));
h72=(-A(2,1)-A(2,4)-A(4,1)-A(4,4))*(B(4,1)+B(4,2)+B(4,3));
h73=(A(1,3)-A(1,4)-A(1,5)+A(2,3)-A(2,4)-A(2,5))*(B(1,1)+B(1,2)-B(1,4)+B(2,4)+B(5,2)-B(5,4));
h74=(A(2,1)-A(2,3)+A(2,4)-A(3,1)+A(3,3)-A(3,4))*(B(4,1)+B(4,2)+B(4,5));
h75=-(A(1,2)+A(1,4)-A(2,2)-A(2,5)-A(3,1)+A(3,2)+A(3,4)+A(3,5)-A(4,1)+A(4,2))*B(2,5);
h76=(A(1,3)+A(3,3))*(-B(1,1)+B(1,4)-B(1,5)+B(2,4)+B(3,4)-B(3,5));

C=zeros(4,5);
C(1,1)=-h10+h12+h14-h15-h16+h53+h5-h66-h7;
C(2,1)=h10+h11-h12+h13+h15+h16-h17-h44+h51;
C(3,1)=h10-h12+h15+h16-h1+h2+h3-h4+h75;
C(4,1)=-h10+h12-h15-h16+h52+h54-h6-h8+h9;
C(1,2)=h13+h15+h20+h21-h22+h23+h25-h43+h49+h50;
C(2,2)=-h11+h12-h13-h15-h16+h17+h18-h19-h21+h43+h44;
C(3,2)=-h16-h19-h21-h28-h29-h38+h42+h44-h47+h48;
C(4,2)=h11-h12-h18+h21-h32+h33-h34-h36+h62-h70;
C(1,3)=h15+h23+h24+h34-h37+h40-h41+h55-h56-h9;
C(2,3)=-h10+h19+h32+h35+h36+h37-h43-h60-h6-h72;
C(3,3)=-h16-h28+h33+h37-h39+h45-h46+h63-h71-h8;
C(4,3)=h10+h15+h16-h33+h34-h35-h37-h54+h6+h8-h9;
C(1,4)=-h10+h12+h14-h16+h23+h24+h25+h26+h5-h66-h7;
C(2,4)=h10+h18-h19+h20-h22-h24-h26-h5-h69+h73;
C(3,4)=-h14+h16-h23-h26+h27+h29+h31+h46-h58+h76;
C(4,4)=h12+h25+h26-h33-h35-h40+h41+h65-h68-h7;
C(1,5)=h15+h24+h25+h27-h28+h30+h31-h4+h61+h64;
C(2,5)=-h10-h18-h2-h30-h38+h42-h43+h46+h67+h74;
C(3,5)=-h10+h12-h15+h28+h29-h2-h30-h3+h46+h4-h75;
C(4,5)=-h12-h29+h30-h34+h35+h39+h3-h45+h57+h59;

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
