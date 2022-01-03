function xy=solveBivarMultilinSys(a,c)
%%SOLVEBIVARMULTILINSYS Given bivariate multilinear equations of the form
%       c(i)=a(1,i)+a(2,i)*x+a(3,i)*y+a(4,i)*x*y
%       for i=1,2, solve for x and y.
%
%INPUTS: a A 4X2 vector of coefficients as defined above.
%        c A length 2 vector of the values of the multilinear equations.
%
%OUTPUTS: xy A 2X2 vector of [x;y] solutions. xy(:,i) is the ith solution.
%
%EXAMPLE:
%We generate a random system and obtain both solutions. One solutions
%agrees with the [x;y] values used to generate the original c values and
%the other solution does not, though it does agree with c values.
% a=100*randn(4,2);
% x=1;
% y=-12;
% c=zeros(2,1);
% for k=1:2
%     c(k)=a(1,k)+a(2,k)*x+a(3,k)*y+a(4,k)*x*y; 
% end
% xy=solveBivarMultilinSys(a,c)%One solution will be [1;-12]
%
%April 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

c1=c(1);
c2=c(2);
a11=a(1,1);
a21=a(2,1);
a31=a(3,1);
a41=a(4,1);
a12=a(1,2);
a22=a(2,2);
a32=a(3,2);
a42=a(4,2);

commonTerm=-(a12*a41-a11*a42+a42*c1-a41*c2);
addTermX=a21*a32-a22*a31+commonTerm;
addTermY=a22*a31-a21*a32+commonTerm;
radTerm=sqrt(-4*(a22*a41-a21*a42)*(a12*a31-a11*a32+a32*c1-a31*c2)+(a22*a31-a21*a32+a12*a41-a11*a42+a42*c1-a41*c2)^2);

xDenom=2*(a22*a41-a21*a42);
yDenom=2*(a32*a41-a31*a42);

x1=(addTermX-radTerm)/xDenom;
y1=(addTermY+radTerm)/yDenom;
x2=(addTermX+radTerm)/xDenom;
y2=(addTermY-radTerm)/yDenom;

xy=[x1,x2;
    y1,y2];

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
