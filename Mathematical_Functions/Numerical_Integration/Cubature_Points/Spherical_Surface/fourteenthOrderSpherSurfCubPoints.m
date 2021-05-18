function [xi,w]=fourteenthOrderSpherSurfCubPoints(numDim)
%%FOURTEENTHORDERSPHERSURFCUBPOINTS Generate fourteenth-order cubature
%               points for integration over the surface of a unit (hyper-)
%               sphere (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. Currently, only numDim=3 is supported.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%Algorithm U3 14-1 in [1], pg. 302, using 72 points, is used.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(numDim~=3)
   error('Only 3D points are supported'); 
end

V=2*pi^(numDim/2)/gamma(numDim/2);
r=sqrt((5-sqrt(5))/10);
s=sqrt((5+sqrt(5))/10);

B=V*125/10080;
C=V*143/10080;

y=roots([2556125;-5112250;3578575;-1043900;115115;-3562;9]);
z=sqrt(y);

u=zeros(5,1);
v=zeros(5,1);
w=zeros(5,1);

u(1)=(-z(3)+z(4))/(2*s);
u(2)=(-z(5)+z(2))/(2*s);
u(3)=(-z(2)+z(6))/(2*s);
u(4)=(-z(6)+z(3))/(2*s);
u(5)=(-z(4)+z(5))/(2*s);

v(1)=(z(5)+z(6))/(2*s);
v(2)=(z(6)+z(4))/(2*s);
v(3)=(z(3)+z(5))/(2*s);
v(4)=(z(4)+z(2))/(2*s);
v(5)=(z(2)+z(3))/(2*s);

w(1)=(z(1)+z(2))/(2*s);
w(2)=(z(1)+z(3))/(2*s);
w(3)=(z(1)+z(4))/(2*s);
w(4)=(z(1)+z(5))/(2*s);
w(5)=(z(1)+z(6))/(2*s);

xi=zeros(3,72);
xi(:,1:4)=PMCombos([r;s;0]);
xi(:,5:8)=PMCombos([0;r;s]);
xi(:,9:12)=PMCombos([s;0;r]);
cur2Add=13;
for i=1:5
    xi(:,cur2Add)=[u(i);v(i);w(i)];
    xi(:,cur2Add+1)=[u(i);-v(i);-w(i)];
    xi(:,cur2Add+2)=[-u(i);-v(i);w(i)];
    xi(:,cur2Add+3)=[-u(i);v(i);-w(i)];
    xi(:,cur2Add+4)=[v(i);w(i);u(i)];
    xi(:,cur2Add+5)=[v(i);-w(i);-u(i)];
    xi(:,cur2Add+6)=[-v(i);-w(i);u(i)];
    xi(:,cur2Add+7)=[-v(i);w(i);-u(i)];
    xi(:,cur2Add+8)=[w(i);u(i);v(i)];
    xi(:,cur2Add+9)=[w(i);-u(i);-v(i)];
    xi(:,cur2Add+10)=[-w(i);-u(i);v(i)];
    xi(:,cur2Add+11)=[-w(i);u(i);-v(i)];
    
    cur2Add=cur2Add+12;
end

w=[B*ones(12,1);C*ones(60,1)];

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
