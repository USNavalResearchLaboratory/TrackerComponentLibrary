function [xi,w]=seventhOrderWeightEllipsCubPoints(a,b)
%%SEVENTHORDERWEIGHTELLIPSCUBPOINTS Generating seventh-order cubature
%               points for integrating over a 2D elliptical region (aligned
%               with the coordinate axes, x and y) having weighting
%               function w(x)=1/(sqrt((x-c)^2+y^2)*sqrt((x+c)^2+y^2) where
%               c=sqrt(a^2-b^2) and a and b are the parameters of the
%               ellipse such that x^2/a^2+y^2/b^2=1 with b<=a.
%
%INPUTS: a,b The positive, real, scalar parameters of the ellipse such that
%            x^2/a^2+y^2/b^2=1 defines the ellipse and b<=a.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This uses algorithm ELP 7-1 in [1], pg. 337, 12 points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(b>a)
    error('b must be <=a')
end

c=sqrt(a^2-b^2);
V=2*pi*log((a+b)/c);

I20=(pi/2)*(a*b+(1/(2*pi))*c^2*V);
I02=(pi/2)*(a*b-(1/(2*pi))*c^2*V);
I40=(3*pi/8)*((1/2)*a*b^3+(5/4)*c^2*a*b+(3/(8*pi))*c^4*V);
I22=(pi/8)*((1/2)*a*b^3+(1/4)*c^2*a*b-(1/(8*pi))*c^4*V);
I04=(3*pi/8)*((1/2)*a*b^3-(3/4)*c^2*a*b+(3/(8*pi))*c^4*V);
I60=(15*pi/48)*((1/3)*a*b^5+(13/12)*c^2*a*b^3+(11/8)*c^4*a*b+(5/(16*pi))*c^6*V);
I42=(3*pi/48)*((1/3)*a*b^5+(7/12)*c^2*a*b^3+(1/8)*c^4*a*b-(1/(16*pi))*c^6*V);
I24=(3*pi/48)*((1/3)*a*b^5+(1/12)*c^2*a*b^3-(1/8)*c^4*a*b+(1/(16*pi))*c^6*V);
I06=(15*pi/48)*((1/3)*a*b^5-(5/12)*c^2*a*b^3+(5/8)*c^4*a*b-(5/(16*pi))*c^6*V);

u=sqrt(I42/I22);
v=sqrt(I24/I22);
C=I22^3/(4*I42*I24);
k1=2*(V-4*C)/3;

D1=k1*(I40-4*C*u^4)-(I20-4*C*u^2)^2;

a0=((I20-4*C*u^2)*(I60-4*C*u^6)-(I40-4*C*u^4)^2)/D1;
a1=((I20-4*C*u^2)*(I40-4*C*u^4)-k1*(I60-4*C*u^6))/D1;

%Get the squared values of r1 and r2 in the paper.
rRoots=roots([1;a1;a0]);
r1=rRoots(1);
r2=rRoots(2);

A1=(I20-4*C*u^2-k1*r2)/(2*(r1-r2));
A2=(I20-4*C*u^2-k1*r1)/(2*(r2-r1));

k2=(V-4*C)/3;
D2=k2*(I04-4*C*v^4)-(I02-4*C*v^2)^2;

b0=((I02-4*C*v^2)*(I06-4*C*v^6)-(I04-4*C*v^4)^2)/D2;
b1=((I02-4*C*v^2)*(I04-4*C*v^4)-k2*(I06-4*C*v^6))/D2;

%Get the squared values of s1 and s2 in the paper.
sRoots=roots([1;b1;b0]);
s1=sRoots(1);
s2=sRoots(2);

B1=(I02-4*C*v^2-k2*s2)/(2*(s1-s2));
B2=(I02-4*C*v^2-k2*s1)/(2*(s2-s1));

r1=sqrt(r1);
r2=sqrt(r2);
s1=sqrt(s1);
s2=sqrt(s2);

xi=[[r1;0],[-r1;0],[r2;0],[-r2;0],[0;s1],[0;-s1],[0;s2],[0;-s2],PMCombos([u;v])];
w=[A1;A1;A2;A2;B1;B1;B2;B2;C;C;C;C];
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
