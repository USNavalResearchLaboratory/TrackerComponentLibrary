function [xi,w]=fifthOrderWeightEllipsCubPoints(a,b)
%%FIFTHORDERWEIGHTELLIPSCUBPOINTS Generating fifth-order cubature points
%               for integrating over a 2D elliptical region (aligned with
%               the coordinate axes, x and y) having weighting function
%               w(x)=1/(sqrt((x-c)^2+y^2)*sqrt((x+c)^2+y^2) where
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
%This uses algorithm ELP 5-1 in [1], pg. 336, 7 points.
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

r=sqrt((I40*I04-I22^2)/(I20*I04-I02*I22));
s=sqrt(I22/I02);
t=sqrt(I04/I02);

B=(I20*I04-I02*I22)^2/(2*I04*(I40*I04-I22^2));
C=I02^2/(4*I04);
A=V-2*B-4*C;

xi=[[0;0],PMCombos([r;0]),PMCombos([s;t])];
w=[A;B;B;C;C;C;C];

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
