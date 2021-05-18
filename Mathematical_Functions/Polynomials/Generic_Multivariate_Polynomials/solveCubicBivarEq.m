function sols=solveCubicBivarEq(coeffsA,coeffsB,AbsTol)
%%SOLVECUBICBIVAREQ Find the zeros of a set of two cubic bivariate
%                  equations. This assumes that the systems are well
%                  defined so there exists a finite number of solutions.
%
%INPUTS: coeffsA, coeffsB Hypermatrices of the coefficients for the
%               multivariate polynomials whose zeros are desired. These are
%               arranged such that coeffsA(a1,a2) corresponds to the
%               coefficient of an x1^(a1-1)*x2^(a2-1) term in the first
%               polynomial. Entries corresponding to degrees higher than
%               two are ignored. For example coeffsA(4+1,4+1) is ignored.
%        AbsTol A tolerance value for determining whether the sum of the
%               magnitudes of the two equations is zero. If omitted or an
%               empty matrix is passed, the default of 1e-7 is used. There
%               is a loss of precision when repeated roots arise.
%
%OUTPUTS: sols A 2XnumSol set of the solutions. sols(1,:) are the x
%              solutions; sols(2,:) are the corresponding y solutions.
%              There can be up to nine solutions.
%
%The problem is solved by eliminating one variable resulting in a ninth-
%degree univariate polynomial equation that can be solved with the roots
%function. Given the solutions to the 9th degree polynomial, the results
%can be substituted back into the cubic equations. However, the cubics
%produce extra solutions. Normally, one would just choose the solutions
%that agrees with both cubic equations. However, that is problematic if
%there are has repeated roots. Thus, we save all solutions below AbsTol in
%agreement, but we can only return nine solutions, because Bézout's number
%is nine. Thus, if more than nine solutions are saved, we take the lowest-
%cost solution and then sequentially take the others such that they
%maximize the minimum distance to previous solutions. In the end, however,
%finite precision issues can occasionally cause problems when repeated
%roots are present.
%
%EXAMPLE:
%This is an example of obtaining solutions to a cubic bivariate system of
%equations.
% coeffsA=zeros(3+1,3+1);
% coeffsA(3+1,0+1)=1;
% coeffsA(1+1,2+1)=-3;
% coeffsA(1+1,0+1)=6;
% coeffsA(0+1,1+1)=-3;
% coeffsA(0+1,0+1)=5;
% 
% coeffsB=zeros(3+1,3+1);
% coeffsB(0+1,3+1)=1;
% coeffsB(2+1,1+1)=-3;
% coeffsB(0+1,1+1)=-6;
% coeffsB(1+1,0+1)=-3;
% coeffsB(0+1,0+1)=-7;
% sols=solveCubicBivarEq(coeffsA,coeffsB);
% 
% costVals=zeros(9,1);
% for curVal=1:9
%     costVals(curVal)=abs(polyValMultiDim(coeffsA,sols(:,curVal)))+abs(polyValMultiDim(coeffsB,sols(:,curVal)));
% end
% costVals
%One will see that the residual costs of all 9 solutions are all less then
%1e-13.
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(AbsTol))
    AbsTol=1e-7;%Absolute tolerance for declaring equations zero.
end

%Normalize coefficients to reduce the influence of scaling issues.
coeffsA=coeffsA/max(abs(coeffsA(:)));
coeffsB=coeffsB/max(abs(coeffsB(:)));

a1=coeffsA(3+1,0+1);
a2=coeffsA(0+1,3+1);
a3=coeffsA(2+1,1+1);
a4=coeffsA(1+1,2+1);
a5=coeffsA(2+1,0+1);
a6=coeffsA(0+1,2+1);
a7=coeffsA(1+1,1+1);
a8=coeffsA(1+1,0+1);
a9=coeffsA(0+1,1+1);
a10=coeffsA(0+1,0+1);

b1=coeffsB(3+1,0+1);
b2=coeffsB(0+1,3+1);
b3=coeffsB(2+1,1+1);
b4=coeffsB(1+1,2+1);
b5=coeffsB(2+1,0+1);
b6=coeffsB(0+1,2+1);
b7=coeffsB(1+1,1+1);
b8=coeffsB(1+1,0+1);
b9=coeffsB(0+1,1+1);
b10=coeffsB(0+1,0+1);

c0=-a2^3*b10^3-a2*(3*b10*b2*(a6*a9*b10+a10^2*b2)+b10*(a6^2*b10+a10*a9*b2)*b6+(-2*a10*a6+a9^2)*b10*b6^2+a10^2*b6^3)+a2*(3*a10^2*b2*b6+a9*b10*(2*a9*b2+a6*b6)+a10*(a6*b10*b2+a9*b6^2))*b9-a10*a2*(2*a9*b2+a6*b6)*b9^2+b2*(a6^3*b10^2+b2*(-a9^3*b10+a10^3*b2-a10^2*a9*b6+a10*a9^2*b9)+a6*(3*a10*a9*b10*b2+a9^2*b10*b6+a10^2*b6^2-a10*(2*a10*b2+a9*b6)*b9)+a6^2*(-2*a10*b10*b6-a9*b10*b9+a10*b9^2))+a2^2*(b10*(2*a9*b10*b6+a6*b10*b9-a9*b9^2)+a10*(3*b10^2*b2-3*b10*b6*b9+b9^3));
c1=(-3*a2^3*b10^2*b8+a2^2*(b10*(3*a8*b10*b2+2*a9*b10*b4+2*a7*b10*b6+a6*b10*b7-3*a10*b6*b7+6*a10*b2*b8+4*a9*b6*b8)-(b10*(3*a8*b6+2*a9*b7-2*a6*b8)+3*a10*(b10*b4+b6*b8))*b9-(a7*b10-3*a10*b7+a9*b8)*b9^2+a8*b9^3)+b2*(a10^2*(3*a8*b2^2-a9*b2*b4-a7*b2*b6+2*a6*b4*b6-2*a6*b2*b7)+2*a6^3*b10*b8+a10*a9*b2*(-2*a8*b6+a9*b7+2*a7*b9)+a9^2*b2*(-3*a7*b10-a9*b8+a8*b9)+a6*a9*(3*a8*b10*b2+a9*b10*b4+2*a7*b10*b6+a9*b6*b8-a8*b6*b9)-a10*a6*(b6*(-2*a8*b6+a9*b7)-3*b2*(a7*b10+a9*b8)+(4*a8*b2+a9*b4+a7*b6)*b9)-2*a10*a6^2*(b10*b4+b6*b8-b7*b9)-a6^2*(2*a8*b10*b6+a9*b10*b7+a7*b10*b9+a9*b8*b9-a8*b9^2))+a4*(3*a6^2*b10^2*b2-3*a2*a9*b10^2*b2+3*a10*a9*b10*b2^2+a9^2*b10*b2*b6+2*a10*a2*b10*b6^2+a10^2*b2*b6^2+(a2*b10-a10*b2)*(a2*b10+2*a10*b2+a9*b6)*b9-a10*a2*b6*b9^2-2*a6*(b10*(a2*b10+2*a10*b2)*b6+a9*b10*b2*b9-a10*b2*b9^2))+a2*(-a6^2*b10*(b10*b4+2*b6*b8)+3*a10^2*(-b4*b6^2+b2*b6*b7-b2^2*b8+b2*b4*b9)-a10*(6*a8*b10*b2^2+a9*b10*b2*b4+a7*b10*b2*b6+2*a8*b6^3-a9*b6^2*b7+a9*b2*b6*b8-b6*(6*a8*b2+2*a9*b4+a7*b6)*b9+4*a9*b2*b7*b9+2*a7*b2*b9^2)+a6*(4*a10*b10*b4*b6+2*a8*b10*b6^2+a10*b10*b2*b7+a9*b10*b6*b7-6*a9*b10*b2*b8+2*a10*b6^2*b8+(a8*b10*b2+a9*b10*b4-2*a10*b6*b7+a10*b2*b8+a9*b6*b8)*b9-(a10*b4+a8*b6)*b9^2+a7*b10*(-3*b10*b2+b6*b9))-a9*(a8*b10*b2*b6+2*a7*b10*(b6^2-2*b2*b9)+a8*b9*(-b6^2+2*b2*b9)+a9*(2*b10*b4*b6-2*b10*b2*b7+b6^2*b8-2*b2*b8*b9))));
c2=(-3*a2^3*b10*(b10*b5+b8^2)+a2*(-3*a3*a6*b10^2*b2-3*a4*a7*b10^2*b2-6*a10*a5*b10*b2^2-3*a8^2*b10*b2^2+a10*a6*b10*b2*b3+2*a9^2*b10*b2*b3-a10*a7*b10*b2*b4-a8*a9*b10*b2*b4+2*a10*a6*b10*b4^2-a9^2*b10*b4^2-6*a6*a9*b10*b2*b5-3*a10^2*b2^2*b5-a4^2*b10^2*b6-a10*a3*b10*b2*b6-a7*a8*b10*b2*b6-a5*a9*b10*b2*b6+a6*a9*b10*b3*b6+3*a10^2*b2*b3*b6+4*a10*a4*b10*b4*b6+4*a6*a8*b10*b4*b6-4*a7*a9*b10*b4*b6-3*a10^2*b4^2*b6-2*a6^2*b10*b5*b6-a10*a9*b2*b5*b6+2*a5*a6*b10*b6^2-a7^2*b10*b6^2+2*a4*a8*b10*b6^2-2*a3*a9*b10*b6^2+a10*a9*b3*b6^2-6*a10*a8*b4*b6^2+2*a10*a6*b5*b6^2-a9^2*b5*b6^2-2*a10*a5*b6^3-a8^2*b6^3+a10*a4*b10*b2*b7+a6*a8*b10*b2*b7+4*a7*a9*b10*b2*b7+a6*a9*b10*b4*b7+3*a10^2*b2*b4*b7+a6*a7*b10*b6*b7+a4*a9*b10*b6*b7+6*a10*a8*b2*b6*b7+2*a10*a9*b4*b6*b7+a10*a7*b6^2*b7+a8*a9*b6^2*b7-2*a10*a9*b2*b7^2-a10*a6*b6*b7^2-6*a6*a7*b10*b2*b8-6*a4*a9*b10*b2*b8-6*a10*a8*b2^2*b8-2*a6^2*b10*b4*b8-a10*a9*b2*b4*b8-a10*a7*b2*b6*b8-a8*a9*b2*b6*b8+4*a10*a6*b4*b6*b8-2*a9^2*b4*b6*b8+2*a10*a4*b6^2*b8+2*a6*a8*b6^2*b8-2*a7*a9*b6^2*b8+a10*a6*b2*b7*b8+2*a9^2*b2*b7*b8+a6*a9*b6*b7*b8-3*a6*a9*b2*b8^2-a6^2*b6*b8^2-2*a4*a6*b10*(b10*b4+2*b6*b8)+(a5*a6*b10*b2+2*a7^2*b10*b2+a4*a8*b10*b2+4*a3*a9*b10*b2-4*a10*a9*b2*b3+a4*a9*b10*b4+6*a10*a8*b2*b4+a10*a9*b4^2+a10*a6*b2*b5+2*a9^2*b2*b5+a3*a6*b10*b6+3*a8^2*b2*b6-2*a10*a6*b3*b6+2*a8*a9*b4*b6+a6*a9*b5*b6+a10*a3*b6^2+a5*b6*(6*a10*b2+a9*b6)-4*a8*a9*b2*b7-2*a10*a6*b4*b7-2*a10*a4*b6*b7-2*a6*a8*b6*b7+(a10*a4*b2+a6*a8*b2+a6*a9*b4+a4*a9*b6)*b8+a7*(a6*b10*b4+a4*b10*b6+2*a10*b4*b6+a8*b6^2-4*a10*b2*b7+4*a9*b2*b8+a6*b6*b8))*b9-(2*a10*a3*b2+2*a7*a8*b2+2*a5*a9*b2+a10*a4*b4+a6*a8*b4+a5*a6*b6+a4*a8*b6)*b9^2)+a2^2*(2*a7*b10^2*b4+6*a10*b10*b2*b5+2*a3*b10^2*b6-3*a10*b10*b3*b6+4*a9*b10*b5*b6+a4*b10^2*b7-3*a10*b10*b4*b7-3*a8*b10*b6*b7-a9*b10*b7^2+6*a8*b10*b2*b8+4*a9*b10*b4*b8+4*a7*b10*b6*b8-3*a10*b6*b7*b8+3*a10*b2*b8^2+2*a9*b6*b8^2-(2*a9*b10*b3+3*a8*b10*b4+3*a10*b5*b6+2*a7*b10*b7-3*a10*b7^2+(-2*a4*b10+3*a10*b4+3*a8*b6+2*a9*b7)*b8)*b9-(a3*b10-3*a10*b3+a9*b5-3*a8*b7+a7*b8)*b9^2+a6*(b10^2*b3+2*b10*b7*b8+2*b10*b5*b9+b8^2*b9)+a5*(3*b10^2*b2-3*b10*b6*b9+b9^3))+b2*(3*a6*a7*a8*b10*b2+3*a5*a6*a9*b10*b2-3*a7^2*a9*b10*b2-3*a3*a9^2*b10*b2-a6^2*a9*b10*b3-2*a6^2*a8*b10*b4+2*a6*a7*a9*b10*b4+2*a6^3*b10*b5-a9^3*b2*b5-2*a5*a6^2*b10*b6+a6*a7^2*b10*b6+2*a3*a6*a9*b10*b6-a8^2*a9*b2*b6+a6*a9^2*b5*b6+a6*a8^2*b6^2+a10^2*(3*a5*b2^2+a6*(-2*b2*b3+b4^2)-b2*(a7*b4+a3*b6))-a6^2*a7*b10*b7+a8*a9^2*b2*b7-a6*a8*a9*b6*b7+3*a6*a8*a9*b2*b8-3*a7*a9^2*b2*b8+a6*a9^2*b4*b8-2*a6^2*a8*b6*b8+2*a6*a7*a9*b6*b8-a6^2*a9*b7*b8+a6^3*b8^2-(a3*a6^2*b10-a9*(2*a7*a8+a5*a9)*b2+a6*(2*a8^2*b2+a8*a9*b4+a7*a8*b6+a5*a9*b6)+a6^2*(a9*b5-2*a8*b7+a7*b8))*b9+a5*a6^2*b9^2+a10*(3*a8^2*b2^2+a9^2*b2*b3+3*a6*a9*b2*b5-2*a5*a9*b2*b6-a6*a9*b3*b6-2*a6^2*b5*b6+2*a5*a6*b6^2+2*a7*a9*b2*b7-a6*a9*b4*b7-a6*a7*b6*b7+a6^2*b7^2-2*a8*(a9*b2*b4+a7*b2*b6-2*a6*b4*b6+2*a6*b2*b7)+3*a6*a7*b2*b8-2*a6^2*b4*b8+(-4*a5*a6*b2+a7^2*b2+2*a6^2*b3-a6*a7*b4)*b9+a3*(3*a6*b10*b2+2*a9*b2*b9-a6*b6*b9))+a4^2*(3*a6*b10^2-a9*b10*b9+a10*(-2*b10*b6+b9^2))-a4*(-a9^2*b10*b4-2*a7*a9*b10*b6+2*a6*a9*b10*b7+a10^2*(-2*b4*b6+2*b2*b7)-6*a6^2*b10*b8-a9^2*b6*b8+2*a6*(a7*b10+a9*b8)*b9+a10*(-3*a7*b10*b2+4*a6*b10*b4-2*a8*b6^2+a9*b6*b7-3*a9*b2*b8+4*a6*b6*b8+(4*a8*b2+a9*b4+a7*b6-4*a6*b7)*b9)+a8*(-3*a9*b10*b2+4*a6*b10*b6+a9*b6*b9-2*a6*b9^2))));
c3=(-a2^3*(3*b1*b10^2+6*b10*b5*b8+b8^3)+a2^2*(a4*b10^2*b3+2*a3*b10^2*b4+6*a8*b10*b2*b5+4*a9*b10*b4*b5+4*a9*b1*b10*b6-3*a8*b10*b3*b6+4*a7*b10*b5*b6-2*a9*b10*b3*b7-3*a8*b10*b4*b7+2*a6*b10*b5*b7-3*a5*b10*b6*b7-a7*b10*b7^2+6*a5*b10*b2*b8+2*a6*b10*b3*b8+4*a7*b10*b4*b8+4*a3*b10*b6*b8+4*a9*b5*b6*b8+2*a4*b10*b7*b8-3*a8*b6*b7*b8-a9*b7^2*b8+3*a8*b2*b8^2+2*a9*b4*b8^2+2*a7*b6*b8^2+a6*b7*b8^2+(2*a6*b1*b10-2*a7*b10*b3-3*a5*b10*b4+2*a4*b10*b5-3*a8*b5*b6-2*a3*b10*b7-2*a9*b5*b7+3*a8*b7^2-(2*a9*b3+3*a8*b4-2*a6*b5+3*a5*b6+2*a7*b7)*b8+a4*b8^2)*b9-(a9*b1-3*a8*b3+a7*b5-3*a5*b7+a3*b8)*b9^2+a1*(3*b10^2*b2-3*b10*b6*b9+b9^3)+a10*(-3*b5*b6*b7+b7^3+6*b2*b5*b8+b1*(6*b10*b2-3*b6*b9)-3*(b10*b3*b4+b3*b6*b8+b4*b7*b8+b4*b5*b9-2*b3*b7*b9)))+a2*(-3*a10^2*b1*b2^2-6*a1*a10*b10*b2^2-6*a5*a8*b10*b2^2+a10*a4*b10*b2*b3+4*a7*a9*b10*b2*b3-a4^2*b10^2*b4-a7*a8*b10*b2*b4-a5*a9*b10*b2*b4+3*a10^2*b2*b3*b4+2*a10*a4*b10*b4^2-2*a7*a9*b10*b4^2-a10^2*b4^3-6*a4*a9*b10*b2*b5-6*a10*a8*b2^2*b5-a10*a9*b2*b4*b5-a10*a9*b1*b2*b6-a5*a7*b10*b2*b6-a1*a9*b10*b2*b6+a4*a9*b10*b3*b6+6*a10*a8*b2*b3*b6-2*a7^2*b10*b4*b6+4*a4*a8*b10*b4*b6+2*a10*a9*b3*b4*b6-6*a10*a8*b4^2*b6-a10*a7*b2*b5*b6-a8*a9*b2*b5*b6-2*a9^2*b4*b5*b6-a9^2*b1*b6^2+2*a4*a5*b10*b6^2+a10*a7*b3*b6^2+a8*a9*b3*b6^2-6*a10*a5*b4*b6^2-3*a8^2*b4*b6^2+2*a10*a4*b5*b6^2-2*a7*a9*b5*b6^2-2*a1*a10*b6^3-2*a5*a8*b6^3+2*a7^2*b10*b2*b7+a4*a8*b10*b2*b7-4*a10*a9*b2*b3*b7+a4*a9*b10*b4*b7+6*a10*a8*b2*b4*b7+a10*a9*b4^2*b7+2*a9^2*b2*b5*b7+a4*a7*b10*b6*b7+6*a10*a5*b2*b6*b7+3*a8^2*b2*b6*b7+2*a10*a7*b4*b6*b7+2*a8*a9*b4*b6*b7+a7*a8*b6^2*b7+a5*a9*b6^2*b7-2*a10*a7*b2*b7^2-2*a8*a9*b2*b7^2-a10*a4*b6*b7^2-6*a4*a7*b10*b2*b8-6*a10*a5*b2^2*b8-3*a8^2*b2^2*b8+2*a9^2*b2*b3*b8-a10*a7*b2*b4*b8-a8*a9*b2*b4*b8-a9^2*b4^2*b8-2*a4^2*b10*b6*b8-a7*a8*b2*b6*b8-a5*a9*b2*b6*b8+4*a10*a4*b4*b6*b8-4*a7*a9*b4*b6*b8-a7^2*b6^2*b8+2*a4*a8*b6^2*b8+a10*a4*b2*b7*b8+4*a7*a9*b2*b7*b8+a4*a9*b6*b7*b8-3*a4*a9*b2*b8^2-a6^2*(2*b10*(b4*b5+b1*b6)+b8*(2*b5*b6+b4*b8))-a3*(3*a4*b10^2*b2+b10*b6*(a8*b2+4*a9*b4+2*a7*b6)-4*a9*b10*b2*b7+2*a9*b6^2*b8+a10*(b10*b2*b4-b6^2*b7+b2*b6*b8))+(2*a9^2*b1*b2+4*a3*a7*b10*b2-4*a10*a7*b2*b3+6*a10*a5*b2*b4+3*a8^2*b2*b4+a10*a7*b4^2+6*a1*a10*b2*b6+6*a5*a8*b2*b6+2*a10*a3*b4*b6+2*a7*a8*b4*b6+a5*a7*b6^2+a3*a8*b6^2-4*a10*a3*b2*b7-4*a7*a8*b2*b7+2*a7^2*b2*b8+a9*(a8*(-4*b2*b3+b4^2)+4*a7*b2*b5+2*a5*b4*b6+a4*b5*b6+a1*b6^2-4*a5*b2*b7+4*a3*b2*b8+a4*b4*b8)+a4*(a5*b10*b2+a7*b10*b4+a10*b2*b5+a3*b10*b6-2*a10*b3*b6-2*a10*b4*b7-2*a8*b6*b7+a8*b2*b8+a7*b6*b8))*b9-(2*a5*a7*b2+2*a3*a8*b2+2*a1*a9*b2+a4*a8*b4+a4*a5*b6)*b9^2+a6*(-6*a7*b10*b2*b5+a7*b10*b3*b6+4*a5*b10*b4*b6-4*a4*b10*b5*b6+4*a10*b4*b5*b6+2*a10*b1*b6^2+2*a1*b10*b6^2+a5*b10*b2*b7+a7*b10*b4*b7+a10*b2*b5*b7+a3*b10*b6*b7-2*a10*b3*b6*b7-a10*b4*b7^2-6*a3*b10*b2*b8+a10*b2*b3*b8-4*a4*b10*b4*b8+2*a10*b4^2*b8+2*a5*b6^2*b8+a7*b6*b7*b8-3*a7*b2*b8^2-2*a4*b6*b8^2+(a10*b1*b2+a1*b10*b2+a3*b10*b4-2*a10*b3*b4+a7*b5*b6-2*a5*b6*b7+a5*b2*b8+a7*b4*b8+a3*b6*b8)*b9-(a5*b4+a1*b6)*b9^2+a9*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9)+a8*(b10*b2*b3+2*b10*b4^2+2*b5*b6^2-b6*b7^2+4*b4*b6*b8+b2*b7*b8+b2*b5*b9-2*b3*b6*b9-2*b4*b7*b9)))+b2*(a4^3*b10^2+2*a6^3*(b1*b10+b5*b8)+a4*(3*a7*a8*b10*b2+3*a5*a9*b10*b2+2*a7*a9*b10*b4+a10^2*(-2*b2*b3+b4^2)+a7^2*b10*b6+2*a3*a9*b10*b6+a9^2*b5*b6+a8^2*b6^2-a8*a9*b6*b7+3*a8*a9*b2*b8+a9^2*b4*b8+2*a7*a9*b6*b8-(2*a8^2*b2+a8*a9*b4+a7*a8*b6+a5*a9*b6)*b9+a10*(3*a3*b10*b2+3*a9*b2*b5-a9*b3*b6+4*a8*b4*b6+2*a5*b6^2-4*a8*b2*b7-a9*b4*b7-a7*b6*b7+3*a7*b2*b8-(4*a5*b2+a7*b4+a3*b6)*b9))-a4^2*(2*a8*b10*b6+a9*b10*b7+a7*b10*b9+a9*b8*b9-a8*b9^2+2*a10*(b10*b4+b6*b8-b7*b9))-a6^2*(-6*a4*b10*b5+2*a10*b4*b5+2*a10*b1*b6+2*a1*b10*b6+2*a8*b5*b6+a3*b10*b7-2*a10*b3*b7+a9*b5*b7-a8*b7^2+a9*b3*b8+2*a8*b4*b8-3*a4*b8^2+a9*b1*b9-2*a8*b3*b9+a3*b8*b9-a1*b9^2+a7*(b10*b3+b7*b8+b5*b9)+2*a5*(b10*b4+b6*b8-b7*b9))+a6*(3*a3*a8*b10*b2+3*a1*a9*b10*b2-2*a4*a9*b10*b3+a7^2*b10*b4-4*a4*a8*b10*b4+2*a3*a9*b10*b4+3*a8*a9*b2*b5+a9^2*b4*b5+a9^2*b1*b6+2*a3*a7*b10*b6-a8*a9*b3*b6+2*a8^2*b4*b6+2*a7*a9*b5*b6-2*a4*a7*b10*b7-2*a8^2*b2*b7-a8*a9*b4*b7-a7*a8*b6*b7+6*a4^2*b10*b8+3*a7*a8*b2*b8+2*a7*a9*b4*b8+a7^2*b6*b8-4*a4*a8*b6*b8+2*a3*a9*b6*b8-2*a4*a9*b7*b8-(2*a3*a4*b10+a7*a8*b4+2*a4*a9*b5+a3*a8*b6+a1*a9*b6-4*a4*a8*b7+2*a4*a7*b8)*b9+a10*(3*a9*b1*b2-4*a8*b2*b3-a9*b3*b4+2*a8*b4^2+3*a7*b2*b5-a7*b3*b6+4*a5*b4*b6-4*a4*b5*b6+2*a1*b6^2-4*a5*b2*b7-a7*b4*b7-a3*b6*b7+2*a4*b7^2+3*a3*b2*b8-4*a4*b4*b8-4*a1*b2*b9+4*a4*b3*b9-a3*b4*b9)+a5*(3*a7*b10*b2-4*a4*b10*b6+2*a8*b6^2-a9*b6*b7+3*a9*b2*b8-(4*a8*b2+a9*b4+a7*b6)*b9+2*a4*b9^2))+b2*(-a9^3*b1-a7^3*b10+(3*a1*a10^2+6*a10*a5*a8+a8^3)*b2-a10^2*a3*b4+a9^2*(a8*b3-3*a7*b5+a5*b7-3*a3*b8+a1*b9)+a7*a8*(-a8*b6+a7*b9)+a10*(-2*a3*a8*b6+a7^2*b7-2*a7*(a8*b4+a5*b6-a3*b9))-a9*(a8^2*b4+2*a5*a8*b6+2*a10*(-a7*b3+a5*b4+a1*b6)-2*a7*a8*b7+3*a7^2*b8-2*a5*a7*b9-2*a3*(-3*a7*b10+a10*b7+a8*b9)))));
c4=(2*a4^3*b10*b2*b8-3*a2^3*(b5*b8^2+b10*(b5^2+2*b1*b8))-a4^2*(-6*a6*b10*b2*b5+2*a5*b10*b2*b6+2*a2*b10*b5*b6+2*a10*b2*b5*b6+a7*b10*b2*b7-a10*b2*b7^2+2*a2*b10*b4*b8+2*a10*b2*b4*b8-3*a6*b2*b8^2+a2*b6*b8^2+b2*(a3*b10-2*a10*b3+a7*b8)*b9-a5*b2*b9^2+a9*b2*(b10*b3+b7*b8+b5*b9)+2*a8*b2*(b10*b4+b6*b8-b7*b9))+a2^2*(6*a5*b10*b2*b5+2*a6*b10*b3*b5+4*a7*b10*b4*b5+3*a10*b2*b5^2+4*a7*b1*b10*b6-3*a5*b10*b3*b6+4*a3*b10*b5*b6-3*a10*b3*b5*b6+2*a6*b1*b10*b7-2*a7*b10*b3*b7-3*a5*b10*b4*b7-3*a10*b4*b5*b7-3*a10*b1*b6*b7-3*a1*b10*b6*b7-a3*b10*b7^2+3*a10*b3*b7^2+6*a10*b1*b2*b8+6*a1*b10*b2*b8+4*a3*b10*b4*b8-3*a10*b3*b4*b8+4*a7*b5*b6*b8+2*a6*b5*b7*b8-3*a5*b6*b7*b8-a7*b7^2*b8+3*a5*b2*b8^2+a6*b3*b8^2+2*a7*b4*b8^2+2*a3*b6*b8^2-(2*a3*b10*b3-3*a10*b3^2+3*a10*b1*b4+3*a1*b10*b4-a6*b5^2+3*a5*b5*b6+2*a7*b5*b7-3*a5*b7^2-2*a6*b1*b8+2*a7*b3*b8+3*a5*b4*b8+3*a1*b6*b8+2*a3*b7*b8)*b9-(a7*b1-3*a5*b3+a3*b5-3*a1*b7)*b9^2-a9*(b10*b3^2-4*b1*b10*b4+b5*(-2*b5*b6+b7^2)-4*b4*b5*b8-4*b1*b6*b8+2*b3*b7*b8+2*b3*b5*b9+2*b1*b7*b9)+a8*(-3*b5*b6*b7+b7^3+6*b2*b5*b8+b1*(6*b10*b2-3*b6*b9)-3*(b10*b3*b4+b3*b6*b8+b4*b7*b8+b4*b5*b9-2*b3*b7*b9)))+a2*(-3*a5^2*b10*b2^2-6*a1*a8*b10*b2^2+2*a7^2*b10*b2*b3+4*a3*a9*b10*b2*b3-a5*a7*b10*b2*b4-a3*a8*b10*b2*b4-a1*a9*b10*b2*b4-a7^2*b10*b4^2-2*a3*a9*b10*b4^2-3*a8^2*b2^2*b5+2*a9^2*b2*b3*b5-a8*a9*b2*b4*b5-a9^2*b4^2*b5-a8*a9*b1*b2*b6-a3*a5*b10*b2*b6-a1*a7*b10*b2*b6+3*a8^2*b2*b3*b6-2*a9^2*b1*b4*b6-4*a3*a7*b10*b4*b6+2*a8*a9*b3*b4*b6-3*a8^2*b4^2*b6-a7*a8*b2*b5*b6-a5*a9*b2*b5*b6-4*a7*a9*b4*b5*b6-2*a7*a9*b1*b6^2-a3^2*b10*b6^2+a7*a8*b3*b6^2+a5*a9*b3*b6^2-6*a5*a8*b4*b6^2-a7^2*b5*b6^2-2*a3*a9*b5*b6^2-a5^2*b6^3-2*a1*a8*b6^3+2*a9^2*b1*b2*b7+4*a3*a7*b10*b2*b7-4*a8*a9*b2*b3*b7+3*a8^2*b2*b4*b7+a8*a9*b4^2*b7+4*a7*a9*b2*b5*b7+6*a5*a8*b2*b6*b7+2*a7*a8*b4*b6*b7+2*a5*a9*b4*b6*b7+a5*a7*b6^2*b7+a3*a8*b6^2*b7+a1*a9*b6^2*b7-2*a7*a8*b2*b7^2-2*a5*a9*b2*b7^2-6*a5*a8*b2^2*b8+4*a7*a9*b2*b3*b8-a7*a8*b2*b4*b8-a5*a9*b2*b4*b8-2*a7*a9*b4^2*b8-a5*a7*b2*b6*b8-a3*a8*b2*b6*b8-a1*a9*b2*b6*b8-2*a7^2*b4*b6*b8-4*a3*a9*b4*b6*b8-2*a3*a7*b6^2*b8+2*a7^2*b2*b7*b8+4*a3*a9*b2*b7*b8-a6^2*(b5*(b5*b6+2*b4*b8)+2*b1*(b10*b4+b6*b8))+(2*a3^2*b10*b2-4*a5*a9*b2*b3+6*a5*a8*b2*b4+a5*a9*b4^2+2*a7^2*b2*b5+4*a3*a9*b2*b5+3*a5^2*b2*b6+6*a1*a8*b2*b6+2*a3*a8*b4*b6+2*a1*a9*b4*b6+a3*a5*b6^2-4*(a3*a8+a1*a9)*b2*b7+a7*(4*a9*b1*b2-4*a8*b2*b3+a8*b4^2+2*a5*b4*b6+a1*b6^2-4*a5*b2*b7+4*a3*b2*b8))*b9-2*(a3*a5+a1*a7)*b2*b9^2-a10*(a9*(2*b2*b3^2+b1*b2*b4-b3*b4^2)+2*a8*(3*b1*b2^2-3*b2*b3*b4+b4^3)+6*a5*b2^2*b5+a7*b2*b4*b5+a7*b1*b2*b6-6*a5*b2*b3*b6-2*a7*b3*b4*b6+6*a5*b4^2*b6+a3*b2*b5*b6-a3*b3*b6^2+6*a1*b4*b6^2+4*a7*b2*b3*b7-6*a5*b2*b4*b7-a7*b4^2*b7-6*a1*b2*b6*b7-2*a3*b4*b6*b7+2*a3*b2*b7^2+6*a1*b2^2*b8+a3*b2*b4*b8+4*a3*b2*b3*b9-6*a1*b2*b4*b9-a3*b4^2*b9)+a6*(-6*a3*b10*b2*b5+a10*b2*b3*b5+2*a10*b4^2*b5-3*a9*b2*b5^2+a3*b10*b3*b6-a10*b3^2*b6+4*a10*b1*b4*b6+4*a1*b10*b4*b6+a9*b3*b5*b6+4*a8*b4*b5*b6+2*a8*b1*b6^2+a10*b1*b2*b7+a1*b10*b2*b7+a3*b10*b4*b7-2*a10*b3*b4*b7+a8*b2*b5*b7+a9*b4*b5*b7+a9*b1*b6*b7-2*a8*b3*b6*b7-a8*b4*b7^2-6*a9*b1*b2*b8+a8*b2*b3*b8+a9*b3*b4*b8+2*a8*b4^2*b8+2*a1*b6^2*b8+a3*b6*b7*b8-3*a3*b2*b8^2+(a8*b1*b2+a9*b1*b4-2*a8*b3*b4+a3*b5*b6-2*a1*b6*b7+a1*b2*b8+a3*b4*b8)*b9-a1*b4*b9^2+a7*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9)+a5*(b10*b2*b3+2*b10*b4^2+2*b5*b6^2-b6*b7^2+4*b4*b6*b8+b2*b7*b8+b2*b5*b9-2*b3*b6*b9-2*b4*b7*b9)))+a4*(6*a6^2*b2*(b1*b10+b5*b8)+b2*(3*a5*a7*b10*b2+3*a3*a8*b10*b2+3*a1*a9*b10*b2+a7^2*b10*b4+2*a3*a9*b10*b4+3*a8*a9*b2*b5+a9^2*b4*b5+a9^2*b1*b6+2*a3*a7*b10*b6-a8*a9*b3*b6+2*a8^2*b4*b6+2*a7*a9*b5*b6+2*a5*a8*b6^2-2*a8^2*b2*b7-a8*a9*b4*b7-a7*a8*b6*b7-a5*a9*b6*b7+3*a7*a8*b2*b8+3*a5*a9*b2*b8+2*a7*a9*b4*b8+a7^2*b6*b8+2*a3*a9*b6*b8-(a7*a8*b4+a3*a8*b6+a1*a9*b6+a5*(4*a8*b2+a9*b4+a7*b6))*b9+a10*(3*a9*b1*b2-4*a8*b2*b3-a9*b3*b4+2*a8*b4^2+3*a7*b2*b5-a7*b3*b6+4*a5*b4*b6+2*a1*b6^2-4*a5*b2*b7-a7*b4*b7-a3*b6*b7+3*a3*b2*b8-4*a1*b2*b9-a3*b4*b9))+a2^2*(2*b10*(b5*b7+b3*b8+b1*b9)+b8*(b7*b8+2*b5*b9))-2*a6*(2*a2*b10*b4*b5+2*a10*b2*b4*b5+2*a2*b1*b10*b6+2*a10*b1*b2*b6+2*a1*b10*b2*b6+2*a8*b2*b5*b6+a3*b10*b2*b7-2*a10*b2*b3*b7+a9*b2*b5*b7-a8*b2*b7^2+a9*b2*b3*b8+2*a8*b2*b4*b8+2*a2*b5*b6*b8+a2*b4*b8^2+a9*b1*b2*b9-2*a8*b2*b3*b9+a3*b2*b8*b9-a1*b2*b9^2+a7*b2*(b10*b3+b7*b8+b5*b9)+2*a5*b2*(b10*b4+b6*b8-b7*b9))+a2*(-6*a7*b10*b2*b5+a7*b10*b3*b6+4*a5*b10*b4*b6+4*a10*b4*b5*b6+2*a10*b1*b6^2+2*a1*b10*b6^2+a5*b10*b2*b7+a7*b10*b4*b7+a10*b2*b5*b7+a3*b10*b6*b7-2*a10*b3*b6*b7-a10*b4*b7^2-6*a3*b10*b2*b8+a10*b2*b3*b8+2*a10*b4^2*b8+2*a5*b6^2*b8+a7*b6*b7*b8-3*a7*b2*b8^2+(a10*b1*b2+a1*b10*b2+a3*b10*b4-2*a10*b3*b4+a7*b5*b6-2*a5*b6*b7+a5*b2*b8+a7*b4*b8+a3*b6*b8)*b9-(a5*b4+a1*b6)*b9^2+a9*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9)+a8*(b10*b2*b3+2*b10*b4^2+2*b5*b6^2-b6*b7^2+4*b4*b6*b8+b2*b7*b8+b2*b5*b9-2*b3*b6*b9-2*b4*b7*b9)))+b2*(a6^3*(b5^2+2*b1*b8)-a6*(-3*a1*a7*b10*b2-a9^2*b1*b4+a8^2*(2*b2*b3-b4^2)-3*a5*a9*b2*b5-2*a7*a9*b4*b5-2*a7*a9*b1*b6-a3^2*b10*b6+a5*a9*b3*b6-a7^2*b5*b6-a5^2*b6^2+a5*a9*b4*b7+a5*a7*b6*b7+a1*a9*b6*b7-3*a5*a7*b2*b8-3*a1*a9*b2*b8-a7^2*b4*b8+(2*a5^2*b2+a5*a7*b4+a1*a9*b4+a1*a7*b6)*b9+a8*(-3*a9*b1*b2+a9*b3*b4-3*a7*b2*b5+a7*b3*b6-4*a5*b4*b6-2*a1*b6^2+4*a5*b2*b7+a7*b4*b7+a3*b6*b7-3*a3*b2*b8+4*a1*b2*b9+a3*b4*b9)-a3*(3*a5*b10*b2+2*a7*b10*b4+2*a9*b5*b6+2*a9*b4*b8+2*a7*b6*b8-a5*b6*b9))-a6^2*(a9*b3*b5+2*a8*b4*b5+2*a8*b1*b6+2*a5*b5*b6+a9*b1*b7-2*a8*b3*b7+a7*b5*b7-a5*b7^2+a7*b3*b8+2*a5*b4*b8+a7*b1*b9-2*a5*b3*b9+a3*(b10*b3+b7*b8+b5*b9)+2*a1*(b10*b4+b6*b8-b7*b9))+b2*(-3*a3^2*a9*b10+3*a5*a8^2*b2-a3*a8^2*b6+a9*(a5*a9*b3-2*a5*a8*b4-a5^2*b6-2*a1*a8*b6+a1*a9*b7)-a7^3*b8+a7^2*(-3*a3*b10-3*a9*b5+a8*b7+a5*b9)+a3*a9*(-3*a9*b5+2*a8*b7+2*a5*b9)+a7*(-3*a9^2*b1+2*a9*(a8*b3+a5*b7-3*a3*b8+a1*b9)-a8*(a8*b4+2*a5*b6-2*a3*b9)))+a10*(a6^2*(b3^2-2*b1*b4)+a6*(3*a7*b1*b2-4*a5*b2*b3-a7*b3*b4+2*a5*b4^2+3*a3*b2*b5-a3*b3*b6+4*a1*b4*b6-4*a1*b2*b7-a3*b4*b7)+b2*(3*a5^2*b2+(a7^2+2*a3*a9)*b3-2*a5*(a7*b4+a3*b6)-2*a1*(-3*a8*b2+a9*b4+a7*b6)+a3*(-2*a8*b4+2*a7*b7+a3*b9)))));
c5=(-3*a2^3*(b5^2*b8+b1*(2*b10*b5+b8^2))+b2*(3*a10*a3*a6*b1*b2+3*a6*a7*a8*b1*b2+3*a5*a6*a9*b1*b2-3*a7^2*a9*b1*b2-3*a3*a9^2*b1*b2+3*a1*a3*a6*b10*b2-3*a3^2*a7*b10*b2+6*a1*a10*a5*b2^2+3*a5^2*a8*b2^2+3*a1*a8^2*b2^2-a6^2*a9*b1*b3-4*a1*a10*a6*b2*b3+2*a10*a3*a7*b2*b3-4*a5*a6*a8*b2*b3+a7^2*a8*b2*b3+2*a5*a7*a9*b2*b3+2*a3*a8*a9*b2*b3+a1*a9^2*b2*b3+a6^2*a8*b3^2-2*a6^2*a8*b1*b4+2*a6*a7*a9*b1*b4+a3^2*a6*b10*b4-2*a10*a3*a5*b2*b4-2*a1*a10*a7*b2*b4-2*a5*a7*a8*b2*b4-a3*a8^2*b2*b4-a5^2*a9*b2*b4-2*a1*a8*a9*b2*b4-a10*a3*a6*b3*b4-a6*a7*a8*b3*b4-a5*a6*a9*b3*b4+2*a1*a10*a6*b4^2+2*a5*a6*a8*b4^2+2*a6^3*b1*b5+3*a5*a6*a7*b2*b5-a7^3*b2*b5+3*a3*a6*a8*b2*b5+3*a1*a6*a9*b2*b5-6*a3*a7*a9*b2*b5-a6^2*a7*b3*b5-2*a5*a6^2*b4*b5+a6*a7^2*b4*b5+2*a3*a6*a9*b4*b5-2*a5*a6^2*b1*b6+a6*a7^2*b1*b6+2*a3*a6*a9*b1*b6-2*a1*a10*a3*b2*b6-a5^2*a7*b2*b6-2*a3*a5*a8*b2*b6-2*a1*a7*a8*b2*b6-2*a1*a5*a9*b2*b6-a5*a6*a7*b3*b6-a3*a6*a8*b3*b6-a1*a6*a9*b3*b6+2*a5^2*a6*b4*b6+4*a1*a6*a8*b4*b6-2*a1*a6^2*b5*b6+2*a3*a6*a7*b5*b6+2*a1*a5*a6*b6^2-a6^2*a7*b1*b7+a10*a3^2*b2*b7-2*a5^2*a6*b2*b7+a5*a7^2*b2*b7-4*a1*a6*a8*b2*b7+2*a3*a7*a8*b2*b7+2*a3*a5*a9*b2*b7+2*a1*a7*a9*b2*b7+2*a5*a6^2*b3*b7-a5*a6*a7*b4*b7-a3*a6*a8*b4*b7-a1*a6*a9*b4*b7-a3*a6^2*b5*b7-a3*a5*a6*b6*b7-a1*a6*a7*b6*b7+a1*a6^2*b7^2+3*a3*a5*a6*b2*b8+3*a1*a6*a7*b2*b8-3*a3*a7^2*b2*b8-3*a3^2*a9*b2*b8-a3*a6^2*b3*b8-2*a1*a6^2*b4*b8+2*a3*a6*a7*b4*b8+a3^2*a6*b6*b8+a4^3*(2*b10*b5+b8^2)+(a3^2*a8*b2+2*a3*(a5*a7+a1*a9)*b2+a1*(-4*a5*a6*b2+a7^2*b2+2*a6^2*b3-a6*a7*b4)-a3*a6*(a6*b1+a5*b4+a1*b6))*b9-a4^2*(a7*b10*b3+2*a5*b10*b4+2*a10*b4*b5+2*a10*b1*b6+2*a1*b10*b6+2*a8*b5*b6+a3*b10*b7-2*a10*b3*b7+a9*b5*b7-a8*b7^2+a9*b3*b8+2*a8*b4*b8+2*a5*b6*b8+a7*b7*b8-6*a6*(b1*b10+b5*b8)+(a9*b1-2*a8*b3+a7*b5-2*a5*b7+a3*b8)*b9-a1*b9^2)+a4*(3*a3*a5*b10*b2+3*a1*a7*b10*b2-2*a3*a6*b10*b3+a9^2*b1*b4-4*a1*a6*b10*b4+2*a3*a7*b10*b4+a8^2*(-2*b2*b3+b4^2)+3*a5*a9*b2*b5-2*a6*a9*b3*b5+2*a7*a9*b4*b5+3*a6^2*b5^2+2*a7*a9*b1*b6+a3^2*b10*b6-a5*a9*b3*b6-4*a5*a6*b5*b6+a7^2*b5*b6+2*a3*a9*b5*b6+a5^2*b6^2-2*a6*a9*b1*b7-a5*a9*b4*b7-2*a6*a7*b5*b7-a5*a7*b6*b7-a1*a9*b6*b7+2*a5*a6*b7^2+a10*(3*a7*b1*b2-4*a5*b2*b3+2*a6*b3^2-4*a6*b1*b4-a7*b3*b4+2*a5*b4^2+3*a3*b2*b5-a3*b3*b6+4*a1*b4*b6-4*a1*b2*b7-a3*b4*b7)+6*a6^2*b1*b8+3*a5*a7*b2*b8+3*a1*a9*b2*b8-2*a6*a7*b3*b8-4*a5*a6*b4*b8+a7^2*b4*b8+2*a3*a9*b4*b8-4*a1*a6*b6*b8+2*a3*a7*b6*b8-2*a3*a6*b7*b8-(2*a5^2*b2+a5*a7*b4+a1*a9*b4+a3*a5*b6+a1*a7*b6+2*a6*(a7*b1-2*a5*b3+a3*b5-2*a1*b7))*b9+a8*(3*a9*b1*b2-a9*b3*b4+3*a7*b2*b5-4*a6*b4*b5-4*a6*b1*b6-a7*b3*b6+4*a5*b4*b6+2*a1*b6^2-4*a5*b2*b7+4*a6*b3*b7-a7*b4*b7-a3*b6*b7+3*a3*b2*b8-4*a1*b2*b9-a3*b4*b9)))+a2^2*(-a7*b10*b3^2+4*a7*b1*b10*b4+6*a10*b1*b2*b5+6*a1*b10*b2*b5+2*a4*b10*b3*b5+4*a3*b10*b4*b5-3*a10*b3*b4*b5+3*a8*b2*b5^2+2*a9*b4*b5^2+4*a3*b1*b10*b6-3*a10*b1*b3*b6-3*a1*b10*b3*b6+4*a9*b1*b5*b6-3*a8*b3*b5*b6+2*a7*b5^2*b6+2*a4*b1*b10*b7-2*a3*b10*b3*b7+3*a10*b3^2*b7-3*a10*b1*b4*b7-3*a1*b10*b4*b7-2*a9*b3*b5*b7-3*a8*b4*b5*b7+a6*b5^2*b7-3*a8*b1*b6*b7-a9*b1*b7^2+3*a8*b3*b7^2-a7*b5*b7^2+6*a8*b1*b2*b8-a9*b3^2*b8+4*a9*b1*b4*b8-3*a8*b3*b4*b8+4*a7*b4*b5*b8+4*a7*b1*b6*b8+4*a3*b5*b6*b8-2*a7*b3*b7*b8+2*a4*b5*b7*b8-3*a1*b6*b7*b8-a3*b7^2*b8+3*a1*b2*b8^2+a4*b3*b8^2+2*a3*b4*b8^2-(2*a9*b1*b3-3*a8*b3^2+3*a8*b1*b4+2*a7*b3*b5-a4*b5^2+3*a1*b5*b6+2*a7*b1*b7+2*a3*b5*b7-3*a1*b7^2-2*a4*b1*b8+2*a3*b3*b8+3*a1*b4*b8)*b9+(-a3*b1+3*a1*b3)*b9^2+2*a6*(b3*b5*b8+b1*(b10*b3+b7*b8+b5*b9))+a5*(-3*b5*b6*b7+b7^3+6*b2*b5*b8+b1*(6*b10*b2-3*b6*b9)-3*(b10*b3*b4+b3*b6*b8+b4*b7*b8+b4*b5*b9-2*b3*b7*b9)))+a2*(-6*a10*a5*b1*b2^2-3*a8^2*b1*b2^2-6*a1*a5*b10*b2^2+a10*a6*b1*b2*b3+2*a9^2*b1*b2*b3+a1*a6*b10*b2*b3-2*a10*a7*b2*b3^2-2*a8*a9*b2*b3^2-a10*a7*b1*b2*b4-a8*a9*b1*b2*b4-a1*a7*b10*b2*b4+6*a10*a5*b2*b3*b4+3*a8^2*b2*b3*b4-a10*a6*b3^2*b4+2*a10*a6*b1*b4^2-a9^2*b1*b4^2+2*a1*a6*b10*b4^2+a10*a7*b3*b4^2+a8*a9*b3*b4^2-2*a10*a5*b4^3-a8^2*b4^3-6*a6*a9*b1*b2*b5-6*a1*a10*b2^2*b5-6*a5*a8*b2^2*b5+a6*a8*b2*b3*b5+4*a7*a9*b2*b3*b5-a7*a8*b2*b4*b5-a5*a9*b2*b4*b5+a6*a9*b3*b4*b5+2*a6*a8*b4^2*b5-2*a7*a9*b4^2*b5-3*a6*a7*b2*b5^2-a6^2*b4*b5^2-a7*a8*b1*b2*b6-a5*a9*b1*b2*b6+a6*a9*b1*b3*b6+6*a1*a10*b2*b3*b6+6*a5*a8*b2*b3*b6-a6*a8*b3^2*b6+4*a6*a8*b1*b4*b6-4*a7*a9*b1*b4*b6+2*a7*a8*b3*b4*b6+2*a5*a9*b3*b4*b6-6*a1*a10*b4^2*b6-6*a5*a8*b4^2*b6-2*a6^2*b1*b5*b6-a5*a7*b2*b5*b6-a1*a9*b2*b5*b6+a6*a7*b3*b5*b6+4*a5*a6*b4*b5*b6-2*a7^2*b4*b5*b6+2*a5*a6*b1*b6^2-a7^2*b1*b6^2+a5*a7*b3*b6^2+a1*a9*b3*b6^2-3*a5^2*b4*b6^2-6*a1*a8*b4*b6^2+2*a1*a6*b5*b6^2-2*a1*a5*b6^3+a6*a8*b1*b2*b7+4*a7*a9*b1*b2*b7-4*a7*a8*b2*b3*b7-4*a5*a9*b2*b3*b7+a6*a9*b1*b4*b7+6*a1*a10*b2*b4*b7+6*a5*a8*b2*b4*b7-2*a6*a8*b3*b4*b7+a7*a8*b4^2*b7+a5*a9*b4^2*b7+a5*a6*b2*b5*b7+2*a7^2*b2*b5*b7+a6*a7*b4*b5*b7+a6*a7*b1*b6*b7+3*a5^2*b2*b6*b7+6*a1*a8*b2*b6*b7-2*a5*a6*b3*b6*b7+2*a5*a7*b4*b6*b7+2*a1*a9*b4*b6*b7+a1*a7*b6^2*b7-2*a5*a7*b2*b7^2-2*a1*a9*b2*b7^2-a5*a6*b4*b7^2-a1*a6*b6*b7^2-6*a6*a7*b1*b2*b8-3*a5^2*b2^2*b8-6*a1*a8*b2^2*b8+a5*a6*b2*b3*b8+2*a7^2*b2*b3*b8-2*a6^2*b1*b4*b8-a5*a7*b2*b4*b8-a1*a9*b2*b4*b8+a6*a7*b3*b4*b8+2*a5*a6*b4^2*b8-a7^2*b4^2*b8-a1*a7*b2*b6*b8+4*a1*a6*b4*b6*b8+a1*a6*b2*b7*b8-a4^2*(2*b10*(b4*b5+b1*b6)+b8*(2*b5*b6+b4*b8))+(2*a7^2*b1*b2+3*a5^2*b2*b4+a5*(a6*b1*b2-4*a7*b2*b3-2*a6*b3*b4+a7*b4^2+6*a1*b2*b6)+a7*(a6*b1*b4+2*a1*b4*b6-4*a1*b2*b7)+a1*(-4*a9*b2*b3+6*a8*b2*b4+a9*b4^2+a6*b2*b5-2*a6*b3*b6-2*a6*b4*b7))*b9+a3^2*(-2*b10*b4*b6+2*b10*b2*b7-b6^2*b8+2*b2*b8*b9)+a3*(-a5*b10*b2*b4-6*a4*b10*b2*b5-a10*b2*b4*b5-a10*b1*b2*b6-a1*b10*b2*b6+a4*b10*b3*b6+2*a10*b3*b4*b6-a8*b2*b5*b6-4*a9*b4*b5*b6-2*a9*b1*b6^2+a8*b3*b6^2-4*a10*b2*b3*b7+a4*b10*b4*b7+a10*b4^2*b7+4*a9*b2*b5*b7+2*a8*b4*b6*b7+a5*b6^2*b7-2*a8*b2*b7^2+4*a9*b2*b3*b8-a8*b2*b4*b8-2*a9*b4^2*b8-a5*b2*b6*b8+a4*b6*b7*b8-3*a4*b2*b8^2+(4*a9*b1*b2-4*a8*b2*b3+a8*b4^2+2*a5*b4*b6+a4*b5*b6+a1*b6^2-4*a5*b2*b7+a4*b4*b8)*b9-2*a1*b2*b9^2+a7*(4*b10*b2*b3-2*b10*b4^2-2*b5*b6^2-4*b4*b6*b8+4*b2*b7*b8+4*b2*b5*b9)+a6*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9))+a4*(-4*a6*b1*b10*b4+a10*b2*b3*b5+2*a10*b4^2*b5-3*a9*b2*b5^2-a10*b3^2*b6+4*a10*b1*b4*b6+4*a1*b10*b4*b6+a9*b3*b5*b6+4*a8*b4*b5*b6-2*a6*b5^2*b6+2*a8*b1*b6^2+a10*b1*b2*b7+a1*b10*b2*b7-2*a10*b3*b4*b7+a8*b2*b5*b7+a9*b4*b5*b7+a9*b1*b6*b7-2*a8*b3*b6*b7-a8*b4*b7^2-6*a9*b1*b2*b8+a8*b2*b3*b8+a9*b3*b4*b8+2*a8*b4^2*b8-4*a6*b4*b5*b8-4*a6*b1*b6*b8+2*a1*b6^2*b8+(a8*b1*b2+a9*b1*b4-2*a8*b3*b4-2*a1*b6*b7+a1*b2*b8)*b9-a1*b4*b9^2+a7*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9)+a5*(b10*b2*b3+2*b10*b4^2+2*b5*b6^2-b6*b7^2+4*b4*b6*b8+b2*b7*b8+b2*b5*b9-2*b3*b6*b9-2*b4*b7*b9))));
c6=(-a2^3*(3*b1^2*b10+b5^3+6*b1*b5*b8)+a2*(-6*a1*a10*b1*b2^2-6*a5*a8*b1*b2^2-3*a1^2*b10*b2^2+a10*a4*b1*b2*b3+4*a7*a9*b1*b2*b3+a1*a4*b10*b2*b3-2*a7*a8*b2*b3^2-2*a5*a9*b2*b3^2-2*a4^2*b1*b10*b4-a7*a8*b1*b2*b4-a5*a9*b1*b2*b4+6*a1*a10*b2*b3*b4+6*a5*a8*b2*b3*b4-a10*a4*b3^2*b4+2*a10*a4*b1*b4^2-2*a7*a9*b1*b4^2+2*a1*a4*b10*b4^2+a7*a8*b3*b4^2+a5*a9*b3*b4^2-2*a1*a10*b4^3-2*a5*a8*b4^3-6*a4*a9*b1*b2*b5-3*a5^2*b2^2*b5-6*a1*a8*b2^2*b5+2*a7^2*b2*b3*b5+a4*a8*b2*b3*b5-a5*a7*b2*b4*b5-a1*a9*b2*b4*b5+a4*a9*b3*b4*b5-a7^2*b4^2*b5+2*a4*a8*b4^2*b5-3*a4*a7*b2*b5^2-a5*a7*b1*b2*b6-a1*a9*b1*b2*b6+a4*a9*b1*b3*b6+3*a5^2*b2*b3*b6+6*a1*a8*b2*b3*b6-a4*a8*b3^2*b6-2*a7^2*b1*b4*b6+4*a4*a8*b1*b4*b6+2*a5*a7*b3*b4*b6+2*a1*a9*b3*b4*b6-3*a5^2*b4^2*b6-6*a1*a8*b4^2*b6-a1*a7*b2*b5*b6+a4*a7*b3*b5*b6+4*a4*a5*b4*b5*b6-a4^2*b5^2*b6+2*a4*a5*b1*b6^2+a1*a7*b3*b6^2-6*a1*a5*b4*b6^2+2*a1*a4*b5*b6^2-a1^2*b6^3-a6^2*b1*(2*b4*b5+b1*b6)+2*a7^2*b1*b2*b7+a4*a8*b1*b2*b7-4*a5*a7*b2*b3*b7-4*a1*a9*b2*b3*b7+a4*a9*b1*b4*b7+3*a5^2*b2*b4*b7+6*a1*a8*b2*b4*b7-2*a4*a8*b3*b4*b7+a5*a7*b4^2*b7+a1*a9*b4^2*b7+a4*a5*b2*b5*b7+a4*a7*b4*b5*b7+a4*a7*b1*b6*b7+6*a1*a5*b2*b6*b7-2*a4*a5*b3*b6*b7+2*a1*a7*b4*b6*b7-2*a1*a7*b2*b7^2-a4*a5*b4*b7^2-a1*a4*b6*b7^2-6*a4*a7*b1*b2*b8-6*a1*a5*b2^2*b8+a4*a5*b2*b3*b8-a1*a7*b2*b4*b8+a4*a7*b3*b4*b8+2*a4*a5*b4^2*b8-2*a4^2*b4*b5*b8-2*a4^2*b1*b6*b8+4*a1*a4*b4*b6*b8+a1*a4*b2*b7*b8+(a1*(-4*a7*b2*b3+6*a5*b2*b4+a7*b4^2+3*a1*b2*b6)+a4*(a5*b1*b2+a7*b1*b4-2*a5*b3*b4+a1*b2*b5-2*a1*b3*b6-2*a1*b4*b7))*b9+a6*(a8*b1*b2*b3-a8*b3^2*b4+2*a8*b1*b4^2+a9*b1*(-3*b1*b2+b3*b4)-6*a7*b1*b2*b5+a5*b2*b3*b5+a7*b3*b4*b5+2*a5*b4^2*b5-3*a3*b2*b5^2-2*a4*b4*b5^2+a7*b1*b3*b6-a5*b3^2*b6+4*a5*b1*b4*b6-4*a4*b1*b5*b6+a3*b3*b5*b6+4*a1*b4*b5*b6+2*a1*b1*b6^2+a5*b1*b2*b7+a7*b1*b4*b7-2*a5*b3*b4*b7+a1*b2*b5*b7+a3*b4*b5*b7+a3*b1*b6*b7-2*a1*b3*b6*b7-a1*b4*b7^2-6*a3*b1*b2*b8+a1*b2*b3*b8-4*a4*b1*b4*b8+a3*b3*b4*b8+2*a1*b4^2*b8+a1*b1*b2*b9+a3*b1*b4*b9-2*a1*b3*b4*b9)+a3^2*(2*b10*b2*b3-b10*b4^2-b5*b6^2-2*b4*b6*b8+2*b2*b7*b8+2*b2*b5*b9)+a3*(-2*a10*b2*b3^2-a10*b1*b2*b4-a1*b10*b2*b4+a10*b3*b4^2+4*a9*b2*b3*b5-a8*b2*b4*b5-2*a9*b4^2*b5-a8*b1*b2*b6-4*a9*b1*b4*b6+2*a8*b3*b4*b6-a5*b2*b5*b6-4*a7*b4*b5*b6-2*a7*b1*b6^2+a5*b3*b6^2+4*a9*b1*b2*b7-4*a8*b2*b3*b7+a8*b4^2*b7+4*a7*b2*b5*b7+2*a5*b4*b6*b7+a1*b6^2*b7-2*a5*b2*b7^2+4*a7*b2*b3*b8-a5*b2*b4*b8-2*a7*b4^2*b8-a1*b2*b6*b8+(4*a7*b1*b2-4*a5*b2*b3+a5*b4^2+2*a1*b4*b6-4*a1*b2*b7)*b9+a4*(-6*b1*b10*b2+b10*b3*b4+b5*b6*b7-6*b2*b5*b8+b3*b6*b8+b4*b7*b8+b4*b5*b9+b1*b6*b9)))+a2^2*(2*a4*b1*b10*b3-a3*b10*b3^2+4*a3*b1*b10*b4+a10*(3*b1^2*b2+b3^3-3*b1*b3*b4)+6*a8*b1*b2*b5-a9*b3^2*b5+4*a9*b1*b4*b5-3*a8*b3*b4*b5+3*a5*b2*b5^2+a6*b3*b5^2+2*a7*b4*b5^2+2*a9*b1^2*b6-3*a8*b1*b3*b6+4*a7*b1*b5*b6-3*a5*b3*b5*b6+2*a3*b5^2*b6-2*a9*b1*b3*b7+3*a8*b3^2*b7-3*a8*b1*b4*b7+2*a6*b1*b5*b7-2*a7*b3*b5*b7-3*a5*b4*b5*b7+a4*b5^2*b7-3*a5*b1*b6*b7-a7*b1*b7^2+3*a5*b3*b7^2-a3*b5*b7^2+6*a5*b1*b2*b8+2*a6*b1*b3*b8-a7*b3^2*b8+4*a7*b1*b4*b8-3*a5*b3*b4*b8+2*a4*b3*b5*b8+4*a3*b4*b5*b8+4*a3*b1*b6*b8+2*a4*b1*b7*b8-2*a3*b3*b7*b8+(a6*b1^2-2*a7*b1*b3+3*a5*b3^2-3*a5*b1*b4+2*a4*b1*b5-2*a3*b3*b5-2*a3*b1*b7)*b9+a1*(-3*b5*b6*b7+b7^3+6*b2*b5*b8+b1*(6*b10*b2-3*b6*b9)-3*(b10*b3*b4+b3*b6*b8+b4*b7*b8+b4*b5*b9-2*b3*b7*b9)))+b2*(a6^3*b1^2-a6^2*(a7*b1*b3-a5*b3^2+2*a5*b1*b4-6*a4*b1*b5+a3*b3*b5+2*a1*b4*b5+2*a1*b1*b6+a3*b1*b7-2*a1*b3*b7)+2*a4^3*(b1*b10+b5*b8)+a4*(3*a5*a9*b1*b2+3*a1*a3*b10*b2-4*a5*a8*b2*b3+a3^2*b10*b4-a5*a9*b3*b4+2*a5*a8*b4^2+a10*(3*a3*b1*b2-4*a1*b2*b3-a3*b3*b4+2*a1*b4^2)+3*a3*a8*b2*b5+3*a1*a9*b2*b5+2*a3*a9*b4*b5+2*a3*a9*b1*b6-a3*a8*b3*b6-a1*a9*b3*b6+2*a5^2*b4*b6+4*a1*a8*b4*b6+2*a1*a5*b6^2+a7^2*(b4*b5+b1*b6)-2*a5^2*b2*b7-4*a1*a8*b2*b7-a3*a8*b4*b7-a1*a9*b4*b7-a3*a5*b6*b7+3*a3*a5*b2*b8+a3^2*b6*b8-(4*a1*a5*b2+a3*a5*b4+a1*a3*b6)*b9+a7*(3*a8*b1*b2+2*a9*b1*b4-a8*b3*b4+3*a5*b2*b5-a5*b3*b6+2*a3*b5*b6-a5*b4*b7-a1*b6*b7+3*a1*b2*b8+2*a3*b4*b8-a1*b4*b9))-a6*(-3*a1*a9*b1*b2+2*a4*a9*b1*b3+4*a1*a8*b2*b3-2*a4*a8*b3^2-a7^2*b1*b4+4*a4*a8*b1*b4+a1*a9*b3*b4-2*a1*a8*b4^2+a5^2*(2*b2*b3-b4^2)-3*a1*a7*b2*b5+2*a4*a7*b3*b5-3*a4^2*b5^2+a1*a7*b3*b6+4*a1*a4*b5*b6-a1^2*b6^2+2*a4*a7*b1*b7+a1*a7*b4*b7-2*a1*a4*b7^2+a5*(-3*a7*b1*b2+a7*b3*b4-3*a3*b2*b5+4*a4*b4*b5+4*a4*b1*b6+a3*b3*b6-4*a1*b4*b6+4*a1*b2*b7-4*a4*b3*b7+a3*b4*b7)-6*a4^2*b1*b8+4*a1*a4*b4*b8-a3^2*(b5*b6+b4*b8)+2*a1*(a1*b2-2*a4*b3)*b9+a3*(-3*a8*b1*b2-2*a9*b1*b4+a8*b3*b4-2*a7*b4*b5-2*a7*b1*b6+2*a4*b5*b7+a1*b6*b7-3*a1*b2*b8+2*a4*b3*b8+2*a4*b1*b9+a1*b4*b9))-a4^2*(-a10*b3^2+2*a10*b1*b4+2*a1*b10*b4+a9*b3*b5+2*a8*b4*b5+2*a8*b1*b6+2*a5*b5*b6+a9*b1*b7-2*a8*b3*b7+a7*b5*b7-a5*b7^2+a7*b3*b8+2*a5*b4*b8+2*a1*b6*b8+a7*b1*b9-2*a5*b3*b9-2*a1*b7*b9+a3*(b10*b3+b7*b8+b5*b9))+b2*(-a7^3*b1-a3^3*b10+(3*a1^2*a10+a5^3+6*a1*a5*a8)*b2-a1*a9*(2*a5*b4+a1*b6)+a7^2*(a5*b3-3*a3*b5+a1*b7)-a3*(a5*(-2*a9*b3+2*a8*b4)+a5^2*b6+2*a1*(a10*b4+a8*b6-a9*b7))+a3^2*(a10*b3-3*a9*b5+a8*b7+a5*b9)-a7*(a5^2*b4+2*a1*(-a9*b3+a8*b4+a5*b6)+3*a3^2*b8-2*a3*(-3*a9*b1+a8*b3+a5*b7+a1*b9)))));
c7=(-3*a2^3*b1*(b5^2+b1*b8)+a4^3*b2*(b5^2+2*b1*b8)-a4^2*(a9*b1*b2*b3-a8*b2*(b3^2-2*b1*b4)-6*a6*b1*b2*b5+a7*b2*b3*b5+2*a5*b2*b4*b5+a2*b4*b5^2+2*a5*b1*b2*b6+2*a2*b1*b5*b6+2*a1*b2*b5*b6+a7*b1*b2*b7-2*a5*b2*b3*b7+a3*b2*b5*b7-a1*b2*b7^2+a3*b2*b3*b8+2*a2*b1*b4*b8+2*a1*b2*b4*b8+a3*b1*b2*b9-2*a1*b2*b3*b9)+a2^2*(a9*b1*(-b3^2+2*b1*b4)+a8*(3*b1^2*b2+b3^3-3*b1*b3*b4)+6*a5*b1*b2*b5+2*a6*b1*b3*b5-a7*b3^2*b5+4*a7*b1*b4*b5-3*a5*b3*b4*b5+3*a1*b2*b5^2+2*a3*b4*b5^2+2*a7*b1^2*b6-3*a5*b1*b3*b6+4*a3*b1*b5*b6-3*a1*b3*b5*b6+a6*b1^2*b7-2*a7*b1*b3*b7+3*a5*b3^2*b7-3*a5*b1*b4*b7-2*a3*b3*b5*b7-3*a1*b4*b5*b7-3*a1*b1*b6*b7-a3*b1*b7^2+3*a1*b3*b7^2+6*a1*b1*b2*b8-a3*b3^2*b8+4*a3*b1*b4*b8-3*a1*b3*b4*b8-2*a3*b1*b3*b9+3*a1*b3^2*b9-3*a1*b1*b4*b9)-a2*(6*a1*a8*b1*b2^2-2*a7^2*b1*b2*b3-4*a3*a9*b1*b2*b3+2*a3*a8*b2*b3^2+2*a1*a9*b2*b3^2+a6^2*b1^2*b4+a3*a8*b1*b2*b4+a1*a9*b1*b2*b4-6*a1*a8*b2*b3*b4+a7^2*b1*b4^2+2*a3*a9*b1*b4^2-a3*a8*b3*b4^2-a1*a9*b3*b4^2+2*a1*a8*b4^3+a5^2*(3*b1*b2^2-3*b2*b3*b4+b4^3)-4*a3*a7*b2*b3*b5+a1*a7*b2*b4*b5+2*a3*a7*b4^2*b5+a1*a7*b1*b2*b6+4*a3*a7*b1*b4*b6-2*a1*a7*b3*b4*b6+a1*a3*b2*b5*b6+2*a3^2*b4*b5*b6+a3^2*b1*b6^2-a1*a3*b3*b6^2+3*a1^2*b4*b6^2-4*a3*a7*b1*b2*b7+4*a1*a7*b2*b3*b7-a1*a7*b4^2*b7-2*a3^2*b2*b5*b7-3*a1^2*b2*b6*b7-2*a1*a3*b4*b6*b7+2*a1*a3*b2*b7^2-a6*(a7*b1*(-3*b1*b2+b3*b4)+a5*(b1*b2*b3-b3^2*b4+2*b1*b4^2)-6*a3*b1*b2*b5+a1*b2*b3*b5+a3*b3*b4*b5+2*a1*b4^2*b5+a3*b1*b3*b6-a1*b3^2*b6+4*a1*b1*b4*b6+a1*b1*b2*b7+a3*b1*b4*b7-2*a1*b3*b4*b7)+a5*(a7*(2*b2*b3^2+b1*b2*b4-b3*b4^2)+6*a1*(b2^2*b5+b4^2*b6-b2*(b3*b6+b4*b7))+a3*(b2*(b4*b5+b1*b6+4*b3*b7)-b4*(2*b3*b6+b4*b7)))+3*a1^2*b2^2*b8-2*a3^2*b2*b3*b8+a1*a3*b2*b4*b8+a3^2*b4^2*b8-2*a3^2*b1*b2*b9+4*a1*a3*b2*b3*b9-3*a1^2*b2*b4*b9-a1*a3*b4^2*b9)+b2*(a1*(a6^2*(b3^2-2*b1*b4)+b2*(3*a5^2*b2+3*a1*a8*b2+a7^2*b3-2*a5*a7*b4-a1*a9*b4-a1*a7*b6)+a6*(3*a7*b1*b2-4*a5*b2*b3-a7*b3*b4+2*a5*b4^2+2*a1*b4*b6-2*a1*b2*b7))-a3*(3*a7^2*b1*b2+a6^2*b1*b3-2*a1*a9*b2*b3+a5^2*b2*b4+2*a1*a8*b2*b4-3*a1*a6*b2*b5+a1*a6*b3*b6+a5*(-3*a6*b1*b2-2*a7*b2*b3+a6*b3*b4+2*a1*b2*b6)+a1*a6*b4*b7-2*a7*(a6*b1*b4+a1*b2*b7))-a3^3*b2*b8+a3^2*(-3*a9*b1*b2+a8*b2*b3-3*a7*b2*b5+a6*b4*b5+a6*b1*b6+a5*b2*b7+a1*b2*b9))+a4*(3*a6^2*b1^2*b2-2*a6*(a7*b1*b2*b3-a5*b2*(b3^2-2*b1*b4)+a3*b2*b3*b5+2*a2*b1*b4*b5+2*a1*b2*b4*b5+a2*b1^2*b6+2*a1*b1*b2*b6+a3*b1*b2*b7-2*a1*b2*b3*b7)+a2*(a9*b1*(-3*b1*b2+b3*b4)+a8*(b1*b2*b3-b3^2*b4+2*b1*b4^2)-6*a7*b1*b2*b5+a5*b2*b3*b5+a7*b3*b4*b5+2*a5*b4^2*b5-3*a3*b2*b5^2+a7*b1*b3*b6-a5*b3^2*b6+4*a5*b1*b4*b6+a3*b3*b5*b6+4*a1*b4*b5*b6+2*a1*b1*b6^2+a5*b1*b2*b7+a7*b1*b4*b7-2*a5*b3*b4*b7+a1*b2*b5*b7+a3*b4*b5*b7+a3*b1*b6*b7-2*a1*b3*b6*b7-a1*b4*b7^2-6*a3*b1*b2*b8+a1*b2*b3*b8+a3*b3*b4*b8+2*a1*b4^2*b8+a1*b1*b2*b9+a3*b1*b4*b9-2*a1*b3*b4*b9)+a2^2*(b3*(b5^2+2*b1*b8)+b1*(2*b5*b7+b1*b9))+b2*(3*a1*a9*b1*b2-4*a1*a8*b2*b3+a7^2*b1*b4-a1*a9*b3*b4+2*a1*a8*b4^2+a5^2*(-2*b2*b3+b4^2)+3*a1*a7*b2*b5-a1*a7*b3*b6+a1^2*b6^2-a1*a7*b4*b7+a5*(3*a7*b1*b2-a7*b3*b4+3*a3*b2*b5-a3*b3*b6+4*a1*b4*b6-4*a1*b2*b7-a3*b4*b7)+a3^2*(b5*b6+b4*b8)-2*a1^2*b2*b9+a3*(3*a8*b1*b2+2*a9*b1*b4-a8*b3*b4+2*a7*b4*b5+2*a7*b1*b6-a1*b6*b7+3*a1*b2*b8-a1*b4*b9))));
c8=(-3*a2^3*b1^2*b5+a2^2*(a6*b1^2*b3-a7*b1*b3^2+2*a7*b1^2*b4+a5*(3*b1^2*b2+b3^3-3*b1*b3*b4)+6*a1*b1*b2*b5+2*a4*b1*b3*b5-a3*b3^2*b5+4*a3*b1*b4*b5-3*a1*b3*b4*b5+2*a3*b1^2*b6-3*a1*b1*b3*b6+(a4*b1^2-2*a3*b1*b3+3*a1*b3^2-3*a1*b1*b4)*b7)+a2*(-a4^2*b1*(2*b4*b5+b1*b6)+a1*(a6*b1*b2*b3-2*a7*b2*b3^2-a7*b1*b2*b4-a6*b3^2*b4+2*a6*b1*b4^2+a7*b3*b4^2-2*a5*(3*b1*b2^2-3*b2*b3*b4+b4^3)-3*a1*b2^2*b5+3*a1*b2*b3*b6-3*a1*b4^2*b6+3*a1*b2*b4*b7)+a4*(a5*b1*b2*b3-2*a6*b1^2*b4-a5*b3^2*b4+2*a5*b1*b4^2+a7*b1*(-3*b1*b2+b3*b4)+a1*b2*b3*b5+2*a1*b4^2*b5-a1*b3^2*b6+4*a1*b1*b4*b6+a1*b1*b2*b7-2*a1*b3*b4*b7)+a3*(4*a7*b1*b2*b3-2*a5*b2*b3^2-a5*b1*b2*b4-2*a7*b1*b4^2+a5*b3*b4^2+a6*b1*(-3*b1*b2+b3*b4)-6*a4*b1*b2*b5-a1*b2*b4*b5+a4*b3*b4*b5-a1*b1*b2*b6+a4*b1*b3*b6+2*a1*b3*b4*b6-4*a1*b2*b3*b7+a4*b1*b4*b7+a1*b4^2*b7)+a3^2*(-b4*(b4*b5+2*b1*b6)+2*b2*(b3*b5+b1*b7)))+b2*(2*a4^3*b1*b5+a3^2*(-3*a7*b1*b2+a5*b2*b3+a6*b1*b4-a3*b2*b5)+a1^2*(3*a5*b2^2+a6*(-2*b2*b3+b4^2)-b2*(a7*b4+a3*b6))+a1*a3*(3*a6*b1*b2+2*a7*b2*b3-2*a5*b2*b4-a6*b3*b4+a3*b2*b7)+a4^2*(3*a6*b1^2-a7*b1*b3+a5*b3^2-2*a5*b1*b4-a3*b3*b5-2*a1*b4*b5-2*a1*b1*b6-a3*b1*b7+2*a1*b3*b7)+a4*(a3^2*(b4*b5+b1*b6)+a1*(3*a7*b1*b2-4*a5*b2*b3+2*a6*b3^2-4*a6*b1*b4-a7*b3*b4+2*a5*b4^2+2*a1*b4*b6-2*a1*b2*b7)+a3*(3*a5*b1*b2-2*a6*b1*b3+2*a7*b1*b4-a5*b3*b4+3*a1*b2*b5-a1*b3*b6-a1*b4*b7))));
c9=(-a2^3*b1^3+a2*(b2*(-3*b1*(a3*a4*b1+a1^2*b2)+(2*a3^2+a1*a4)*b1*b3-2*a1*a3*b3^2)-(a4^2*b1^2+a1*b2*(a3*b1-3*a1*b3)+a4*b3*(-a3*b1+a1*b3))*b4+(-a3^2*b1+2*a1*a4*b1+a1*a3*b3)*b4^2-a1^2*b4^3)+a2^2*(b1*(a4*b1*b3-a3*b3^2+2*a3*b1*b4)+a1*(3*b1^2*b2+b3^3-3*b1*b3*b4))+b2*(a4^3*b1^2+b2*(-a3^3*b1+a1^3*b2+a1*a3^2*b3-a1^2*a3*b4)+a4^2*(-a3*b1*b3+a1*(b3^2-2*b1*b4))+a4*(a3^2*b1*b4+a1*a3*(3*b1*b2-b3*b4)+a1^2*(-2*b2*b3+b4^2))));

x1List=roots([c9;c8;c7;c6;c5;c4;c3;c2;c1;c0]);

sols=zeros(2,4*9);
costVals=zeros(4*9,1);
numSol=0;
for curSol=1:9
    x1=x1List(curSol);
    %There are two equations that we can try to solve. Solve the one that
    %has the highest weight on the cubic term.
    if(abs(a2)>abs(b2))
        c0=a1*x1^3+a5*x1^2+a8*x1+a10;
        c1=a3*x1^2+a7*x1+a9;
        c2=a4*x1+a6;
        c3=a2;
    else
        c0=b1*x1^3+b5*x1^2+b8*x1+b10;
        c1=b3*x1^2+b7*x1+b9;
        c2=b4*x1+b6;
        c3=b2;
    end
    x2Sol=roots([c3;c2;c1;c0]);
    
    for k=1:3
        curSolVal=[x1;x2Sol(k)];
        costVal=abs(polyValMultiDim(coeffsA,curSolVal))+abs(polyValMultiDim(coeffsB,curSolVal));
        if(costVal<AbsTol)
            numSol=numSol+1;
            sols(:,numSol)=curSolVal;
            costVals(numSol)=costVal;
        end
    end
end

sols=sols(:,1:numSol);

%If sols contains more than 9 solutions (the maximum allowed), then some of
%We will only keep the one that has the lowest cost and add others the
%maximize the minimum distance from the other solutions.
if(numSol>9)
    distMat=Inf*ones(numSol,numSol);

    for curSol1=1:(numSol-1)
       for curSol2=(curSol1+1):numSol
           diff=sols(:,curSol1)-sols(:,curSol2);
           distMat(curSol1,curSol2)=diff'*diff;
           distMat(curSol2,curSol1)=distMat(curSol1,curSol2);
       end
    end

    solsKeptIdx=zeros(9,1);
    solsKeptBool=false(numSol,1);

    [~,minIdx]=min(costVals);
    %Keep the minimum cost solution.
    solsKeptIdx(1)=minIdx;
    solsKeptBool(minIdx)=true;

    %The next solutions added maximize the minimum distance from the
    %solutions already added.
    numSolFound=1;
    for k=2:9
        maxDist=-Inf;
        maxIdx=0;

        for curSol=1:numSol
            if(solsKeptBool(curSol)==true)
                continue;
            end
            %For this solution, find the previously chosen solution to
            %which it is closest. 
            minDist=Inf;
            minIdx=0;
            for prevSol=1:numSolFound
                prevSolIdx=solsKeptIdx(prevSol);
                if(distMat(curSol,prevSolIdx)<minDist)
                    minDist=distMat(curSol,prevSolIdx);
                    minIdx=curSol;
                end
            end

            if(minDist>maxDist)
               maxDist=minDist;
               maxIdx=minIdx;
            end
        end

        solsKeptIdx(k)=maxIdx;
        solsKeptBool(maxIdx)=true;
        numSolFound=numSolFound+1;
    end
    sols=sols(:,solsKeptIdx);
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
