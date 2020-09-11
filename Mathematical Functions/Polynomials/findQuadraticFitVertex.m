function x=findQuadraticFitVertex(x1,x2,param1,param2,param3)
%%FINDQUADRATICFIRVERTEX Find the minimum of a parabola that fits given
%           scalar data when presented with one of two parameterizations.
%           the first parameterization uses two function values and one
%           derivaitve value. The second parameterization uses just two
%           derivative values.
%
%INPUTS: x1, x2 The two points at which the values are taken.
% param1,param2,param3 These thee values can be given in two ways, which
%               determine their meaning. The first is:
%               f1 The value of the function at x1.
%               f2 The value of the function at x2.
%               g1 The derivative of the function at x1.
%               The other parameterization is
%               g1 The gradient at x1.
%               g2 The gradient at x2.
%               param3 is either omitted or an empty matrix is passed when
%               using the second formulation.
%
%OUTPUTS: x The interpolated critical point (minimum or maximum) of the
%           fitted parabola.
%
%For the first formulation of the problem, we are given equations with
%unknown coefficients (a,b,c):
% f1==a*x1^2+b*x1+c
% f2==a*x2^2+b*x2+c
% gl==b+2*a*x1
%Solving for (a,b,c) one gets
%a=(-f1+f2+g1*(x1-x2))/(x1-x2)^2
%b=(2*f1*x1-x1*(2*f2+g1*x1)+g1*x2^2)/(x1-x2)^2
%c=(f2*x1^2+g1*x1*(x1-x2)*x2+f1*x2*(-2*x1+x2))/(x1-x2)^2
%Substitutign into the expression for the derivative and setting the
%derivative equal to zero, one gets a value of x of 
%x=x1+(g1*(x2-x1)^2)/(2*(f1-f2+g1*(x2-x1)))
%
%For the second formulation of the problem, we are given equations with
%unknown coefficients (a,b,c):
% f1==a*x1^2+b*x1+c
% g1==b+2*a*x1
% g2==b+2*a*x2
%Solving for (a,b,c) one gets
% a=(g1-g2)/(2*x1-2*x2)
% b=(g2*x1-g1*x2)/(x1-x2)
% c=f1+(x1*(-(g1+g2)*x1+2*g1*x2))/(2*(x1-x2))
%Substituting into the expression for the derivative and setting the
%derivative equal to zero, one gets a value of x of
% x=x2+g2*(x1-x2)/(g2-g1);
%which does not require the value f1.
%
%EXAMPLE:
%Here, we use the function f=3*(x-12)*(x-14)=3*x^2-78*x+504;
% x1=-10;
% x2=24;
% f1=3*x1^2-78*x1+504;
% f2=3*x2^2-78*x2+504;
% g1=2*3*x1-78;
% g2=2*3*x2-78;
% xSol1=findQuadraticFitVertex(x1,x2,f1,f2,g1)
% xSol2=findQuadraticFitVertex(x1,x2,g1,g2)
%One will find that both xSol1 and xSol2 =13, which is the minimum point of
%the parabola.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>4&&~isempty(param3))
    %Given the first type of parameterization.
    f1=param1;
    f2=param2;
    g1=param3;

    
    diff=x2-x1;
    x=x1+(g1*diff^2)/(2*(f1-f2+g1*diff));
else
    %Given the second type of parameterization.
    g1=param1;
    g2=param2;
    
    x=x2+g2*(x1-x2)/(g2-g1);
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
