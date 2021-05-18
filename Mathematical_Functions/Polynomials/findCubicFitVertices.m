function x=findCubicFitVertices(x1,x2,f1,f2,g1,g2)
%%FINDCUBICFITVERTICES Given two real points as well a real scalar function
%        values and the first derivatives of the functions at those points,
%        fit a cubic polynomial and find the local minimum and maximum of
%        the polynomial. This can be useful for interpolation in a step of
%        an algorithm to minimize a function. This function only returns
%        real solutions or an empty matrix in none exist.
%
%INPUTS: x1,x2 The real scalar points at which the function values and its
%              first derivatives are available.
%        f1,f2 Values of the real function at x1 and x2.
%        g1,g2 The first derivative of the real function at x1 and x2. 
%
%OUTPUTS: x A 2X1 vector where x(1) is the local minimum of the cubic
%           function and x(2) is the local maximum of the function. If the
%           function does not have a local minimum or maximum, then this is
%           an empty matrix.
%
%If x1=-1 and x2=1, then the four equations for the function values f and
%the first derivatives g (using unknwon polynomial coefficients (a,b,c,d))
%are:
% fl==-a+b-c+d;
% ft==a+b+c+d;
% gl==3*a-2*b+c;
% gt==3*a+2*b+c;
%Solving for the unknown polynomial coefficients (a,b,c,d), one gets
% a=1/4*(f1-f2+g1+g2);
% b=1/4*(-g1+g2);
% c=1/4*(-3*f1+3*f2-g1-g2);
% d=1/4*(2*f1+2*f2+g1-g2);
%in an equation fo the form f=a*x^3+b*x^2+c*x+d. The maximum and minimum
%points are when the derivative is zero. This is when
%0=3*a*x^2+2*b*x+c
%Solving for x we get the two solutions:
%x=(-b(+/-)sqrt(b^2-3*a*c))/(3*a)
%However, in general x1 ill not be -1 and x2 will not be 1, so we use the
%transformation 
%x=((x2-x1)/2)*y+(x2-x1)/2
%where y goes from -1 to 1 and thus x goes from x1 to x2. This
%transformation means that g1 and g2 must be scaled, and then the final
%value solved must be transformed back.
%
%A second dervative test can easily tell us which is the minimum (the
%minimum is when the second derivative is positive. When the argument of
%the square root is zero or negative, then there are no local minima or
%maxima.
%
%EXAMPLE:
%We use the function f(x)=(x-4)*(x+5)*(x-12)=x^3-11*x^2-32*x+240
% x1=-6;
% x2=1;
% f1=x1^3-11*x1^2-32*x1+240;
% f2=x2^3-11*x2^2-32*x2+240;
% g1=3*x1^2-22*x1-32;
% g2=3*x2^2-22*x2-32;
% x=findCubicFitVertices(x1,x2,f1,f2,g1,g2)
%We find x is about [8.5770; -1.2436]
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

scalFact=(x2-x1)/2;
offset=(x1+x2)/2;
g1=g1*scalFact;
g2=g2*scalFact;

a=f1-f2+g1+g2;
b=g2-g1;
c=3*(f2-f1)-g1-g2;

discriminant=b^2-3*a*c;
if(discriminant<=0)
    %There is no real maximum or minimum.
    x=[];
    return;
end

gamma=sqrt(discriminant);

x(1)=(-b+gamma)/(3*a);
x(2)=(-b-gamma)/(3*a);

%If the first zero is not the minimum, then swap the ordering so that the
%solution at the local function minimum comes first. This is a second
%derivaitve test.
if(3*a*x<-b)
    temp=x(1);
    x(1)=x(2);
    x(2)=temp;
end

x=scalFact*x+offset;

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
