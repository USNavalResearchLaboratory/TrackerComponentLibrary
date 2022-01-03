function [x,t,intersects]=constrainedLSLineSpher(a,b,alphaVal)
%%CONSTRAINEDLSLINESPHER This function solves the real optimization problem
%       arg min _{x,t} norm(a+b*t-x) such that norm(x)=alphaVal
%       where t is a scalar and x is a vector. That is, it finds the
%       closest point on the unit sphere to a line in an arbitrary number
%       of dimensions. If the line does not intersect the sphere, there
%       will be 1 solution. If the line intersects the sphere, there can be
%       1 or two solutions. All solutions and inputs are real.
%
%INPUTS: a A real xDimX1 vector.
%        b A real xDimX1 vector.
% alphaVal The norm of the solution. If this is omittd or ane mpty matrix
%          is passed, the default of 1 is used.
%
%OUTPUTS: x An xDimX1 or xDimX2 set of solutions.
%         t The t value(s) associated with the solution(s) in x. If b=0,
%           then t will be 0, though its value does not matter.
% intersects If the line intersects the ellipsoid this is true. Otherwise,
%           this is false.
%
%We want to solve arg min _{x,t} norm(a+b*t-x)^2 such that
%norm(x)^2=alphaVal. First, we consider two special cases. The first is if
%b is all zero, in which case t does not matter and we are finding the unit
%vector x such that solves min _{x} norm(a-x). The solution is just
%x=sqrt(alphaVal)*(a/norm(a)). Note that if a=all zeros, then a bunch of
%NaNs will be returned for x --every point on the sphere is the "closest".
%
%The second special case is if the line intersects the unit sphere. We do
%not know a priori if the line does intersect the sphere, so we solve the
%problem and if the solutions are complex, then it does not intersect the
%sphere. In this case, we just substitute the line equation into the
%constraint:
%norm(a+b*t)^2=alphaVal
%This is a quadratic equation and can be solved using the quadratic
%formula.
%
%Finally, if the solutions to the above are complex, we solve the general
%problem. In this instance, the Lagrangian for equality constrained
%optimization is
%L=norm(a+b*t-x)^2+lambda*(norm(x)^2-alphaVal)
%where lambda is the Lagrgian parameter. Taking the gradient with respect
%to x and the derivative with respect to t and setting it equal to 0 leads
%to the following linear system of equations:
%[-b',                   norm(b)^2;
%  (1+lambda)*eye(xDim),        -b] *[x;t] == [-a'*b;a]
%The matrix on the left can be explicitly inverted as
%1/(lambda*norm(b)^2)*[b, (1/(1+lambda)*(b*b'+lambda*norm(b)^2*eye(xDim));
%                       1+lambda,  b'];
%This leads to solutions for t and x as
%t=-a'b/norm(b)^2
%x=(1/(1+lambda))*(a-(1/norm(b)^2)*(a'*b)*b)
%Substituting x into the constraint leads to the solutions for lambda:
%lambda=-1+norm(a-(1/norm(b)^2)*(a'*b)*b)
%lambda=-1-norm(a-(1/norm(b)^2)*(a'*b)*b)
%The solution for lambda to us is the one that minimizes norm(a+b*t-x) for
%the given t.
%
%EXAMPLE:
%In this example, the unit circle is plotted along with some lines. The
%points on the circle closest to the lines are plotted in red and black
%circles are used to mark the location on each line that is closest to the
%circle. Dashed lines connect the closest points on the lines to the
%closest points on the circle when the lines do not intersect the circle.
%Generate a unit circle to plot.
% figure(1)
% clf
% hold on
% %Draw a circle.
% drawEllipse([0;0],eye(2,2),1,'linewidth',2)
% axis([-2,2,-2,2])
% axis square
% 
% t=linspace(-10,10,500);
% %Consider 4 lines.
% a=[1,0,0, -2;
%    2,0,-2, -1];
% b=[4,1,-2,3;
%    3,1,-1,-1];
% numLines=size(a,2);
% 
% for k=1:numLines
%     xy=bsxfun(@plus,a(:,k),bsxfun(@times,b(:,k),t));
%     plot(xy(1,:),xy(2,:),'linewidth',2)
%     [x,ts]=constrainedLSLineSpher(a(:,k),b(:,k));
%     
%     scatter(x(1,:),x(2,:),200,'.r')
%     xySol=bsxfun(@plus,a(:,k),bsxfun(@times,b(:,k),ts));
%     scatter(xySol(1,:),xySol(2,:),100,'ok')
%     if(length(ts)==1)
%         %If the line either grazes the unit circle or does not intercept
%         %the unit circle, connect the point on the line to the nearest
%         %point on the circle.
%         plot([x(1),xySol(1)],[x(2),xySol(2)],'--k')
%     end
% end
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(alphaVal))
    alphaVal=1; 
end

%Check for the special case of b being 0. In this instance,
%x=sqrt(alphaVal)*a/norm(a).
if(all(b==0))
    normA=norm(a);
    intersects=normA==1;
    
    x=sqrt(alphaVal)*(a/normA);
    t=0;%t does not matter.
    return; 
end

xDim=size(a,1);

%Deal with the case where the line intersects the unit sphere. This means
%solve for the intersection points and see if they are real. If they are
%complex, then the line does not intersect the unit sphere.
aMag2=dot(a,a);
bMag2=dot(b,b);
ab=dot(a,b);
t=roots([bMag2,2*ab,aMag2-alphaVal]).';
if(~isempty(t)&&isreal(t))
    %The line intersects the sphere of radius alphaVal
    numSol=length(t);
    x=zeros(xDim,numSol);
    
    for curSol=1:numSol
       x(:,curSol)=a+b*t(curSol);
       %Force normalization (deal with any finite precision issues).
       x(:,curSol)=sqrt(alphaVal)*(x(:,curSol)/norm(x(:,curSol)));
    end
    intersects=true;
    return;
end

%In this instance, the line does not intersect the sphere.
intersects=false;
t=-ab/bMag2;

magVal=norm(a-b*(ab/bMag2))/sqrt(alphaVal);
lambda=zeros(2,1);
lambda(1)=-1+magVal;
lambda(2)=-1-magVal;

%For the return values.
x=zeros(xDim,2);
for k=1:2
    x(:,k)=(1/(1+lambda(k)))*a-(1/(bMag2*(1+lambda(k))))*ab*b;
end

linePt=a+b*t;

if(norm(linePt-x(:,1))<norm(linePt-x(:,2)))
    x=x(:,1); 
else
    x=x(:,2);
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
