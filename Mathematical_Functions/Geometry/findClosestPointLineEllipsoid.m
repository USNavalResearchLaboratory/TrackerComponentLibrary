function [x,xL,t,intersects]=findClosestPointLineEllipsoid(xc,A,x0,a,aType,epsConst)
%%FINDCLOSESTPOINTLINEELLIPSOID Given a 2D or 3D ellipsoid of the form
%           (x-xc)'*A*(x-xc)=1, and a line of the form x=x0+s*t, where t is
%           a scalar parameteric parameter, find the point on the ellipsoid
%           that is closest to the line as well as the point on the line
%           that is closest to the ellipsoid. If the line intersects the
%           ellipsoid, then more than once point can be returned. All
%           solutions and inputs are real.
%
%INPUTS: xc The xDimX1 center of the ellipsoid.
%         A The xDimXxDim symmetric positive definite matrix defining the
%           above equation for an ellipsoid.
%     x0, a These xDimX1 values define a line and the type of line is given
%           by aType. If aType=0, then the line is parameteric as x=x0+a*t
%           where t is a scalar value. If aType=1, then x0 and a are both
%           points on the line. 
%     aType As mentioned above, this specifies how the line is
%           parameterized. The default if omitted or an empty matrix is
%           passed is 0.
%  epsConst In the rooting formulation of the solution, is is possible that
%           some candidate solutions are invalid (they do not satisfy the
%           constraint (x-xc)'*A*(x-xc)=1. Thus, the absolute value of a
%           transformed version of abs((x-xc)'*A*(x-xc)-1) must be
%           <=epsConst to be valid. The default if omitted or an empty
%           matrix is passed is 1e-9.
%
%OUTPUTS: x The xDimXnumSol set of closest points on the ellipsoid. If the
%           line and ellipsoid do not intersect, numSol=1.
%        xL The xDimXnumSol set of closest points on the line. These
%           correspond to those in x.
%         t The 1XnumSol set of parametric values such that the kth one is
%           xL(:,k)=x0+s*t(k).
% intersects If the line intersects the ellipsoid, this is true. Otherwise,
%           it is false.
%
%This function utilizes a change of coordiantes to solve the problem.
%We have to solve the optimization problem.
%minimize norm(x0+s*t-x)^2
%such that
%(x-xc)'*A*(x-xc)=1
%We perform an eigenvalue decomposition of A as A=U*D*U', where D is a
%diagonal matrix and U is in a unitary matrix (because A is symmetric and
%positive definite). Next, perform the change of variables
%x=U*y+xc
%s=U*c1        or c1=U'*s
%c=x0-xc=U*c2  or c2=U'*(x0-xc)
%c=U*c2
%Because U'*U=I, the optimization problem becomes:
%norm(c1*t-y+c2)^2
%Subject to 
%y'*D*y=1
%Thus, we solve for y and transform back into x. The method of solution is
%Lagrangian relaxation. Chapter 4.1 or [1] goes over Lagrangian
%multipliers. The Lagrangian for the problem is
%L=norm(c1*t-y+c2)^2+lambda*(y'*D*y-1)
%Setting the derivatives of L with respect to t and y equal to 0, one gets
%a solution in terms of lambda (the Lagrangian parameter). One can then
%choose lambda to satisfy the constraint y'*D*y=1. In this instance, that
%can be formulated as involves solving a polynomial in terms of lambda.
%That is what this function does. Due to possible singularities in a ratio
%that is reduced to the polynomial, it is possible that invalid solutiosn
%will be found. Thus, the epsConst parameter is introduced to throw out
%possible solutions that do not sufficiently satisfy the condition. Of the
%multiple solutions found for lambda, only the one minimizing the cost
%function is kept.
%
%EXAMPLE 1:
%This example is done in 2D, so it can be easily plotted. Multiple lines
%are drawn. Two intersect an ellipse and two do not intersect the ellipse.
%The closest points on the ellipse are drawn in red. Those on the lines at
%black circles. Dashed black lines connect them.
% xc=[2;2];
% A=[0.5,-0.1;
%   -0.1, 0.2];
% figure(1)
% clf
% hold on
% drawEllipse(xc,A,1,'linewidth',2)
% axis([-2, 6, -2, 6])
% axis square
% t=linspace(-10,10,500);
% %Consider 4 lines.
% a=[1,0,0, -2;
%    2,0,-2, -1];
% b=[4,1,-2,3;
%    3,1,-1,-1;];
% numLines=size(a,2);
% for k=1:numLines
%     xy=bsxfun(@plus,a(:,k),bsxfun(@times,b(:,k),t));
%     plot(xy(1,:),xy(2,:),'linewidth',2)
%     [x,xL]=findClosestPointLineEllipsoid(xc,A,a(:,k),b(:,k));
%     
%     scatter(x(1,:),x(2,:),200,'.r')
%     scatter(xL(1,:),xL(2,:),100,'ok')
%     if(size(xL,2)==1)
%         %If the line either grazes the unit circle or does not intercept 
%         %the unit circle, connect the point on the line to the nearest
%         %point on the circle.
%         plot([x(1),xL(1)],[x(2),xL(2)],'--k')
%     end
% end
%
%EXAMPLE 2:
%This example is similar to the first one, but is in 3D.
% xc=[2;2;0];
% A=[0.5,-0.1, 0.15;
%   -0.1, 0.2, 0.15;
%    0.15,0.15,0.3];
% figure(1)
% clf
% hold on
% drawEllipse(xc,A,1,'edgecolor','none')
% light()
% %axis([-2, 6, -2, 6])
% %axis square
% t=linspace(-10,10,500);
% %Consider 4 lines.
% a=[1,0,0, -2;
%    2,0,-2, -1;
%    1,0,1, -1];
% b=[4,1,-2,3;
%    3,1,-1,-1;
%    1,2,0,0];
% numLines=size(a,2);
% for k=1:numLines
%     xyz=bsxfun(@plus,a(:,k),bsxfun(@times,b(:,k),t));
%     plot3(xyz(1,:),xyz(2,:),xyz(3,:),'linewidth',2)
%     [x,xL]=findClosestPointLineEllipsoid(xc,A,a(:,k),b(:,k));
%     scatter3(x(1,:),x(2,:),x(3,:),200,'.r')
%     scatter3(xL(1,:),xL(2,:),xL(3,:),100,'ok')
%     if(size(xL,2)==1)
%         %If the line either grazes the unit circle or does not intercept 
%         %the unit circle, connect the point on the line to the nearest
%         %point on the circle.
%         plot3([x(1),xL(1)],[x(2),xL(2)],[x(3),xL(3)],'--k')
%     end
% end
% axis([-6, 6, -6, 6,-6,6])
% axis square
% view(-8,30)
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 2nd ed. Belmont, MA:
%    Athena Science, 1999.
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(epsConst))
    epsConst=1e-9;
end

if(nargin<5||isempty(aType))
    aType=0;
end

if(aType==1)
    s=a-x0;
else
    s=a;
end

%The special case of s being 0. This means that we just find the closest
%point on the ellipsoid to the point x0.
if(all(s==0))
    x=nearestPointOnEllipsoid(xc,A,x0,1);
    xL=x0;
    t=0;
    intersects=all(x==xL);
    return; 
end

%Check the special case of the line intercepting the ellipsoid. In this
%instance,
[x,t]=findEllipsLineIntersect(xc,A,x0,s,0);
if(~isempty(x))
    t=t';
    xL=x;
    intersects=true;
    return; 
end

%There is a line and it does not intercept the ellipsoid. Find the closest
%point.
intersects=false;

[~,D,U]=eig(A);
c1=U'*s;
c2=U'*(x0-xc);

[y,t]=solveDiagA(D,c2,c1,epsConst);

%Transform back to x
x=bsxfun(@plus,U*y,xc);
xL=bsxfun(@plus,x0,bsxfun(@times,s,t));

end

function [x,t]=solveDiagA(A,a,b,epsConst)

numDim=size(a,1);
switch(numDim)
    case 2
        a11=A(1,1);
        a22=A(2,2);
        a1=a(1);
        a2=a(2);
        b1=b(1);
        b2=b(2);
        
        firstTerm=-a11*a22*(b1^2+b2^2)*(a11*b1^2+a22*b2^2);
        sqrtTerm=sqrt(a11^3*a22^3*(a2*b1-a1*b2)^2*(b1^2+b2^2)^2*(a11*b1^2+a22*b2^2));
        denom=a11^2*a22^2*(b1^2+b2^2)^2;
        
        if(imag(sqrtTerm)~=0)
            %If there is some finite precision error.
            x=[];
            t=[];
            return 
        end

        lambda(1)=(firstTerm+sqrtTerm)/denom;
        lambda(2)=(firstTerm-sqrtTerm)/denom;

        x=zeros(2,2);
        t=zeros(1,2);
        cost=zeros(2,1);
        constErr=zeros(2,1);
        for k=1:2
            denom=a22*b2^2+a11*(b1^2+a22*(b1^2+b2^2)*lambda(k));
            
            x(:,k)=[a22*b2*(a1*b2-a2*b1);a11*b1*(a2*b1-a1*b2)]/denom;
            t(k)=-(a2*a22*b2*(1+a11*lambda(k))+a1*a11*b1*(1+a22*lambda(k)))/denom;
             
            cost(k)=norm(a+b*t(k)-x);
            constErr(k)=abs(x(:,k)'*A*x(:,k)-1);
        end
        sel=(constErr<=epsConst);
        if(sum(sel)==0)
            %If nothing met the constraint.
            x=[];
            t=[];
            return 
        end
        x=x(:,sel);
        t=t(sel);
        cost=cost(sel);
        
        [~,idx]=min(cost);
        
        x=x(:,idx);
        t=t(idx);
    case 3
        a11=A(1,1);
        a22=A(2,2);
        a33=A(3,3);
        
        a1=a(1);
        a2=a(2);
        a3=a(3);
        b1=b(1);
        b2=b(2);
        b3=b(3);

        c0=(a11*b1^2+a22*b2^2+a33*b3^2)*(a22*(-1+a3^2*a33)*b2^2-2*a2*a22*a3*a33*b2*b3+(-1+a2^2*a22)*a33*b3^2+a11*((-1+a2^2*a22+a3^2*a33)*b1^2-2*a1*b1*(a2*a22*b2+a3*a33*b3)+a1^2*(a22*b2^2+a33*b3^2)));
        c1=2*(a11*b1^2+a22*b2^2+a33*b3^2)*(-a11*a33*(b1^2+b3^2)-a22*a33*(b2^2+b3^2)+a11*a22*((-1+(a2^2+a3^2)*a33)*b1^2+(-1+(a1^2+a3^2)*a33)*b2^2-2*a2*a3*a33*b2*b3+(a1^2+a2^2)*a33*b3^2-2*a1*a33*b1*(a2*b2+a3*b3)));
        c2=(-a22^2*a33^2*(b2^2+b3^2)^2+a11^2*(-a33^2*(b1^2+b3^2)^2+a22^2*((-1+a3^2*a33)*(b1^2+b2^2)^2-2*a3*a33*(a1*b1+a2*b2)*(b1^2+b2^2)*b3+a33*(a1*b1+a2*b2)^2*b3^2)+a22*a33*((-4+a2^2*a33)*b1^4-2*a1*a2*a33*b1^3*b2+2*a1*a33*b1*b2*b3*(a3*b2-a2*b3)+b3^2*((-2+a3^2*a33)*b2^2-2*a2*a3*a33*b2*b3+a2^2*a33*b3^2)+b1^2*((-4+a1^2*a33)*b2^2-2*a2*a3*a33*b2*b3+2*(-2+a2^2*a33)*b3^2)))+a11*a22*a33*(-2*a33*(b1^2*b2^2+2*(b1^2+b2^2)*b3^2+2*b3^4)+a22*(-2*a1*a33*b1*(a2*b2+a3*b3)*(b2^2+b3^2)+(b2^2+b3^2)*((-4+a1^2*a33)*b2^2+a1^2*a33*b3^2)+b1^2*((-4+a2^2*a33)*b2^2+2*a2*a3*a33*b2*b3+(-2+a3^2*a33)*b3^2))));
        c3=-2*a11*a22*a33*(b1^2+b2^2+b3^2)*(a11*a22*(b1^2+b2^2)+a11*a33*(b1^2+b3^2)+a22*a33*(b2^2+b3^2));
        c4=-a11^2*(a22^2)*(a33^2)*((b1^2+b2^2+b3^2)^2);

        lambda=roots([c4;c3;c2;c1;c0]);
        lambda=lambda(imag(lambda)==0);
        numSol=length(lambda);
        
        if(numSol==0)
            %If there is some finite precision error. There should be two
            %solutions.
            x=[];
            t=[];
            return 
        end
        
        x=zeros(3,numSol);
        t=zeros(1,numSol);
        cost=zeros(numSol,1);
        %Due to singularities, some solutions might not be correct. Thus,
        %we also evaluate the relative error for enforcing the constraint
        %and if the error is too high, those solutions are discarded.
        constErr=zeros(numSol,1);
        for k=1:numSol
            lambdaCur=lambda(k);
            denom=a22*b2^2+a33*b3^2+a22*a33*(b2^2+b3^2)*lambdaCur+a11*b1^2*(1+a22*lambdaCur)*(1+a33*lambdaCur)+a11*lambdaCur*(a33*b3^2+a22*(b2^2+a33*(b2^2+b3^2)*lambdaCur));
            x(:,k)=[-a3*a33*b1*b3*(1+a22*lambdaCur)-a2*a22*b1*b2*(1+a33*lambdaCur)+a1*(a22*b2^2+a33*b3^2+a22*a33*(b2^2+b3^2)*lambdaCur);
                     a11*b1*(a2*b1-a1*b2)+a33*b3*(-a3*b2+a2*b3)+a11*a33*(-b2*(a1*b1+a3*b3)+a2*(b1^2+b3^2))*lambdaCur;
                     a11*b1*(a3*b1-a1*b3)+a22*b2*(a3*b2-a2*b3)+a11*a22*(a3*(b1^2+b2^2)-(a1*b1+a2*b2)*b3)*lambdaCur]/denom;
            t(k)=-(a1*a11*b1*(1+a22*lambdaCur)*(1+a33*lambdaCur)+(1+a11*lambdaCur)*(a3*a33*b3*(1+a22*lambdaCur)+a2*a22*b2*(1+a33*lambdaCur)))/denom;            
            cost(k)=norm(a+b*t(k)-x);
            constErr(k)=abs(x(:,k)'*A*x(:,k)-1);
        end
        sel=(constErr<=epsConst);
        if(sum(sel)==0)
            %If nothing met the constraint.
            x=[];
            t=[];
            return 
        end
        x=x(:,sel);
        t=t(sel);
        cost=cost(sel);

        [~,idx]=min(cost);
        x=x(:,idx);
        t=t(idx);
    otherwise
        error('An unsupported dimensionality is being used.')
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
