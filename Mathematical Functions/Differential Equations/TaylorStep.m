function xNew=TaylorStep(deltaT,x,order,a,papt,papx,p2apxpt,p2apt2,p2apxpx)
%TAYLORSTEP Perform a single step of a Taylor series method to integrate a
%           deterministic differential equation of the form dxdt=a(x,t).
%
%INPUTS: deltaT The time increment over which the integration is performed.
%             x The xDimX1 initial vector value of the argument being
%               integrated.
%         order The order of the Taylor series expansion to use. This can
%               go from 1 to 3. Higher orders require more inputs. All
%               orders require a. Order 2 also requires papt and papx.
%               Order 3 requires the same as order 2 plus
%               p2apxpt, p2apt2, p2apxpx.
%             a The xDimX1 vector of values of the derivative of x with
%               respect to t.
%          papt The xDimX1 vector of values of the partial derivative of a
%               with respect to t.
%          papx The xDimXxDim matrix of partial derivatives if a with
%               respect to the elements of x. papx(:,i) is the derivative
%               with respect to x(i).
%       p2apxpt The xDimXxDim matrix of partial derivatives of papt with
%               respect to the elements of x. p2apxpt(:,i) is the
%               derivative with respect to x(i).
%        p2apt2 The xDimX1 second partial derivative of a with respect to
%               t.
%       p2apxpx The xDimXxDimXxDim matrix of second derivatives of a with
%               respect to the elements of x. p2apxpx(:,i,j) is the partial
%               derivative of a with respect to x(i) and x(j).
%
%A third-order taylor series expansion (covered in most basic calculus
%books) is:
%xNew=x+a*deltaT+(1/2)*dadt*deltaT^2+(1/6)*d2adt2*deltaT^3
%where dadt and d2adt2 are the first and second order total derivatives of
%a. A first order solution is just x+a*deltaT. A second order expansion
%requires the computation of dadt. From the derfinition of the total
%derivative (covered in basic multivariate calculus books), this is
%dadt=papt+papx*a, where the 'p's denote partial derivative operators,
%since a is the derivative of x with respect to t. Note that papx is a
%matrix. The computation of d2adt2 is more complicated.
%d2adt2=d/dt{dadt}=d/dt{papt}+d/dt{papx*a}.
%The first term is:
%d/dt{papt}=p2apt2+sum_i( p2aptpx_i a(i)+papx_i papt(i) )
%The second term can be expanded using the chain rule:
%d/dt{papx*a}=d/dt{papx}*a+papx*dadt
%This just leaves the term d/dt{papx} to evaluate.
%d/dt{papx}=p2apxpt+sum_i(p2apxpx(i)*a(i))
%
%EXAMPLE 1:
%This is an example, where we know the solution. We start with the
%solution:
%x(1)=sin(omega*t^2)^2, x(2)=cos(omega*t^2)^2, x(3)=cos(omega*t^2)
%We then get second order differential equations in terms of these
%variables with x(4), x(5), and x(6) defined as the derivatives of
%x(1), x(2), and x(3), respectively. The full system is epxressed without
%trigonometric functions. We demonstrate that as the order of the Taylor
%step increases, the solution reaches an increasingly high number of digits
%of accuracy.
% numSteps=100;
% omega=2*pi;
% t0=2.04;
% t1=2.08;
% a=@(x,t)[x(4);
%          x(5);
%          x(6);
%          8*t^2*omega^2*x(3)^2+x(4)/t-2*x(6)^2;
%         -8*t^2*omega^2*x(3)^2+x(5)/t+2*x(6)^2;
%         -4*t^2*omega^2*x(3)+x(6)/t];
% papt=@(x,t)[0;
%             0;
%             0;
%             -(x(4)/t^2)+16*t*x(3)^2*omega^2;
%             -(x(5)/t^2)-16*t*x(3)^2*omega^2;
%             -(x(6)/t^2)-8*t*x(3)*omega^2];
% p2apt2=@(x,t)[0;
%               0;
%               0;
%               (2*x(4))/t^3+16*x(3)^2*omega^2;
%               (2*x(5))/t^3-16*x(3)^2*omega^2;
%               (2*x(6))/t^3-8*x(3)*omega^2];
% p2aptpx=@(x,t)[0,0,0,                   0,       0,        0;
%                0,0,0,                   0,       0,        0;
%                0,0,0,                   0,       0,        0;
%                0,0,32*t*x(3)*omega^2,  -(1/t^2), 0,        0;
%                0,0,-32*t*x(3)*omega^2,  0,       -(1/t^2), 0;
%                0,0,-8*t*omega^2,        0,       0,        -(1/t^2)]; 
% 
% papx=@(x,t)[0,0,0,                      1,      0,      0;
%             0,0,0,                      0,      1,      0;
%             0,0,0,                      0,      0,      1;
%             0,0,16*t^2*x(3)*omega^2,    1/t,    0,      -4*x(6);
%             0,0,-16*t^2*x(3)*omega^2,   0,      1/t,    4*x(6);
%             0,0,-4*t^2*omega^2,         0,      0,      1/t];
% tVals=linspaceNoEnd(t0,t1,numSteps);
% deltaT=tVals(2)-tVals(1);
% exactSol=@(t)[sin(omega*t^2)^2;
%        cos(omega*t^2)^2;
%        cos(omega*t^2);
%        4*omega*t*cos(omega*t^2)*sin(omega*t^2);
%        -4*omega*t*cos(omega*t^2)*sin(omega*t^2);
%        -2*omega*t*sin(omega*t^2)];
% xInit=exactSol(t0);
% x=repmat(xInit,[1,3]);
% for curStep=1:numSteps
%     %Get the function values at the point.
%     t=tVals(curStep);
%     p2apxpxCur=zeros(6,6,6);
%     %Second derivatives with respect to x(1), x(2), x(4), and x(5) are
%     %all zero.
%     %Second derivatives with respect to x(3).
%     p2apxpxCur(:,:,3)=[0,0,0,               0,0,0;
%                        0,0,0,               0,0,0;
%                        0,0,0,               0,0,0;
%                        0,0,16*t^2*omega^2,  0,0,0;
%                        0,0,-16*t^2*omega^2, 0,0,0;
%                        0,0,0,               0,0,0];
%     %Second derivatives with respect to x(6).
%     p2apxpxCur(:,:,6)=[0,0,0,0,0,0;
%                        0,0,0,0,0,0;
%                        0,0,0,0,0,0;
%                        0,0,0,0,0,-4;
%                        0,0,0,0,0,4;
%                        0,0,0,0,0,0];
%     for order=1:3
%         xCur=x(:,order);
%         aCur=a(xCur,t);
%         paptCur=papt(xCur,t);
%         papxCur=papx(xCur,t);
%         p2aptpxCur=p2aptpx(xCur,t);
%         p2apt2Cur=p2apt2(xCur,t);
%         x(:,order)=TaylorStep(deltaT,x(:,order),order,aCur,paptCur,papxCur,p2aptpxCur,p2apt2Cur,p2apxpxCur);
%     end
% end
% solEnd=exactSol(t1);
% solErr=bsxfun(@minus,x,solEnd);
% errNorm=sqrt(sum(solErr.*solErr,1))
%One will see that as the order goes from 1 to 3, the norm of the error
%across the whole estiamted vector goes from about 3e-1 to 4e-3 to 3e-6.
%
%EXAMPLE 2:
%In this example, we demonstrate that one will get an exact solution in a
%single step when the differential equation a is just in terms of t, t is a
%second order polynomial, and one uses the third order approximation.
% order=3;
% t0=0;
% t1=1/2;
% a=@(x,t)(t^2);
% papt=@(x,t)(2*t);
% p2apt2=2;
% papx=0;
% p2apxpt=0;
% p2apxpx=0;
% 
% deltaT=t1-t0;
% xInit=-1;
% aCur=a(xInit,t0);
% paptCur=papt(xInit,t0);
% %Take a single step.
% x=TaylorStep(deltaT,xInit,order,aCur,paptCur,papx,p2apxpt,p2apt2,p2apxpx);
% exactSol=-23/24;
% RelErr=(x-exactSol)/exactSol
%In this instance, the relative error is zero. 
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(order))
    order=1;
end

n=length(x);

switch(order)
    case 1
        xNew=x+a*deltaT;
    case 2
        dadt=papt+papx*a;
        
        xNew=x+deltaT*(a+(1/2)*dadt*deltaT);
    case 3
        dadt=papt+papx*a;
        
        dpaptdt=p2apt2;
        sumVal=zeros(n,1);
        for i=1:n
            sumVal=sumVal+p2apxpt(:,i)*a(i)+papx(:,i)*papt(i);
        end
        dpaptdt=dpaptdt+sumVal;
        
        dpapxdt=p2apxpt;
        for i=1:n
            dpapxdt=dpapxdt+p2apxpx(:,:,i)*a(i);
        end

        d2adt2=dpaptdt+dpapxdt*a+papx*dadt;
        
        xNew=x+deltaT*(a+deltaT*((1/2)*dadt+(1/6)*d2adt2*deltaT));
    otherwise
        error('An unsupported order was requested.')
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
