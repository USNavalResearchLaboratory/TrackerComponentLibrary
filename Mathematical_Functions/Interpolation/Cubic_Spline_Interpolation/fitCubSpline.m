function C=fitCubSpline(tau,y,boundValStart,boundValEnd,boundValTypes)
%%FITCUBSPLINE Fit a series of piecewise cubic splines to a set of provided
%       points. Interpolation of the function and its derivatives can then
%       be performed using the evalCubSpline function. Unlike Matlab's
%       spline function, one can specify second derivatives values at the
%       endpoints, if desired.
%
%INPUTS: tau The NX1 or 1XN set of values of the independent variable at
%            which the dependent variable, y, is evaluated.
%          y The XN1 or 1XN dependent variables values evaluated at the
%            points in tau.
% boundValStart, boundValEnd If the first or second derivative at either
%            endpoint is available, then it can be specified here as the
%            "bounds". If no information about the derivative at either end
%            is available, then these can be omitted or an empty matrix can
%            be passed.
% boundValTypes This is a 2X1 vector holding [bound type start; bound type
%            end]. If either or both of boundValStart and boundValEnd are
%            provided, then boundValTypes(1) specifies the bound type in
%            boundValStart and boundValTypes(2) specified the bound type in
%            boundValEnd. If either the start or end bound is not given,
%            then the corresponding entry in boundValTypes is ignored. THe
%            default types if this is omitted or an empty matrix is passed
%            is 1 for both values. The types can be
%            1 The bound contains the first derivative.
%            2 The bound contrains the second derivative.
%
%OUTPUTS: C A 4XN matrix. C(1,:) is the same as y(1:N). C(i-1,1:(N-1))
%           holds the ith derivatives evaluated at tau(1:(N-1)). The
%           values in C(2:4,N) are not used for interpolation.
%
%Piecewise cubic spline fitting with edge constraints is derived in Chapter
%4 of [1]. This function implements the algorithm described as CUBSPL
%therein.
%
%EXAMPLE:
%This is an example of approximating the exponential function. Since the
%derivatives of the function are the same as the fucntion, it is easy to
%see the effect of including them. We plot the function, the interpolation
%of the function with no endpoint info (the "not-a-knot" condition in [1]),
%the function given the true endpoint second derivatives, and also the
%function given zero second derivatives at the endpoints, which is
%sometimes known as a "natural spline".
% x=[-1;0;1;2];
% y=exp(x);
% C=fitCubSpline(x,y);
% CWith2Derivs=fitCubSpline(x,y,y(1),y(end),[2;2]);
% CNatural=fitCubSpline(x,y,0,0,[2;2]);
% 
% figure(1)
% clf
% hold on
% xp=linspace(-2,3,1000);
% plot(xp,exp(xp),'-k','linewidth',6);
% yInterp=evalCubSpline(xp,C,x);
% plot(xp,yInterp,'-b','linewidth',2)
% yWithGrad=evalCubSpline(xp,CWith2Derivs,x);
% plot(xp,yWithGrad,'-y','linewidth',2)
% yNatural=evalCubSpline(xp,CNatural,x);
% plot(xp,yNatural,'--r','linewidth',2)
% scatter(x,y,400,'.g')
% xlabel('x')
% ylabel('y')
% legend('True Function','No End Info','Given Second Derivatives','Zero Second Derivatives','Fitted Points','location','northwest')
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(tau);

if(N<2)
   error('At least two points must be provided.') 
end

C=NaN(4,N);
C(1,:)=y;

if(nargin<5||isempty(boundValTypes))
    boundValTypes=[1;1];
end

if(nargin<3||isempty(boundValStart))
    boundaryCondStart=0;
else
    boundaryCondStart=boundValTypes(1);
    C(2,1)=boundValStart;
end

if(nargin<3||isempty(boundValEnd))
    boundaryCondEnd=0;
else
    boundaryCondEnd=boundValTypes(2);
    C(2,N)=boundValEnd;
end

%As noted in [1], a tridiagonal linear system solving for the unknown
%slopes at the points in tau is generated and solved by Gaussian
%elimination. The slopes are placed in the second row of C. The third and
%forth rows of C are used for temporary storage during this process and
%then in the end, they are updated to hold the interpolating polynomial
%coefficients.

%Put differences between tau values into C(3,:) and put the first divided
%differences of the data into C(4,:).
for m=2:N
    C(3,m)=tau(m)-tau(m-1);
    C(4,m)=(C(1,m)-C(1,m-1))/C(3,m);
end

if(boundaryCondStart==0)%Not a knot condition.
    if(N>2)
        %Not-a-knot condition at the left end and N>2.
        C(4,1)=C(3,3);
        C(3,1)=C(3,2)+C(3,3);
        C(2,1)=((C(3,2)+2*C(3,1))*C(4,2)*C(3,3)+C(3,2)^2*C(4,3))/C(3,1);
    else
        %N==2 and no condition at the left end. 
        C(4,1)=1;
        C(3,1)=1;
        C(2,1)=C(4,2);
    end
elseif(boundaryCondStart==1)
    %Slope at left end is given.
    C(4,1)=1;
    C(3,1)=0;
else
    %Second derivative at the left end is given.
    C(4,1)=2;
    C(3,1)=1;
    C(2,1)=3*C(4,2)-C(3,2)/2*C(2,1);
end

if(N>2)
    %If interior knots are present, perform the forward pass of Gaussian
    %elimination.
    for m=2:(N-1)
        g=-C(3,m+1)/C(4,m-1);
        C(2,m)=g*C(2,m-1)+3*(C(3,m)*C(4,m+1)+C(3,m+1)*C(4,m));
        C(4,m)=g*C(3,m-1)+2*(C(3,m)+C(3,m+1));
    end

    %This section constructs the second boundary condition, if given.
	if(boundaryCondEnd==0)
        if(N==3&&boundaryCondStart==0)
            %N=3 and a not-a-knot condition on the left.
            C(2,N)=2*C(4,N);
            C(4,N)=1;
            g=-1/C(4,N-1);
           
            %Complete forward pass of Gaussian elimination
            C(4,N)=g*C(3,N-1)+C(4,N);
            C(2,N)=(g*C(2,N-1)+C(2,N))/C(4,N);
        else
            %Not-a-knot condition on the left and N>3.
            g=C(3,N-1)+C(3,N);
            C(2,N)=((C(3,N)+2*g)*C(4,N)*C(3,N-1)+C(3,N)^2*(C(1,N-1)-C(1,N-2))/C(3,N-1))/g;
            g=-g/C(4,N-1);
            C(4,N)=C(3,N-1);
           
            %Complete forward pass of Gaussian elimination
            C(4,N)=g*C(3,N-1)+C(4,N);
            C(2,N)=(g*C(2,N-1)+C(2,N))/C(4,N);
        end
    elseif(boundaryCondEnd==2)
        %Second derivative prescribed at right endpoint.
        C(2,N)=3*C(4,N)+(C(3,N)/2)*C(2,N);
        C(4,N)=2;
        g=-1/C(4,N-1);
        
        %Complete forward pass of Gaussian elimination
        C(4,N)=g*C(3,N-1)+C(4,N);
        C(2,N)=(g*C(2,N-1)+C(2,N))/C(4,N);
	end
else%N==2.
    %This section constructs the second boundary condition, if given.
    if(boundaryCondEnd==0)
        if(boundaryCondStart>0)
            %N=2 and a not-a-knot condition on the left.
            C(2,N)=2*C(4,N);
            C(4,N)=1;
            g=-1/C(4,N-1);
            
            %Complete forward pass of Gaussian elimination
            C(4,N)=g*C(3,N-1)+C(4,N);
            C(2,N)=(g*C(2,N-1)+C(2,N))/C(4,N);
        else
            %Not-a-knot condition on the left and right and N=2.
            C(2,N)=C(4,N);
        end
	elseif(boundaryCondEnd==2)
        %Second derivative prescribed at right endpoint.
        C(2,N)=3*C(4,N)+(C(3,N)/2)*C(2,N);
        C(4,N)=2;
        g=-1/C(4,N-1);
        
        %Complete forward pass of Gaussian elimination
        C(4,N)=g*C(3,N-1)+C(4,N);
        C(2,N)=(g*C(2,N-1)+C(2,N))/C(4,N);
    end
end

%Perform back substitution.
for j=(N-1):-1:1
    C(2,j)=(C(2,j)-C(3,j)*C(2,j+1))/C(4,j); 
end

%Generate the cubic coefficients in each interval.
for i=2:N
   dTau=C(3,i);
   divDif1=(C(1,i)-C(1,i-1))/dTau;
   divDif3=C(2,i-1)+C(2,i)-2*divDif1;
   C(3,i-1)=2*(divDif1-C(2,i-1)-divDif3)/dTau;
   C(4,i-1)=(divDif3/dTau)*(6/dTau);
end

%The C matrix is returned.

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
