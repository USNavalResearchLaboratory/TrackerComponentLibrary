function [x,exitCode]=LSEstLMarquardt(fun,x0,hasJacobian,TolG,TolX,delta,deltaAbs,maxIter,maxTries)
%%LSESTLMARQUARDT Perform minimization of function over x having the
%           quadratic form F(x)=(1/2)f(x)'*f(x), where x is a vector and
%           f(x) is a vector function, using the Levenberg-Marquardt
%           algorithm. If no Jacobian for f is available, then a secant
%           version of the algorithm is used to estimate the Jacobian. This
%           algorithm is useful for solving least squares problems. The
%           algorithm typically has better convergence without the need for
%           second derivatives than methods made for more general
%           (non-least squares) optimization problems.
%
%INPUTS: fun A handle to the function f related to the cost as 
%            F(x)=(1/2)f(x)'*f(x). If a Jacobian is available, fun returns
%            f(x) and the associated Jacobian J(x). If no Jacobian is
%            available, then fun just returns f(x). The Jacobian matrix
%            holds in row i and column j the partial derivative of the ith
%            component of f with respect to the jth component of x.
%         x0 The initial estimate of the optimal value of x. This is an NX1
%            vector.
% hasJacobian An optional boolean value indicating whether the function fun
%            returns a Jacobian. If omitted or an empty matrix is passed, a
%            default value of false is used indicating that no Jacobian is
%            present.
%       TolG The absolute convergence tolerance for the gradient of F.
%            Convergence is declared if norm(g,Inf)<=TolG. If omitted or an
%            empty matrix is passed, the default value of 1e-6 is used.
%       TolX The relative convergence tolerance for the magnitude of the
%            descent direction vector h. Convergence is declared if
%            norm(h)<=TolX*(norm(x)+TolX) where h is nominally the negative
%            gradient of F. If omitted or an empty matrix is passed, the
%            default value of 1e-9 is used.
%      delta This parameter is only used if hasJacobian=false. This is the
%            relative stepsize used for finite differencing when using the
%            secant method. The stepsize is chosen such that
%            epsilon=max(delta*abs(x),deltaAbs); The deltaAbs limit is for
%            when elements of x are near zero. If omitted or an empty
%            matrix is passed, the default value of 1e-7 is used.
%   deltaAbs An absolute lower bound for finite differencing, which is only
%            used when the Jacobian is not provided. If omitted or an empty
%            matrix is passed, the default value is delta^2.
%    maxIter The maximum number of iterations of the algorithm to perform
%            before declaring a failure to converge. If omitted or an empty
%            matrix is passed, the default value of 100+10*length(x0) is
%            used.
%   maxTries The Levenberg-Marquart algorithm tries to adjust the damping
%            parameter each iteration before declaring a failure to obtain
%            a damping parameter sufficient to get a positive definite
%            matrix for determining a descent direction. If omitted or an
%            empty matrix is passed, the default value of 100 is used.
%
%OUTPUTS: x The optimal value of x as obtained by the algorithm, or an
%           empty matrix if an error occurred.
%  exitCode A parameter indicating how the algorithm terminated. Possible
%           values are:
%          -2 A nonfinite number was encountered.
%          -1 Solving for the descent direction failed. This can be due to
%             maxTries being reached.
%           0 The maximum number of iterations maxIter was reached.
%           1 Convergence according to TolG was determined.
%           2 Convergence according to TolX was determined.
%           3 The step size became so small that finite precision
%             limitations preclude changes in x.
%
%The Levenberg-Marquart algorithm is a dampened form of the Gauss-Newton
%method for least squares optimization. Unlike the Gauss-Newton methods, it
%does not use a line search to determine the stepsize. The original
%Levenberg-Marquardt algorithm is described in [1]. The implementation here
%follows that of Algorithm 3.16 in [2]. When the Jacobian is not provided,
%the secant version of the algorithm, partially summarized as Algorithm
%3.34 in [2], is used. Algorithm 3.34 in 2 omits a number of steps.
%However, the missing steps coincide with those of Algorithm 3.31, so it is
%not difficult to add the Jacobian update steps of the secant method to
%Algorithm 3.31. In the comments to the code, all equation numbers refer to
%those in [2].
%
%The algorithm can be demonstrated using the following scenario:
% %The true parameters that are to be estimated.
% xTrue=[20;10;1;50];
% 
% t=(1:1:100)';
% %The measurement function.
% fMeas=@(x)(x(1)*exp(-t/x(2))+x(3)*t.*exp(-t/x(4)));
% 
% %Generate noisy measurements. In this instance, they are not noisy, so
% %that one can see that the algorithm converges to the true value.
% z=fMeas(xTrue);
% 
% %The error function. The total cost is (1/2)*f(x)'*f(x).
% f=@(x)(z-fMeas(x));
% 
% %The Jacobian matrix of the error function.
% J=@(x)[-exp(-t/x(2)),-(t*x(1)/x(2)^2).*exp(-t/x(2)),-t.*exp(-t/x(4)),-(x(3)/x(4)^2)*t.^2.*exp(-t/x(4))];
% %Thus, the total function to give to this routine is
% fun=@(x)deal(f(x),J(x));
% %The deal function puts f(x) and J(x) into two outputs.
% 
% %Use an initial estimate that is off by a bit.
% x0=[5;2;0.2;10];
% %Run the algorithm to estimate x.
% [x,exitCode]=LSEstLMarquardt(fun,x0,false)
% %One should get back the true results.
% %On the other hand, if one did not wish to compute the Jacobian, the
% %numerical differentiation approach will also work:
% [x,exitCode]=LSEstLMarquardt(f,x0,true)
% %In both instances, the exitCode should be 1.
%
%REFERENCES:
%[1] D. W. Marquardt, "An algorithm for least-squares estimation of non-
%    linear parameters," Journal of the Society for Industrial and Applied
%    Mathematics, vol. 11, no. 2, pp. 431-441, Jun. 1963.
%[2] K. Madsen, H. B. Nielsen, and O. Tingleff, "Methods for non-linear
%    least squares problems," Informatics and Mathematical Modelling,
%    Technical University of Denmark, Tech. Rep., Apr. 2004. [Online].
%    Available: http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(hasJacobian))
    hasJacobian=false;
end

if(nargin<4||isempty(TolG))
    TolG=1e-6;
end

if(nargin<5||isempty(TolX))
   TolX=1e-9; 
end

if(nargin<6||isempty(delta))
    delta=1e-7;
end

if(nargin<7||isempty(deltaAbs))
   deltaAbs=delta^2; 
end

N=length(x0);

if(nargin<8||isempty(maxIter))
	maxIter=100+10*N;
end

if(nargin<9||isempty(maxTries))
    maxTries=100;
end

%Default value of tau as suggested on page 25 for a mediocre initial
%estimate.
tau=1e-3;

%Initial estimate, as supplied by the user.
x=x0;

%Initial function value f and the Jacobian of the function value J.
if(hasJacobian==false)
    f=fun(x);
    
    %Get the Jacobian via numerical differentiation. Rather than using the
    %forward difference of Eq. 3.28, a central difference is used via the
    %numDiff function, because central differences are more numerically
    %stable. The same criterion for choosing the step size epsilon is used
    %as in the forward difference routine for updating J later.
    epsilon=max(delta*abs(x),deltaAbs);
    J=numDiff(x,fun,length(f),1,epsilon);
else
    [f,J]=fun(x);
end

%If numerical problems arose.
if(any(~isfinite(f))||any(~isfinite(J(:))))
    exitCode=-2;
    return;
end

A=J'*J;%The approximate Hessian.

%The initial value of mu as given by Equation 3.14.
mu=tau*max(diag(A));

%The value of the gradient is from Eq. 3.4a.
g=J'*f;
nu=2;

j=0;%Used by the secant method if no Jacobian is provided.

%If numerical problems arose.
if(any(~isfinite(A(:)))||any(~isfinite(g)))
    exitCode=-2;
    x=[];
    return;
end

%If the starting estimate meets the convergence criterion.
if(norm(g,Inf)<=TolG)
    exitCode=1;
    return;
end

for curIter=1:maxIter
    %Solve (A+mu*I)*h=-g to find h and modify mu, if necessary to avoid
    %numerical problems.
    [h,mu]=solveForH(A,g,mu,maxTries);
    %If solving for h failed.
    if(isempty(h))
        exitCode=-1;
        x=[];
        return;
    end

    %Check for convergence
    if(norm(h)<=TolX*(norm(x)+TolX))
       exitCode=2;
       return;
    end

    %Update the Jacobian if the Secant method is used.
    if(hasJacobian==false)
        j=mod(j,N)+1;
        if(abs(h(j))<0.8*norm(h))
            %Update the Jacobian estimate J by 3.32 using xNew=x+eta*ej
            %with eta as described by the comments below Algorithm 3.34 on
            %page 44. This is essentially using a forward difference to
            %update J. The use of a forward difference over a central
            %difference avoids having to perform one additional function
            %evaluation.

            %First, find xNew.
            if(x(j)==0)
                eta=deltaAbs;
            else
                eta=delta*abs(x(j));
            end
            xNew=x;
            xNew(j)=xNew(j)+eta;

            %Next, perform the update
            fNew=fun(xNew);

            %If numerical problems arose.
            if(any(~isfinite(fNew)))
                exitCode=-2;
                return;
            end

            %This is from Definition 3.32 in page 42: Broyden's rank one
            %update.
            hNew=xNew-x;
            u=(1/(hNew'*hNew))*(fNew-f-J*hNew);
            J=J+u*hNew';%Update B.
        end
    end
    
    %Compute the next estimate of x.
    xNew=x+h;

    if(xNew==x)
        exitCode=3;
        return;
    end

    %If numerical problems arose.
    if(any(~isfinite(xNew)))
       exitCode=-2;
       x=[];
       return;
    end
    
    if(hasJacobian==false)%When using the secant method.
        fNew=fun(xNew);
        
        %Perform the update of J as per Eq. 3.32.
        u=(1/(h'*h))*(fNew-f-J*h);
        JNew=J+u*h';%Update J.
    else
        [fNew,JNew]=fun(xNew);
    end
    
    %If numerical problems arose.
    if(any(~isfinite(fNew))||any(~isfinite(JNew(:))))
        exitCode=-2;
        x=[];
        return;
    end
    
    %Compute the gain ratio from the equation after Equation 3.14.
    %First, we have the numerator computed in a manner that should be
    %more numerically robust than evaluating f'*f-fNew'*fNew.
    deltaF=(f-fNew)'*(f+fNew)/2;
    
    %The denominator
    deltaL=(1/2)*h'*(mu*h-g);
    rho=deltaF/deltaL;%The gain ratio

    %The check for finiteness of rho deals with deltaL being zero, in
    %which case the gain ratio cannot be used. deltaL should theoretically
    %always be positive or zero, so we can just check deltaF as being
    %positive instead of checking rho.
    if(isfinite(rho)&&(deltaF>0)&&(deltaL>0))
        x=xNew;
        J=JNew;
        f=fNew;

        A=J'*J;
        g=J'*f;
        
        %If numerical problems arose.
        if(any(~isfinite(A(:)))||any(~isfinite(g(:))))
            exitCode=-2;
            x=[];
            return;
        end      
        
        %Check for convergence
        if(norm(g,Inf)<=TolG)
            exitCode=1;
            return;
        end
        %The 1/9 is an arbitrary factor to shrink the mu parameter.
        mu=mu*max(1/9,1-(2*rho-1)^3);
        nu=2;
    else%The step was not acceptable.
        mu=mu*nu;
        nu=2*nu;
    end
    
    %Check for numerical problems
    if(~isfinite(mu)||~isfinite(nu))
        exitCode=-2;
        x=[];
        return;
    end
end
    
%If we are here, then it reached the maximum number of iterations without
%convergence.
exitCode=0;
end

function [h,mu]=solveForH(A,g,mu,maxTries)
%This finds h in Equation 3.13. That is, it just solves (A+mu*I)*h=-g
%where A=J'*J is the squared Jacobian. However, one cannot directly solve
%the equation, because it is possible that mu is so small that (A+mu*I) is
%nearly singular. Thus, this tries to 

%How h is determined depends on A.
maxA=max(abs(A(:)));
if(maxA==0)%If A is a zero matrix.
    %This is the approximation for large values of mu as given on page 25.
    h=-g/mu;
    
    %Check for numerical problems and return empty matrices if there was a
    %problem.
    if(~isfinite(h))
       h=[];
       mu=[];
       return;
    end
else
    N=size(A,1);

    foundMu=false;
    for curTry=1:maxTries
        %Perform a Cholesky decompoostion without generating an error if
        %the matrix in question is singular. If the matrix is singular,
        %very close to singular, or contains NaNs, then p will be nonzero.
        [R,p]=chol(A+mu*eye(N),'upper');
        
        %If p is zero, make sure that the matrix is not nearly singular.
        %The "nearly singular" designation is based on the inverse
        %condition number being larger than the amount that would make
        %Matlab generate a warning when attempting to invert the matrix.
        if(p==0)
            p=(rcond(R)<eps);
            
            if(p==0)
                foundMu=true;
                break;
            end
        end

        %If we are here, then p~=0 and we will have to increase mu to make
        %the matrix nonsingular. The scaling factor of 10 is suggested in
        %[1].
        mu=max(10*mu,eps*maxA);
    end

    %If the above loop did not converge in a reasonable number of
    %iterations, then return empty matrices to indicate a failure.
    if(foundMu==false)
        h=[];
        mu=[];
        return;
    end
    
    %If we get here, then we can assume that R is nonsingular and we can
    %solve Equation 3.9, which is in a form using the square root of
    %(A+mu*I). Specifically, we are solving (R'*R)h=-g.
    
    h=-R\(R'\(g));
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
