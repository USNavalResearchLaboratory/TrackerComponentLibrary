function [xMin,fMin,exitCode]=quasiNetwonBFGS(f,x0,D0,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter)
%%QUASINEWTONBFGS Perform unconstrained nonlinear optimization using the
%                 Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm. The
%                 algorithm performs unconstrained minimization of a
%                 nonlinear function without one having to provide a
%                 Hessian matrix. Note that if a minimization over a
%                 least-squares problem is desired, then the
%                 Levenberg-Marquardt algorithm in LSEstLMarquardt is
%                 often preferable. The limited memory version of thie
%                 algorithm, quasiNewtonLBFGS, is more appropriate for use
%                 with very large matrices. The L-BFGS algorithm also
%                 supports an additional L1 norm term in the objective
%                 function. If one wishes to zero a vector (not a scalar),
%                 then NewtonsMethod is more appropriate as this function
%                 assumes the Hessian matrix is symmetric.
%
%INPUTS: f A handle to the function (and its gradient) over which the
%          minimization is to be performed. The function [fVal,gVal]=f(x)
%          takes the NX1 x vector and returns the real scalar function
%          value fVal and gradient gVal at the point x.
%       x0 The NX1-dimensional point from which the minimization starts.
%       D0 An estimate of the inverse Hessian matrix at x0. If omitted or
%          an empty matrix is passed, then the identity matrix is used.
%  epsilon The parameter determining the accuracy of the desired solution
%          in terms of the gradient. The function terminates when
%          norm(g) < epsilon*max([1, norm(x)])
%          where g is the gradient. The default if omitted or an empty
%          matrix is passed is 1e-6.
% deltaTestDist The number of iterations back to use to compute the
%          decrease of the objective function if a delta-based convergence
%          test is performed. If zero, then no delta-based convergence
%          testing is done. The default if omitted or an empty matrix is
%          passed is zero.
%    delta The delta for the delta convergence test. This determines the
%          minimum rate of decrease of the objective function. Convergence
%          is determined if (f'-f)<=delta*f, where f' is the value of the
%          objective function f deltaTestDist iterations ago,and f is the
%          current objective function value. The default if this parameter
%          is omitted or an empty matrix is passed is 0.
% lineSearchParams An optional structure whose members specify tolerances
%          for the line search. The parameters are described as in the
%          lineSearch function.
%   scaleD A boolean parameter indicating whether the inverse Hessian
%          estimate should be scaled as in Equation 1.201 of Section 1.7 of
%          [1]. The default if omitted or an empty matrix is passed is
%          false.
%  maxIter The maximum number of iterations to use for the overall BFGS
%          algorithm. The default if this parameter is omitted or an empty
%          matrix is passed is 1000.
%
%OUTPUTS: xMin The value of x at the minimum point found. If exitCode is
%              negative, then this might be an empty matrix.
%         fMin The cost function value at the minimum point found. If
%              exitCode is negative, then this might be an empty
%              matrix.
%     exitCode A value indicating the termination condition of the
%              algorithm. Nonnegative values indicate success; negative
%              values indicate some type of failure. Possible values are:
%                  0 The algorithm termiated successfully based on the
%                    gradient criterion.
%                  1 The algorithm terminated successfully based on the
%                    accuracy criterion.
%                 -1 The maximum number of overall iterations was reached.
%              Other negative values correspond to a failure in lineSearch
%              and correspond to the exitCode returned by the lineSearch
%              function. 
%
%The algorithm is implemented based on the description in Chapter 1.7 of
%[1] with the function lineSearch used to perform the line search.
%
%EXAMPLE 1:
%The first example is that used in the lineSearch file. 
% f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
%            [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
%            (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
% %Note that the deal function is used to make an anonymous function have
% %two outputs.
% x0=[0.5;0.25];
% [xMin,fMin,exitCode]=quasiNetwonBFGS(f,x0)
%The optimum point found is such that sum(xMin) is approximately
%1.73539450 with a minimum function value of approximately -2.1860756.
%
%EXAMPLE 2:
%The second example requires that the cost function be placed in a separate
%file. The cost function is
% function [fx,g]=objFun(x)
%    n=length(x);
%    fx=0;
%   
%    g=zeros(n,1);
%    for i=0:(n-1)
%        if(mod(i,2)==1)
%            continue;
%        end
%        t1=1-x(i+1);
%        t2=10*(x(i+1+1)-x(i+1)^2);
%        g(i+1+1)=20*t2;
%        g(i+1)=-2*(x(i+1)*g(i+1+1)+t1);
%        fx =fx+t1^2+t2^2;
%    end
% end
%It is used as
% n=100;
% x0=zeros(n,1);
% for i=0:(n-1)
%    if(mod(i,2)==1)
%        continue;
%    end
%
%    x0(i+1)=-1.2;
%    x0(i+1+1)=1;
% end
% f=@(x)objFun(x);
% [xMin,fMin,exitCode]=quasiNetwonBFGS(f,x0)
%whereby the optimal solution is all ones with a minimum function value of
%zero. This second example is the same as that provided with the L-BFGS
%library in C.
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 2nd ed. Belmont, MA: Athena
%    Science, 1999.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(x0,1);
if(nargin<3||isempty(D0))
    D=eye(xDim,xDim);
else
    D=D0;
end

if(nargin<4||isempty(epsilon))
    epsilon=1e-6;
end

if(nargin<5||isempty(deltaTestDist))
    deltaTestDist=0;
end

if(nargin<6||isempty(delta))
    delta=0;
end

if(nargin<7)
    lineSearchParams=[];
end

if(nargin<8||isempty(scaleD))
   scaleD=false; 
end

if(nargin<9||isempty(maxIter))
   maxIter=1000; 
end

xPrev=x0;
[fValPrev,gradFPrev]=f(x0);

if(deltaTestDist>0)
    pastFVals=zeros(deltaTestDist,1);
    pastFVals(1)=fValPrev;
end
for curIter=1:maxIter
    %Equation 1.181 for the descent direction.
    d=-D*gradFPrev;
    
    %Perform a line search in the given descent direction.
    [xCur,fValCur,gradFCur,~,~,exitCode]=lineSearch(f,xPrev,d,[],[],lineSearchParams);

    if(isempty(xCur))
        xMin=[];
        fMin=[];
        return;
    end
 
    %Check for convergence based on the gradient.
    if(norm(gradFCur)<epsilon*max([1, norm(xCur)]))
        xMin=xCur;
        fMin=fValCur;
        exitCode=0;
        return;
    end
    
    %Check for convergence based on the actual function value.
    if(deltaTestDist~=0)
        if(pastFVals(end)-fValCur<=delta*fValCur)
            xMin=xCur;
            fMin=fValCur;
            exitCode=1;
            return; 
        end
        pastFVals=circshift(pastFVals,[1,0]);
        pastFVals(1)=fValCur;
    end

    %After taking the step, update the inverse Hessian approximation.
    D=updateHessianApprox(D,xCur,xPrev,gradFCur,gradFPrev,scaleD);

    xPrev=xCur;
    gradFPrev=gradFCur;
end

%The maximum number of iterations elapsed without convergence
xMin=xCur;
fMin=fValCur;
exitCode=-1;
end

function D=updateHessianApprox(D,xCur,xPrev,gradFCur,gradFPrev,scaleD)
%This is the update for the inverse Hessian estimate in the BFGS algorithm.
%All equation numbers refer to Section 1.7 of [1].

q=gradFCur-gradFPrev;%Equation 1.184
p=xCur-xPrev;%Equation 1.183

%A test is added so that D does not change if p'*q is very close to zero.
%Otherwise, D would just be full of NaN and Inf values. The extra realmin
%test avoids denormalized numbers.
if(p'*q>max(realmin,max(eps(p),eps(q))))
    if(scaleD==true)
    %If D should be scaled as in Equation 1.201. As noted, this can improve
    %the condition number of D and is often recommended after the first
    %iteration.
        D=(p'*q/(q'*D*q))*D;
    end

    %The BFGS update from Problem 1.7.2 with the correction from the errata
    %for the denominator of the p*p' term. 
    D=D+(1+q'*D*q/(p'*q))*(p*p')/(p'*q)-(D*q*p'+p*q'*D)/(p'*q);
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
