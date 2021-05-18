function [xMin,fMin,exitCode]=conjugateGradMethod(f,x0,epsilon,deltaTestDist,delta,lineSearchParams,resetEvery,gammaVal,maxIter)
%%CONJUGATEGRADMETHOD Perform unconstrained nonlinear optimization using
%              the conjugate gradient algorithm. The algorithm performs
%              unconstrained minimization of a nonlinear function without
%              having to provide a Hessian matrix. On large problems, this
%              algorithm can be faster than Newton's method, which has a
%              large matrix inversion, and steepest ascent, which usually
%              converges more slowly.
%
%INPUTS: f A handle to the function (and its gradient) over which the
%          minimization is to be performed. The function [fVal,gVal]=f(x)
%          takes the NX1 x vector and returns the real scalar function
%          value fVal and NX1 gradient gVal at the point x.
%       x0 The NX1-dimensional point from which the minimization starts.
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
%          lineSearch function, except if -1 is passed instead of a
%          parameter structure, then no line search is performed.
% resetEvery In Chapter 2.1 of [1], when handling non-quadratic function,
%          it is suggested that the subgradients be reset every N
%          iterations. By reset, this means that the past direction is
%          discarded and the next step is a steepest descent step. This is
%          the number of iterations after which a reset will occur.
%          resetEvery>=1. The default if omitted or an empty matrix is
%          passed is N.
% gammaVal In Chapter 2.1 of [1], it is suggested that the descent
%          direction algorithm be reset if
%          abs(g'*gPrev)>gammaVal*norm(gPrev)^2
%          where 0<gammaVal<1. The notion is that the previous and next
%          gradient directions should be nearly orthogonal, so if they are
%          not, within some tolerance, then a steepest descent step is
%          taken. The default if this parameter is omitted or an empty
%          matrix is passed is 0.1.
%  maxIter The maximum number of iterations to use for the algorithm. The
%          default if this parameter is omitted or an empty matrix is
%          passed is 4*N.
%
%OUTPUTS: xMin The value of x at the minimum point found. If exitCode is
%              negative, then this might be an empty matrix.
%         fMin The cost function value at the minimum point found. If
%              exitCode is negative, then this might be an empty matrix.
%     exitCode A value indicating the termination condition of the
%              algorithm. Nonnegative values indicate success; negative
%              values indicate some type of failure. Possible values are:
%                  0 The algorithm termiated successfully based on the
%                    gradient criterion.
%                  1 The algorithm terminated successfully based on the
%                    accuracy criterion.
%                 -1 The maximum number of overall iterations was reached.
%                 -2 A non-finite value was encoutered outside the line
%                    search.
%              Other negative values correspond to a failure in lineSearch
%              and correspond to the exitCode returned by the lineSearch
%              function. 
%
%The algorithm is implemented based on the description in Chapter 2.1 of
%[1].
%
%EXAMPLE:
%The example is that used in the lineSearch file. 
% f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
%            [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
%            (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
% %Note that the deal function is used to make an anonymous function have
% %two outputs.
% x0=[0.5;0.25];
% [xMin,fMin,exitCode]=conjugateGradMethod(f,x0)
%The optimum point found is such that sum(xMin) is approximately
%1.73539450 with a minimum function value of approximately -2.1860756.
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 3rd ed. Belmont, MA: Athena
%    Science, 2016.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

    if(nargin<3||isempty(epsilon))
        epsilon=1e-6;
    end

    if(nargin<4||isempty(deltaTestDist))
        deltaTestDist=0;
    end

    if(nargin<5||isempty(delta))
        delta=0;
    end

    %Take any line search parameters given.
    useLineSearch=true;
    if(nargin>=6&&~isempty(lineSearchParams))
        if(isnumeric(lineSearchParams)&&lineSearchParams==-1)
            %If one just wishes to blindly take steps.
            useLineSearch=false;
        elseif(~isstruct(lineSearchParams)||(isnumeric(lineSearchParams)&&lineSearchParams~=-1))
            error('Unknown lineSearchParams value given')
        end
    else
        lineSearchParams=[];
    end
    
    N=size(x0,1);
    
    if(nargin<7||isempty(resetEvery))
       resetEvery=N; 
    end
    
    if(nargin<8||isempty(gammaVal))
       gammaVal=0.1; 
    end

    if(nargin<9||isempty(maxIter))
        maxIter=4*N; 
    end

    x=x0;
    [fVal,gradF]=f(x0);
    
    %The previous values are not yet set.
    gradFPrev=[];
    dPrev=[];
    
    if(deltaTestDist>0)
        pastFVals=zeros(deltaTestDist,1);
        pastFVals(1)=fVal;
    end
    
    stepsSinceReset=0; 

    for cutIter=1:maxIter
        if(stepsSinceReset>0)
            if(stepsSinceReset>=resetEvery)
                stepsSinceReset=0;
            else
                gradFPrevMag2=gradFPrev'*gradFPrev;

                if(abs(gradF'*gradFPrev)>gammaVal*gradFPrevMag2)
                    stepsSinceReset=0;
                end
            end
        end

        if(stepsSinceReset>0)
            %If it hasn't just been reset.

            %Equation 2.35
            beta=gradF'*(gradF-gradFPrev)/gradFPrevMag2;
            
            %Equation 2.34 for the descent direction.
            d=-gradF+beta*dPrev;
        else
            %Steepest descent direction.
            d=-gradF;
        end

        if(useLineSearch)
            %Perform a line search in the given descent direction.
            [xCur,fValCur,gradFCur,~,~,exitCode]=lineSearch(f,x,d,[],[],lineSearchParams);
            if(isempty(xCur))
                xMin=[];
                fMin=[];
                return;
            end
        else%Take a step without a line search.
            xCur=xPrev+d;
        end
        
        if(any(~isfinite(xCur(:))))
            exitCode=-2;
            xMin=[];
            fMin=[];
            return;
        end
        
        %Check for convergence based on the gradient.
        if(norm(gradFCur)<epsilon*max([1,norm(xCur)]))
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
            pastFVals(1)=fVal;
        end
        
        x=xCur;
        gradFPrev=gradF;
        gradF=gradFCur;
        dPrev=d;
        stepsSinceReset=stepsSinceReset+1;
    end
    
    %The maximum number of iterations elapsed without convergence
    xMin=xCur;
    fMin=fValCur;
    exitCode=-1;
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
