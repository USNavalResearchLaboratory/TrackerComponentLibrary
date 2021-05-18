function [xMin,fMin,exitCode]=NewtonsMethod(f,fHess,x0,epsilon,deltaTestDist,delta,lineSearchParams,cholSemiVal,maxIter)
%%NEWTONSMETHOD Perform unconstrained nonlinear optimization using Newton's
%               method. Unlike quasiNetwonBFGS, the Hessian matrix must be
%               provided. However, while quasiNetwonBFGS can only minimize
%               a function, NewtonsMethod can also zero a vector that is
%               given as the "gradient" in which case the "Hessian" is a
%               Jacobian. While a real Hessian must be symmetric, this
%               Jacobian need not be symmetric.
%
%INPUTS: f A handle to the function and its gradient over which the
%          minimization is to be performed. The function
%          [fVal,grad,Hess]=f(x) takes the NX1 x vector returns the real
%          scalar function value fVal, and the gradient grad, at the point
%          x. If the Hessian is not given separately by fHess, then it can
%          be given as a third output and fHess set to an empty matrix. If
%          one wishes to perform diagonal loading, it should be already in
%          Hess when provided to this function.
%    fHess A handle to the Hessian of the function f. This is only needed
%          if f does not return the Hessian as a third output. This must be
%          positive definite for the algorithm to work. As shown in the
%          example code below, the "Hessian" can actually be a (non-
%          symmetric) Jacobian matrix when the optimization is meant to
%          zero a vector and not minimize a function. If an empty matrix is
%          passed, then it is assumed that f returns the Hessian after the
%          function value and gradient. If one wishes to perform diagonal
%          loading, it should be built into fHess.
%       x0 The NX1-dimensional point from which the minimization starts.
%  epsilon The parameter determining the accuracy of the desired solution
%          in terms of the gradient. This can be a 2X1 or 1X2 vector, where
%          epsilon(1) is an absolute tolerance and epsilon(2) is a relative
%          tolerance, or a single epsilon can be passed for both values to
%          be set equal. Letting maxVal=max(abs(gradFCur)), the function
%          terminates when maxVal<=AbsTol or maxVal<=RelTol*max(abs(xMin))
%          The default if omitted or an empty matrix is passed is 1e-6.
% deltaTestDist The number of iterations back to use to compute the
%          decrease of the objective function if a delta-based convergence
%          test is performed. If zero, then no delta-based convergence
%          testing is done. The default if omitted or an empty matrix is
%          passed is zero.
%    delta The delta for the delta convergence test. This determines the
%          minimum rate of decrease of the objective function. Convergence
%          is determined if
%          (f'-f)<=delta*f, where f' is the value of the
%          objective function f deltaTestDist iterations ago, and f is the
%          current objective function value. The default if this parameter
%          is omitted or an empty matrix is passed is 0.
% lineSearchParams An optional structure whose members specify tolerances
%          for the line search. The parameters are described as in the
%          lineSearch function, except if -1 is passed instead of a
%          parameter structure, then no line search is performed. The line
%          search must often be disabled if one is interested in zeroing
%          a vector (the gradient term being the vector) rather than
%          actually minimizing a function.
% cholSemiVal A value indicating whether a method similar to that suggested
%          in Chapter 1.4 of [1] for using a modified Cholesky
%          decomposition of the Hessian matrix should be used so as to
%          assure that one always obtains a descent direction. If a
%          positive value <1 of cholSemiVal is given, then the Hessian
%          matrix is decomposed with cholSemiDef with epsVal of
%          cholSemiVal. If a negative value is provided, then no
%          decomposition is done and no effort is made to assure that the
%          direction traveled is a descent direction. The default if this
%          parameter is omitted or an empty matrix is passed is 1e-6, This
%          parameter assumes the "Hessian" matrix is symmetric. Thus, this
%          option should be set to -1 if this function is used to zero a
%          vector rather than minimize a function.
%  maxIter The maximum number of iterations to use for the algorithm. The
%          default if this parameter is omitted or an empty matrix is
%          passed is 100.
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
%              -1023 A logical error in the line search code occurred.
%              -1001 A finite precision error occurred or no line-search
%                    step satisfies the sufficient decrease and curvature
%                    conditions.
%              -1000 The line-search step size became less than minStep.
%               -999 The line-search step size became larger than maxStep.
%               -998 The maximum number of line search iterations was
%                    reached.
%               -996 The relative width of the interval of uncertainty in
%                    the line search algorithm is at most xTol.
%               -995 A negative line-search step occurred.
%               -994 The current search direction increases the objective
%                    function.
%
%The algorithm is implemented based on the description in Chapter 1.2 of
%[1].
%
%EXAMPLE 1:
%In this example, we want to use Newton's method to zero a vector and we do
%not care about any scalar function value. Thus, the convergence test
%using the function value is disabled by setting deltaTestDist=0. The line
%search must be disabled for this to work.
% fUnused=@(x)0;
% gVec=@(x)[x(1)^2+x(2)^2+x(3)^2-3;
%           x(1)^2+x(2)^2-x(3)-1;
%           x(1)+x(2)+x(3)-3];
% Jacob=@(x)[2*x(1), 2*x(2), 2*x(3);
%            2*x(1), 2*x(2), -1;
%            1,   1,   1];
% f=@(x)deal(fUnused(x),gVec(x));
% deltaTestDist=0;%Ignore the function value when considering convergence.
% lineSearchParams=-1;%Do not perform a line search.
% epsilon=1e-8;%Tolerance on gVec.
% x0=[1;2;3];%Initial estimate.
% cholSemiVal=-1;%Must be negative as Jacob is not symmetric.
% [xMin,fMin,exitCode]=NewtonsMethod(f,Jacob,x0,epsilon,deltaTestDist,[],lineSearchParams,cholSemiVal)
%
%EXAMPLE 2:
%This is the same as the above example, but in this case, we zero the
%vector by minimizing the squared magnitude of vector. Of course, such an
%approach can introduce many new local minima and is often not ideal.
% func=@(x)(x(1)^2+x(2)^2+x(3)^2-3)^2+(x(1)^2+x(2)^2-x(3)-1)^2+(x(1)+x(2)+x(3)-3)^2;
% grad=@(x)[2*(-3+4*x(1)^3+x(2)+x(3)+x(1)*(-7+4*x(2)^2+2*(x(3)-1)*x(3)));
%           2*(-3+x(1)+4*x(1)^2*x(2)+4*x(2)^3+x(3)+x(2)*(-7+2*(-1+x(3))*x(3)));
%           2*(-2+x(1)+x(2)-4*x(3)+2*x(3)^3+x(1)^2*(-1+2*x(3))+x(2)^2*(-1+2*x(3)))];
% f=@(x)deal(func(x),grad(x));
% Hess=@(x)vech2Mat([-14+24*x(1)^2+8*x(2)^2+4*(-1+x(3))*x(3);
%       2+16*x(1)*x(2);
%       2+x(1)*(-4+8*x(3));
%       -14+8*x(1)^2+24*x(2)^2+4*(-1+x(3))*x(3);
%       2+x(2)*(-4+8*x(3));
%       4*(-2+x(1)^2+x(2)^2+3*x(3)^2)]);
% epsilon=1e-8;%Tolerance on gVec.
% x0=[1;2;3];%Initial estimate.
% [xMin,fMin,exitCode]=NewtonsMethod(f,Hess,x0,epsilon)
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 2nd ed. Belmont, MA: Athena
%    Science, 1999.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<4||isempty(epsilon))
        epsilon=1e-6;
    end

    if(length(epsilon)==2)
        AbsTol=epsilon(1);
        RelTol=epsilon(2);
    else
        AbsTol=epsilon;
        RelTol=epsilon;
    end

    if(nargin<5||isempty(deltaTestDist))
        deltaTestDist=0;
    end

    if(nargin<6||isempty(delta))
        delta=0;
    end

    %Take any line search parameters given.
    useLineSearch=true;
    if(nargin>=7&&~isempty(lineSearchParams))
        if(isnumeric(lineSearchParams)&&lineSearchParams==-1)
            %If one just wishes to blindly take steps.
            useLineSearch=false;
        elseif(~isstruct(lineSearchParams)||(isnumeric(lineSearchParams)&&lineSearchParams~=-1))
            error('Unknown lineSearchParams value given')
        end
    else
        lineSearchParams=[];
    end
    
    if(nargin<8||isempty(cholSemiVal))
        cholSemiVal=1e-6;
    elseif(cholSemiVal>1)
        error('cholSemiVal should be less than 1')
    end
    
    if(nargin<9||isempty(maxIter))
        maxIter=100; 
    end
    
    xPrev=x0;
    if(isempty(fHess))
        [fValPrev,gradFPrev,fHessVal]=f(x0);
    else
        [fValPrev,gradFPrev]=f(x0);
        fHessVal=fHess(xPrev);
    end

    if(deltaTestDist>0)
        pastFVals=zeros(deltaTestDist,1);
        pastFVals(1)=fValPrev;
    end
    
    %Options for upper and lower triangular matrices for use in the
    %linsolve function if cholSemiVal>0.
    optsLT.LT=true;
    optsLT.UT=false;
    optsUT.UT=true;
    optsUT.LT=false;

    for cutIter=1:maxIter
        if(cholSemiVal>0)
            L=cholSemiDef(fHessVal,'lower',0,cholSemiVal);
            
            if(any(~isfinite(L(:))))
                exitCode=-2;
                xMin=[];
                fMin=[];
                return;
            end
            
            %Get the descent direction, taking advantage of the
            %lower-triangular structure of L.
            D=-linsolve(L',linsolve(L,gradFPrev,optsLT),optsUT);
        else
            %The descent direction to use. This option will work with
            %non-symmetric fHessVal, but will not work if fHessVal is not
            %positive definite.
            D=-fHessVal\gradFPrev;
        end

        if(useLineSearch)
            %Perform a line search in the given descent direction.
            [xCur,fValCur,gradFCur,~,~,exitCode]=lineSearch(f,xPrev,D,[],[],lineSearchParams);

            if(isempty(xCur))
                xMin=[];
                fMin=[];
                return;
            end
            
            %If the Hessian is provided by the function, not by a separate
            %function.
            if(isempty(fHess))
                [fValCur,gradFCur,fHessVal]=f(xCur);
            else
                fHessVal=fHess(xCur);
            end
        else%Take a step without a line search.
            xCur=xPrev+D;
            if(isempty(fHess))
                [fValCur,gradFCur,fHessVal]=f(xCur);
            else
                [fValCur,gradFCur]=f(xCur);
                fHessVal=fHess(xCur);
            end
        end
        
        if(any(~isfinite(xCur(:))))
            exitCode=-2;
            xMin=[];
            fMin=[];
            return;
        end
        
        %Check for convergence based on the gradient.
        maxVal=max(abs(gradFCur));
        if(maxVal<=AbsTol||maxVal<=RelTol*max(abs(xCur)))
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
    
        xPrev=xCur;
        gradFPrev=gradFCur;
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
