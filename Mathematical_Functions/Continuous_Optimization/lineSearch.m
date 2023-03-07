function [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearch(func,x,d,fInit,fGradInit,lineSearchParams)
%%LINESEARCH Given a function f taking an NX1 dimensional vector x and a
%            descent direction d, find the positive value alpha such that
%            f(x+alpha*d) is minimized to an extent necessary for
%            satisfying certain convergence criterion (Armijo's rule and
%            possibly Wolfe conditions) that are useful in certain
%            optimization methods, such as quasi-Newton methods.
%
%INPUTS: func A handle to the function (and its gradient) over which the line
%          search minimization is to be performed. The function
%          [fVal,fGrad]=f(x) takes the NX1 x vector returns the real scalar
%          function value fVal and NX1 gradient fGrad at the point x. 
%        x The NX1-dimensional point from which the line search starts.
%        d The NX1 vector defining the initial descent direction. The
%          minimization is performed such that f(x+alpha*d) is minimized
%          where alpha is a scalar parameter.
% fInit, fGradInit The value of func and its gradient at the point x. If
%          these are omitted or empty matrices are, then func is
%          called to obtain these values.
% lineSearchParams A structure whose members determine the accuracy and
%          performance of the algorithms. Possible members names
%          (used, for example, as linSearchParams.fTol) are:
%          algorithm A parameter specifying the line search algorithm to use.
%             Possible values are:
%             0) (The default if omitted or an empty matrix is passed) Use
%               the algorithm proposed by Mor� and Thuente in [1] (modified 
%               as described below). Unless one sets
%               lineSearchParams.maxStepSize=stepSize, this assumes that
%               one has not obtained initial bounds on the stepsize search 
%               region and it firsts finds an upper bound (possibly by
%               expanding the search region by successive factors of
%               initStepMult). Then, quadratic and cubic interpolation are
%               used to reduce the size of the region until the specified
%               conditions are satisfied.
%             1) Use a backtracking algorithm to just satisfy Armijo's
%               rule. Armijo's rule is originally from [2] and in places
%               such at Chapter 1.2 of [5], the expression of the rule
%               implies the backtracking algorithm. The input stepSize is
%               the maximum and initial step size tested. 
%             2) Assume that stepSize is the upper bound on the step and
%               use a bisection approach to satisfy the conditions. The
%               rules for branching in the bisection algorithm are the same
%               as those in [1], where more sophisticated techniques are
%               used.
%          stepSize The initial step size used in the algorithms. The
%               default is 1, which is generally suitable for quasi-Newton
%               methods.
%          fTol The tolerance used for declaring convergence using
%               Armijo's rule. This also has to be satisfied with the Wolfe
%               conditions. Convergence is declared if
%               func(x+alpha*d)<=fInit+fTol*alpha*fGradInit'*d
%               The default value is 1e-6. Note that 0<=fTol<=1.
%          wolfeTol The tolerance used for the weak and/or strong Wolfe
%               conditions for declaring convergence. The strong Wolfe
%               condition is satisfied if
%               abs(fGrad(x+alpha*d)'*d)<=wolfeTol*abs(fGrad'*d)
%               The regular Wolfe condition is satisfied if 
%               fGrad(x+alpha*d)'*d>=-wolfeTol*abs(fGrad'*d)
%               The default is 0.9. Note that 0<=wolfeTol<=1. The Wolfe
%               conditions are from [3] and [4].
%          RelTolAlpha A relative tolerance used in algorithms 0 and 2 such
%               that if the bounded search region is less than RelTolAlpha
%               time the lower bound, then convergence is declared. The
%               default is eps(1);
%          conditionType In algorithms 0 and 2, one can choose different
%               convergence conditions. 0 is only Armijo's rule. 1 is the
%               regular Wolfe condition and 2 is the strong WOlfe
%               condition. For 1 and 2 Armijos' rule must also be
%               satisfied. The default is 2.
%          maxIter The maximum number of function evaluations (not counting
%               obtaining fInit, fGradInit if not provided) that are
%               allowed. The default is 20.
%          minStepSize The smallest that the step is allowed to get. The
%               default is 1e-20.
%          maxStepSize An upper bound on the step size that is only used
%               when algorithm=2. The default is 1e20.
%          initStepMult When algorithm=0 and an upper bound has not been
%               obtained, this is the multiple by which the upper end of
%               the search region is increased. the default if omitted is
%               4. This value must be >1.
%          scalFact The factor by which the search region in the
%               backtracking algorithm is shrunk each step. The default is
%               1/2. This must be >0 and <1.
%
%OUTPUTS: xMin The value of x at the minimum point found. This and the
%              following 3 outputs will be empty matrices if the algorithm
%              has any of the failure cases described in exitCode.
%         fMin The cost function value at the minimum point found.
%         gMin The value of the gradient of the cost function at the
%              minimum point found.
%        alpha The step size found such that xMin=x+alpha*d
%     exitCode A value indicating the termination condition of the
%              algorithm. Negative values indicate a lack of convergence,
%              but except in the case of certain failures, the results are
%              still usable.
%              0 or a positive number: The algorithm terminated
%                successfully. The number is the number of iterations
%                required. 
%              -1001 A complex or non-finite value was encountered. This is
%                    often a failure leading to the other outputs being
%                    empty matrices.
%              -1000 The line search step size was less than
%                    lineSearchParams.minStepSize.
%              -998 The maximum number of iterations was reached.
%              -995 The search direction given by d is not a descent
%                   direction (it increases the cost function). This is a
%                   failure that will lead to the other outputs being empty
%                   matrices.
%              -994 The search direction at the initial point is not a
%                   descent direction. This is a failure that will lead to
%                   the other outputs being empty matrices.
%
%Line searches play a pivotal role in multivariate optimization routines.
%They are related to single variable optimization methods, such as in the
%function goldenSectionSearch or in Matlab fminbnd function, in that they
%try to minimize a function over one parameter (in this case a
%multivariate function being minimized over a line). However, line search
%methods generally try less to bound the minimum than to determine a point
%that sufficiently decreases the objective function to assure convergence
%of the overall multivariate optimization problem over time.
%
%Thus, outside of the context of larger multivariate optimization
%problems, line search algorithms are not always very good at performing
%univariate optimization. On the other hand, univariate optimization
%algorithms can usually be used for performing a line search, because with
%many objective function, one can often use an initial step size of 1 as
%the maximum step size and then find the minimum over the bounded line.
%
%The implementation fo the More-Thuente algorithm (algorithm 0) here
%differs from the description in [1] in that it forces a bisection step
%after only a single failure to decrease the step size. The paper did not
%give a pseudocode listing, so a number of trivial details were filled in
%to create a full algorithm.
%
%EXAMPLE 1:
%As an example, we will take the sample problem used in goldenSectionSearch
%and turn it into a bivariate minimization problem where a particular line
%over which we want to minimize is the same as in the goldenSectionSearch
%example. This will let us see how a line search algorithms helps minimize
%a function, but does not necessarily find the global minimum unless we use
%the correct algorithm and adjust the default aprameters.
% f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
%             [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
%             (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
% %Note that the deal function is used to make an anonymous function have
% %two outputs.
% x=[4;0];
% d=[1;0];
% [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearch(f,x,d)
% %In the above example, the starting point is x=[4;0] and only the x
% %coordinate changes due to the direction of d, so the y coordinate can be
% %effectively ignored. The local minimizing point found is
% %x=[6.091159459757882;0] with a minimum function value of
% %fMin=0.048242331317758. However, this is not the absolute local minimum.
% %The absolute local minimum is at x=[6;0], where fMin=0. Thus, unlike
% %bounding a region and using goldenSectionSearch, this function decreases
% %the cost sufficient according to certain conditions, but it does not find
% %a local minimum. However, if we change the convergence conditions:
% lineSearchParams=[];
% lineSearchParams.algorithm=0;
% lineSearchParams.fTol=0;
% lineSearchParams.wolfeTol=0;
% [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearch(f,x,d,[],[],lineSearchParams)
%One will find that termination was due to hitting the maximum number of
%iteration (20) and that x=[5.999999996806924;0] and
%fMin=6.736154520036488e-32, which is very near the globally optimal
%solution.
%
%EXAMPLE 2:
%In this example, we minimize a scalar function but the initial step size
%used in actually less than what woud be ideal. This is an instance where
%algorithm 0 can be useful, since it searched for an appropriate upper
%bound. Here, when run with tight tolerances to get something close to the
%globally optimal value from a single line search.
% b=2;
% func=@(a)deal(-a/(a^2+b),(a^2-b)/(a^2+b)^2);
% d=1;
% x=0;
% lineSearchParams=[];
% lineSearchParams.algorithm=0;
% lineSearchParams.stepSize=1e-3;
% lineSearchParams.fTol=0;
% lineSearchParams.wolfeTol=0;
% [fInit,fGradInit]=func(x);
% [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearch(func,x,d,fInit,fGradInit,lineSearchParams)
%One will find that xMin is close to sqrt(2) (the minimum value), with a
%difference on the order of 1e-9. sqrt(2) is larger than the initial step
%size. Termination was due to hitting the maximum number of iterations. If
%algorithm 2 were used with the same step size, then it would just return
%the initial step given, because it could never find anything better.
%
%REFERENCES:
%[1] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%[2] L. Armijo, "Minimization of functions having Lipschitz continuous
%    first partial derivatives," Pacific Journal of Mathematics, vol. 16,
%    no. 1, pp. 1-3, Nov. 1966.
%[3] P. Wolfe, "Convergence conditions for ascent methods," SIAM Review,
%    vol. 11, no. 2, pp. 226-235, Apr. 1969.
%[4] P. Wolfe, "Convergence conditions for ascent methods. II: Some
%    corrections," SIAM Review, vol. 13, no. 2, pp. 185-188, Apr. 1971.
%[5] D. P. Bertsekas, Nonlinear Programming, 3nd ed. Belmont, MA: Athena
%    Science, 2016.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(fInit)||nargin<5||isempty(fGradInit))
    [fInit,fGradInit]=func(x);
    
    if(~isfinite(fInit)||~isreal(fInit)||any(~isfinite(fGradInit))||~isreal(fGradInit))
        xMin=[];
        fMin=[];
        fGradMin=[];
        gMin=[];
        alpha=[];
        exitCode=-1001;
        return
    end
end

algorithm=0;
stepSize=1;
maxStepSize=1e20;%The default for algorithm=0.
minStepSize=1e-20;
maxIter=20;
wolfeTol=0.9; 
fTol=1e-6;
RelTolAlpha=eps();
conditionType=2;
initStepMult=4;
scalFact=1/2;

%Extract any parameters if provided.
if(nargin>5&&~isempty(lineSearchParams))
    if(isfield(lineSearchParams,'algorithm')&&~isempty(lineSearchParams.algorithm))
        algorithm=lineSearchParams.algorithm;
    end

    if(isfield(lineSearchParams,'stepSize')&&~isempty(lineSearchParams.stepSize))
        stepSize=lineSearchParams.stepSize;
    end

    if(algorithm~=0)
       maxStepSize=stepSize;%The default for algorithm~=0.
    end

    if(isfield(lineSearchParams,'maxStepSize')&&~isempty(lineSearchParams.maxStepSize))
        minStepSize=lineSearchParams.maxStepSize;
    end
    
    if(isfield(lineSearchParams,'minStepSize')&&~isempty(lineSearchParams.minStepSize))
        minStepSize=lineSearchParams.minStepSize;
    end
    
    if(isfield(lineSearchParams,'maxIter')&&~isempty(lineSearchParams.maxIter))
        maxIter=lineSearchParams.maxIter;
    end

    if(isfield(lineSearchParams,'wolfeTol')&&~isempty(lineSearchParams.wolfeTol))
        wolfeTol=lineSearchParams.wolfeTol;

        if(wolfeTol<0||wolfeTol>1)
            error('It is required that 0<=wolfeTol<=1.') 
        end
    end

    if(isfield(lineSearchParams,'fTol')&&~isempty(lineSearchParams.fTol))
        fTol=lineSearchParams.fTol;
        
        if(fTol<0||fTol>1)
            error('It is required that 0<=fTol<=1.') 
        end
    end

    if(isfield(lineSearchParams,'RelTolAlpha')&&~isempty(lineSearchParams.RelTolAlpha))
        RelTolAlpha=lineSearchParams.RelTolAlpha;
    end

    if(isfield(lineSearchParams,'conditionType')&&~isempty(lineSearchParams.conditionType))
        conditionType=lineSearchParams.conditionType;
    end

    if(isfield(lineSearchParams,'initStepMult')&&~isempty(lineSearchParams.initStepMult))
        initStepMult=lineSearchParams.initStepMult;
        
        if(initStepMult<=1)
            error('It is required that initStepMult>1.')
        end
    end

    if(isfield(lineSearchParams,'scalFact')&&~isempty(lineSearchParams.scalFact))
        scalFact=lineSearchParams.scalFact;
        
        if(scalFact>=1)
            error('It is required that scalFact<1.')
        end
    end
end

if(stepSize<minStepSize)
   error('The provided stepSize is less than minStepSize.') 
end

switch(algorithm)
    case 0
        [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearchMoreThuente(x,func,fInit,fGradInit,d,stepSize,RelTolAlpha,fTol,wolfeTol,conditionType,maxIter, minStepSize,maxStepSize,initStepMult);
    case 1
        [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearchBacktracking(x,func,fInit,fGradInit,d,stepSize,fTol,scalFact,maxIter,minStepSize);
    case 2
        [xMin,fMin,fGradMin,gMin,alpha,exitCode]=lineSearchBisection(x,func,fInit,fGradInit,d,stepSize,RelTolAlpha,fTol,wolfeTol,conditionType,maxIter,minStepSize);
    otherwise
        error('Unknown algorithm specified.')
end
end

function [xT,f,fGrad,g,stepSize,exitCode]=lineSearchBacktracking(x,func,fInit,fGradInit,s,stepSize,fTol,scalFact,maxIter,minStepSize)
%%LINESEARCHBACKTRACKING This function implements the backtracking
%           algorithm implied in Chapter 1.2 of [1] to satisfy Armijo's
%           rule, which is derived in [2].
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 3nd ed. Belmont, MA: Athena
%    Science, 2016.
%[2] L. Armijo, "Minimization of functions having Lipschitz continuous
%    first partial derivatives," Pacific Journal of Mathematics, vol. 16,
%    no. 1, pp. 1-3, Nov. 1966.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

f=fInit;
xT=x;

%The initial gradient in the descent direction.
fGrad=fGradInit;
gInit=fGrad'*s;

%If the search direction is not a descent direction.
if(gInit>0)
    xT=[];
    f=[];
    fGrad=[];
    g=[];
    stepSize=[];
    exitCode=-994;
    return;
elseif(gInit==0)%If the  function is already at a vertex.
    g=gInit;
    stepSize=0;
    %The function could be stuck at a point of inflection.
    exitCode=0;
    return;
end

%If the initial step size is invalid.
if(stepSize<=0)
    xT=[];
    f=[];
    fGrad=[];
    g=[];
    stepSize=[];
    exitCode=-995;
    return
end

curIter=0;
while(1)
    xT=x+stepSize*s;
    
    [f,fGrad]=func(xT);
    if(~isfinite(f)||~isreal(f)||any(~isfinite(fGrad))||~isreal(fGrad))
        xT=[];
        f=[];
        fGrad=[];
        g=[];
        stepSize=[];
        exitCode=-1001;
        return
    end

    curIter=curIter+1;
    %Check Armijo's rule (Chapter 1.2 of [1], pg. 36)
    if(f>fInit+stepSize*fTol*gInit)
        scalFactor=scalFact;
    else
        g=fGrad'*s;
        exitCode=curIter;
        return;
    end

    if(stepSize<=minStepSize)
        g=fGrad'*s;
        exitCode=-1000;
        return;  
    end
    if(curIter>=maxIter)
        g=fGrad'*s;
        exitCode=-998;
        return
    end

    %Shink the step size. This is beta^m*s in the stepsize reduction rule
    %for Armijo's method in Chapter 1.2 of [1] (pg. 36).
    stepSize=stepSize*scalFactor;
end

end

function [xL,fL,fGradL,gL,alphaL,exitCode]=lineSearchMoreThuente(x,func,fInit,fGradInit,s,stepSize,RelTolAlpha,fTol,wolfeTol,conditionType,maxIter, minStepSize,maxStepSize,initStepMult)
%%LINESEARCHMORETHUENTE This function implements the Mor�-Thuente algorithm
%               of [1].
%
%REFERENCES:
%[1] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The current test step size.
alphaT=stepSize;

alphaMin=minStepSize;
alphaMax=maxStepSize;

%Make sure that the initial test step size satisfies the absolute bounds
%given.
alphaT=min(max(alphaT,alphaMin),alphaMax);

%We use the notation in the paper.
mu=fTol;
eta=wolfeTol;

%The values of the step, the functional value, and the directional
%derivative at the best step --they all end in L as in [1].
alphaL=0;
fL=fInit;
fGradL=fGradInit;
gL=fGradL'*s;%The initial gradient in the search direction.
xL=x;%The value of x at the best step thus far.

%The initial gradient is needed for the right-hand side of Equation 1.1 in
%[1].
gInit=gL;

%If the search direction is not a descent direction.
if(gInit>0)
    xL=[];
    fL=[];
    fGradL=[];
    gL=[];
    alphaL=[];
    exitCode=-994;
    return;
elseif(gInit==0)%If the  function is already at a vertex.
    %The function could be stuck at a point of inflection.
    exitCode=0;
    return;
end

%If the initial step size is invalid.
if(stepSize<=0)
    xL=[];
    fL=[];
    fGradL=[];
    gL=[];
    alphaL=[];
    exitCode=-995;
    return
end

%The values of the step, the function value,  at the other uncertainty
%endpoint is just set to the first endpoint to start, because we haven't
%determined a useful upper bound.
alphaU=alphaL;
fU=fL;
fGradU=fGradL;
gU=gL;
%Note that the magnitudes of alphaL and alphaU are not ordered. One is not
%always smaller or larger than the other.

%Initially, we do not have a meaningful upper bound on the step size. The
%upper bound is just maxStepsize (alphaMax). We will use different
%heuristics when we have meaningfully bounded the region. We generally do
%not want to even consider anything near maxStepsize initially.
isBracketed=false;

%There are two stages. The approach of Section 3 of [1] only works when the
%modified function value (psi in [1])<=0 and the gradient in the descent
%direction >0. If that is not the case, then we have to use the modified
%approach of Section 2 until that condition is satisfied.
stage1=true;

curIter=0;
while(1)
    %alphaL and alphaU are not ordered in terms of their magnitude, so we
    %set the bounds for the minimum and maxmimum stepsize by ordering them,
    %unless we haven't bracketed the uncertainty interval.
    if(isBracketed)
        alphaMinCur=min(alphaL,alphaU);
        alphaMaxCur=max(alphaL,alphaU);
    else
        %If we haven't bracketed the interval, then we have a minimum step
        %size, but we do not have a maximum step size. We use the heuristic
        %that the maximum should be the test stepsize plus initStepMult
        %times the distance from the lower bound of the stepsize to the
        %test stepsize.
        
        alphaMinCur=alphaMin;
        alphaMaxCur=alphaT+initStepMult*(alphaT-alphaMinCur);
        
        %Check that the bound is not too big.
        alphaMaxCur=min(alphaMaxCur,alphaMax);
        
        %If the upper bound is the actual value
        if(alphaMaxCur==alphaMax)
            isBracketed=true;
        end
    end
    
    if((alphaMaxCur-alphaMinCur)<=RelTolAlpha*alphaMinCur)
        %The uncertainty region for the line search became too small
        %without satisfying all of the necessary conditions.
        exitCode=-996;
        return;
    end
    
    %Compute the value of x after taking the test step.
    xT=x+alphaT*s;
 
    [fT,fGradT]=func(xT);
    
    gT=fGradT'*s;
    
    %A non-finite or imaginary value was encountered.
    if(~isfinite(fT)||~isfinite(gT)||~isreal(fT)||~isreal(gT))
        xL=[];
        fGradL=[];
        fL=[];
        gL=[];
        alphaL=[];

        exitCode=-1001;
        return
    end

    %fTest is the right-hand side of Equation 1.1 in [1].
    fTest=fInit+alphaT*mu*gInit;
    
    curIter=curIter+1;
    
    %The convergence conditions of Equations 1.1 and 1.2 of [1] are
    %the strong Wolfe conditions. Here, we let one specify whether they
    %only want Armijo's rule, they want the regular Wolfe condition, or they
    %also want the strong Wolfe conditions
    switch(conditionType)
        case 0%Just Armijo's rule
            didConverge=(fT<=fTest);
        case 1%Armijo's rule plus the regular Wolfe condition.
            didConverge=(fT<=fTest)&&(gT>=eta*gInit);
        case 2%The strong Wolfe conditions.
            didConverge=(fT<=fTest)&&(abs(gT)<=-eta*gInit);
        otherwise
            error('Unknown conditionType specified).')
    end

    if(didConverge)
        xL=xT;
        fL=fT;
        fGradL=fGradT;
        gL=gT;
        alphaL=alphaT;
        exitCode=curIter;
        return
    end

    if(curIter>=maxIter)
        %If we get here, then the maximum number of iterations was reached.  
        exitCode=-998;
        return
    end

    %If the conditions mentioned in [1] after Equation 3.3 are satisfied,
    %then we can go to the second stage.
    psiVal=fT-fInit-alphaT*mu*gInit;
    if(stage1&&psiVal<=0&&gT>0)
        stage1=false;
        
        %When going from the optimization rule in Section 2 to that in
        %Section 3, the cost function changes from psi to phi. The check
        %below deals with the possibility that changing the functions
        %switched the ordering of fL and fU.

        if(fL>fU)
            temp=fU;
            fU=fL;
            fL=temp;
            
            temp=gU;
            gU=gL;
            gL=temp;
            
            temp=alphaU;
            alphaU=alphaL;
            alphaL=temp;
            
            temp=fGradL; 
            fGradL=fGradU;
            fGradU=temp;

            xL=x+alphaU*s;
        end
    end

    if(stage1)
        %We use the transformed values of Section 2. In other words, we
        %transform the f values into phi values.
        fLMod=fL-alphaL*mu*gInit-fInit;
        gLMod=gL-mu*gInit;

        fTMod=fT-alphaT*mu*gInit-fInit;
        gTMod=gT-mu*gInit;
        
        fUMod=fU-alphaU*mu*gInit-fInit;
        gUMod=gU-mu*gInit;
        
        [alphaT,alphaL,alphaU,xL,fL,fU,gL,gU,fGradL,fGradU,isBracketed]=trialValSelect(xL,xT,fLMod,fTMod,fUMod,gLMod,gTMod,gUMod,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU,alphaMinCur,alphaMaxCur,isBracketed);
        %Undo the transformation
        fL=fL+alphaL*mu*gInit+fInit;
        gL=gL+mu*gInit;
        
        fU=fU+alphaU*mu*gInit+fInit;
        gU=gU+mu*gInit;
    else
        [alphaT,alphaL,alphaU,xL,fL,fU,gL,gU,fGradL,fGradU,isBracketed]=trialValSelect(xL,xT,fL,fT,fU,gL,gT,gU,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU,alphaMinCur,alphaMaxCur,isBracketed);
    end
end

%We should not get here.
end

function [alphaTNew,alphaLNew,alphaUNew,xLNew,fLNew,fUNew,gLNew,gUNew,fGradLNew,fGradUNew,isBracketedNew]=trialValSelect(xL,xT,fL,fT,fU,gL,gT,gU,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU,alphaMin,alphaMax,isBracketed)
%%TRIALVALSELECT This implements the trial value selection of Section 4 of
%           [1] as well as the updating algorithms in Sections 2 and 3. For
%           the algorithm of Section 2, the inputs to this function must be
%           first transformed.
%
%REFERENCES:
%[1] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Update the bounds. Whether or not this is the modified algorithm of
%Section 3 or the algorithm of Section 2 in [1] depends on whether or not
%the inputs to this function were transformed. The updating of the bounds
%is independent on how the new alphaT is chosen and the choice of the new
%alphaT depends on the old values for the bounds. We update the bounds
%before finding alphaT, so we can decide which approach to use to find
%alphaT. Section 2 of [1] suggests performing a bisection step if
%insufficient progress is made on the last two trials. here, we just use
%the last trial, because it is simpler.
if(fT>fL) %Case a
    alphaUNew=alphaT;
    fUNew=fT;
    gUNew=gT;
    fGradUNew=fGradT;
    
    alphaLNew=alphaL;
    xLNew=xL;
    fLNew=fL;
    gLNew=gL;
    fGradLNew=fGradL;
elseif(gT*(alphaL-alphaT)>0)%Case b (fT<=fL)
    alphaUNew=alphaU;
    fUNew=fU;
    gUNew=gU;
    fGradUNew=fGradU;
    
    alphaLNew=alphaT;
    xLNew=xT;
    fLNew=fT;
    gLNew=gT;
    fGradLNew=fGradT;
else%Case c
    alphaUNew=alphaL;
    fUNew=alphaL;
    gUNew=gL;
    fGradUNew=fGradL;
    
    alphaLNew=alphaT;
    xLNew=xT;
    fLNew=fT;
    gLNew=gT;
    fGradLNew=fGradT;
end

isBracketedNew=isBracketed;

if(isBracketed)
    prevWidth=alphaMax-alphaMin;
    width=abs(alphaUNew-alphaLNew);
    
    %0.66 is the ad-hoc value used in Section 2. Basically, if the last
    %step did not make sufficient progress, then the next step will bisect
    %the region. We can only do this if the region is bracketed. Otherwise,
    %we use the step determination method of Section 4.
    if(width>=0.66*prevWidth)
        alphaTNew=(1/2)*(alphaUNew+alphaLNew);
        return;
    end
end

if(fT>fL)%Case 1: The test value is greater than the last best value. This
    %brackets the minimum. As in [1], if the cubic interpolated value is
    %closer to the previous best value, we take it; otherwise, we take the
    %average of the cubic and the quadratic solution.

    %We know the minimum is between fL and fT. As it is not desirable
    %for the interpolation to place the maximum too close too alphaU, we 
    isBracketedNew=true;
    
    alphaC=findCubicFitVertices(alphaL,alphaT,fL,fT,gL,gT);
    alphaQ=findQuadraticFitVertex(alphaL,alphaT,fL,fT,gL);
    if(isempty(alphaC))
        %This check should take care of the eventuality of finite precision
        %errors making it such that no local minimum to the cubic exists.
        %This case is not expected to occur.
        alphaTNew=alphaQ;
    else
        
        alphaC=alphaC(1);%Only take the minimum.

        if(abs(alphaC-alphaL)<abs(alphaQ-alphaL))
            alphaTNew=alphaC;
        else
            alphaTNew=(1/2)*(alphaQ+alphaC);
        end
    end
    
elseif(gT*gL<0)%Case 2
    %The minimum is again bracketed due to the change in the function
    %derivatives, not due to the values.
    isBracketedNew=true;
    
    alphaC=findCubicFitVertices(alphaL,alphaT,fL,fT,gL,gT);
    alphaS=findQuadraticFitVertex(alphaL,alphaT,fL,gL,gT);
    
    if(isempty(alphaC))
        %This check should take care of the eventuality of finite precision
        %errors making it such that no local minimum to the cubic exists.
        %This case is not expected to occur.
        alphaTNew=alphaS;
    else
        alphaC=alphaC(1);%Only take the minimum.

        if(abs(alphaC-alphaT)<abs(alphaS-alphaT))
            alphaTNew=alphaC;
        else
            alphaTNew=alphaS;
        end
    end
elseif(gT*gL>=0&&abs(gT)<=abs(gL))
    %Case 3: the test value fT is less than the current lower bound fL, the
    %derivatives have the same sign, and 
    alphaC=findCubicFitVertices(alphaL,alphaT,fL,fT,gL,gT);
    alphaS=findQuadraticFitVertex(alphaL,alphaT,gL,gT);

    %If alphaC does not exist, is infinite, or alphaC is in the wrong
    %direction 
    if(isempty(alphaC)||~isfinite(alphaC(1))||...%if alphaC does not exist
       ((alphaT>alphaL)&&(alphaC(1)<alphaT))||...%If alphaC is in the wrong direction.
       ((alphaT<alphaL)&&(alphaC(1)>alphaT)))%If alphaC is in the wrong direction.
        %The secant step should always exist and be in the correct
        %direction.
        alphaTNew=alphaS;
    else
        alphaC=alphaC(1);%Take the minimum.
        if(abs(alphaC-alphaT)<abs(alphaS-alphaT))
            alphaTNew=alphaC;
        else
            alphaTNew=alphaS;
        end
    end
    
    %Make sure that the alphaT value is not too close to the existing upper
    %bound point, but this only applies if the solution has been bracketed.
    if(isBracketed)
        delta=0.66;
        if(alphaT>alphaL)
            alphaTNew=min(alphaT+delta*(alphaU-alphaT),alphaTNew);
        else
            alphaTNew=max(alphaT+delta*(alphaU-alphaT),alphaTNew);
        end
    end
else%Case 4
    %If function is not already bracketed, then  the step is either the
    %minimum or the maximum step.
    if(isBracketed)
        alphaC=findCubicFitVertices(alphaT,alphaU,fT,fU,gT,gU);
        
        %If there is no solution to the cubic, then bisect the region.
        if(isempty(alphaC))
            alphaTNew=(1/2)*(alphaUNew+alphaLNew);
            return
        end

        alphaC=alphaC(1);%Only take the minimum.
        alphaTNew=alphaC;
    else
        if(alphaT>alphaL)
            alphaTNew=alphaMax;
        else
            alphaTNew=alphaMin;
        end
    end
end

%Clip the new alphaT to the bounds of the region.
alphaTNew=min(max(alphaTNew,alphaMin),alphaMax);

if(alphaTNew==alphaLNew||isBracketed&&(alphaTNew==alphaUNew))
    %Make sure that alphaT is not on the lower bound and if not previously
    %bracketed, not on the upper bound. This is a change from [1].
    alphaTNew=(1/2)*(alphaUNew-alphaLNew);
end
end

function [xL,fL,fGradL,gL,alphaL,exitCode]=lineSearchBisection(x,func,fInit,fGradInit,s,maxStepSize,RelTolAlpha,fTol,wolfeTol,conditionType,maxIter,minStepSize)
%%LINESEARCHBISECTION This function implements the Mor�-Thuente algorithm
%               of [1], but assumes that stepSize is an upper bound that is
%               provided and only performs bisection of the uncertainty
%               region, rather than the interpolation used in [1].
%
%REFERENCES:
%[1] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%We use the notation in the paper.
mu=fTol;
eta=wolfeTol;

%The initial bounds are the function value evaluated at the zero step and
%the initial test step. Here is everything at the intial point.
alphaL=0;
fL=fInit;
fGradL=fGradInit;
gL=fGradInit'*s;%The initial gradient in the search direction.
xL=x;

%The initial gradient is needed for the right-hand side of Equation 1.1 in
%[1].
gInit=gL;

%If the search direction is not a descent direction.
if(gInit>0)
    xL=[];
    fL=[];
    fGradL=[];
    gL=[];
    alphaL=[];
    
    exitCode=-994;
    return;
elseif(gInit==0)%If the  function is already at a vertex.
    %The function could be stuck at a point of inflection.
    exitCode=0;
    return;
end

%If the initial step size is invalid.
if(maxStepSize<=0)
    xL=[];
    fL=[];
    fGradL=[];
    gL=[];
    alphaL=[];
    
    exitCode=-995;
    return
end

alphaU=maxStepSize;
xU=x+maxStepSize*s;
[fU,fGradU]=func(xU);
gU=fGradU'*s;

%A non-finite or imaginary value was encountered.
if(~isfinite(fU)||~isfinite(gU)||~isreal(fU)||~isreal(gU))
    xL=[];
    fL=[];
    fGradL=[];
    gL=[];
    alphaL=[];
    exitCode=-1001;
    return
end

%Before starting, check whether the initial test step is sufficient.
%fTest is the right-hand side of Equation 1.1 in [1].
fTest=fInit+alphaU*mu*gInit;
switch(conditionType)
    case 0%Just Armijo's rule
        didConverge=(fU<=fTest);
    case 1%Armijo's rule plus the regular Wolfe condition.
        didConverge=(fU<=fTest)&&(gU>=eta*gInit);
    case 2%The strong Wolfe conditions.
        didConverge=(fU<=fTest)&&(abs(gU)<=-eta*gInit);
    otherwise
        error('Unknown conditionType specified).')
end

%The convergence conditions of Equations 1.1 and 1.2 of [1] are
%the strong Wolfe conditions. Here, we let one specify whether they
%only want Armijo's rule, they want the regular Wolfe condition, or they
%also want the strong Wolfe conditions.
if(didConverge)
    xL=xU;
    fL=fU;
    fGradL=fGradU;
    gL=gU;
    alphaL=alphaU;
    exitCode=0;
    return
end

%We start in stage 1 of the algorithm. Thus, we order upper and lower
%bounds based on the transformed function in Section 2.
fLMod=fL-alphaL*mu*gInit-fInit;
fUMod=fU-alphaU*mu*gInit-fInit;
if(fLMod>fUMod)
    temp=fU;
    fU=fL;
    fL=temp;

    temp=gU;
    gU=gL;
    gL=temp;

    temp=alphaU;
    alphaU=alphaL;
    alphaL=temp;
    
    temp=fGradL;
    fGradL=fGradU;
    fGradU=temp;

    %xU is not needed after this.
    xL=xU;
end

%The test step bisects the region.
alphaT=(1/2)*(alphaL+alphaU);

%There are two stages. The approach of Section 3 of [1] only works when the
%modified function value (psi in [1])<=0 and the gradient in the descent
%direction >0. If that is not the case, then we have to use the modified
%approach of Section 2 until that condition is satisfied.
stage1=true;

curIter=1;
while(1)
    %alphaL and alphaU are not ordered in terms of their magnitude, so we
    %set the bounds for the minimum and maxmimum stepsize by ordering them,
    %unless we haven't bracketed the uncertainty interval.
    alphaMinCur=min(alphaL,alphaU);
    alphaMaxCur=max(alphaL,alphaU);
    
    if((alphaMaxCur-alphaMinCur)<=RelTolAlpha*alphaMinCur)
        %The uncertainty region for the line search became too small
        %without satisfying all of the necessary conditions.
        exitCode=-996;
        return;
    end
    
    %Compute the value of x after taking the test step.
    xT=x+alphaT*s;
 
    [fT,fGradT]=func(xT);
    
    gT=fGradT'*s;
    
    %A non-finite or imaginary value was encountered.
    if(~isfinite(fT)||~isfinite(gT)||~isreal(fT)||~isreal(gT))
        exitCode=-1001;
        return
    end

    %fTest is the right-hand side of Equation 1.1 in [1].
    fTest=fInit+alphaT*mu*gInit;
    
    curIter=curIter+1;
    
    %Check for convergence.
    switch(conditionType)
        case 0%Just Armijo's rule
            didConverge=(fT<=fTest);
        case 1%Armijo's rule plus the regular Wolfe condition.
            didConverge=(fT<=fTest)&&(gT>=eta*gInit);
        case 2%The strong Wolfe conditions.
            didConverge=(fT<=fTest)&&(abs(gT)<=-eta*gInit);
        otherwise
            error('Unknown condition specified).')
    end

    if(didConverge)
        xL=xT;
        fL=fT;
        gL=gT;
        alphaL=alphaT;
        exitCode=curIter;
        return
    end

    if(curIter>=maxIter)
        %If we get here, then the maximum number of iterations was reached.  
        exitCode=-998;
        return
    end
    
    %If alphaT was clipped to the minimum step size.
    if(alphaT<=minStepSize)
        if(fT<fL)
            xL=xT;
            fL=fT;
            gL=gT;
            alphaL=alphaT;
        end

        exitCode=-1000;
        return;
    end
    
    %If the conditions mentioned in [1] after Equation 3.3 are satisfied,
    %then we can go to the second stage.
    psiVal=fT-fInit-alphaT*mu*gInit;
    if(stage1&&psiVal<=0&&gT>0)
        stage1=false;
        
        %When going from the optimization rule in Section 2 to that in
        %Section 3, the cost function changes from psi to phi. The check
        %below deals with the possibility that changing the functions
        %switched the ordering of fL and fU.
        if(fL>fU)
            temp=fU;
            fU=fL;
            fL=temp;
            
            temp=gU;
            gU=gL;
            gL=temp;
            
            temp=alphaU;
            alphaU=alphaL;
            alphaL=temp;
            
            temp=fGradL; 
            fGradL=fGradU;
            fGradU=temp;

            xL=x+alphaU*s;
        end
    end
    
    if(stage1)
        %We use the transformed values of Section 2. In other words, we
        %transform the f values into phi values.
        fLMod=fL-alphaL*mu*gInit-fInit;
        gLMod=gL-mu*gInit;

        fTMod=fT-alphaT*mu*gInit-fInit;
        gTMod=gT-mu*gInit;
        
        fUMod=fU-alphaU*mu*gInit-fInit;
        gUMod=gU-mu*gInit;
        
        [alphaT,alphaL,alphaU,xL,fL,fU,gL,gU,fGradL,fGradU]=trialValSelBisect(xL,xT,fLMod,fTMod,fUMod,gLMod,gTMod,gUMod,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU);
        %Undo the transformation
        fL=fL+alphaL*mu*gInit+fInit;
        gL=gL+mu*gInit;
        
        fU=fU+alphaU*mu*gInit+fInit;
        gU=gU+mu*gInit;
    else
        [alphaT,alphaL,alphaU,xL,fL,fU,gL,gU,fGradL,fGradU]=trialValSelBisect(xL,xT,fL,fT,fU,gL,gT,gU,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU);
    end
    
    alphaT=max(alphaT,minStepSize);
end

%We should not get here.
end

function [alphaTNew,alphaLNew,alphaUNew,xLNew,fLNew,fUNew,gLNew,gUNew,fGradLNew,fGradUNew]=trialValSelBisect(xL,xT,fL,fT,fU,gL,gT,gU,fGradL,fGradT,fGradU,alphaL,alphaT,alphaU)
%%%TRIALVALSELBISECT This implements the updating algorithms in Sections 2
%               and 3 of [1], where the difference between the algorithms
%               depends on how one transforms the values before passing
%               them to this function. Rather than using the rules of
%               Section 4 of [1] to find the stepsize, this function just
%               uses bisection to get a new test point alphaTNew.
%
%REFERENCES:
%[1] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Update the bounds. Whether or not this is the modified algorithm of
%Section 3 or the algorithm of Section 2 in [1] depends on whether or not
%the inputs to this function were transformed.
if(fT>fL) %Case a
    alphaUNew=alphaT;
    fUNew=fT;
    gUNew=gT;
    fGradUNew=fGradT;
    
    alphaLNew=alphaL;
    xLNew=xL;
    fLNew=fL;
    gLNew=gL;
    fGradLNew=fGradL;
elseif(gT*(alphaL-alphaT)>0)%Case b (fT<=fL)
    alphaUNew=alphaU;
    fUNew=fU;
    gUNew=gU;
    fGradUNew=fGradU;
    
    alphaLNew=alphaT;
    xLNew=xT;
    fLNew=fT;
    gLNew=gT;
    fGradLNew=fGradL;
else%Case c
    alphaUNew=alphaL;
    fUNew=alphaL;
    gUNew=gL;
    fGradUNew=fGradL;
    
    alphaLNew=alphaT;
    xLNew=xT;
    fLNew=fT;
    gLNew=gT;
    fGradLNew=fGradT;
end

%Just bisect the new region.
alphaTNew=(1/2)*(alphaUNew+alphaLNew);

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
