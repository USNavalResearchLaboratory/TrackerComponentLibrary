function [xEst,exitCode,path]=homotopyPathFollow(x0,f,df,homF,dHomF,errorTols,stepSizeParams,lambdaTurnAround,iterParams)
%%HOMOTOPYPATHFOLLOW Follow a homotopy path to find a single real zero of a
%          continuous, differentiable multivariate real function f(x). The
%          algorithm is based in the dense normal flow algorithm described
%          in [1] and [2]. However, there are notable modifications,
%          described below, that allow the solver to find solutions where
%          the original algorithms would fail. For example, finding the
%          zero in f=@(x)x^2 would fail in [1] and [2], but works here.
%          Details are described below. Theoretical convergence conditions
%          are given in [3].
%
%INPUTS: x0 The xDimX1 initial estimate. This does not have to be very
%           accurate, but algorithmic success is more likely when it is not
%           orders of mangitude away from the true value.
%         f A handle to the function to be zeroed. f(x) takes an xDImX1
%           input and returns and xDimX1 output.
%        df A handle to a function that returns the xDimXxDim Jacobian of
%           f. The value in row i and column j of the output is the
%           derivative of the ith element of f(x) with respect to the jth
%           element of x.
% homF,dHomF Optional parameters allowing the user to specify the homotopy
%           to be used. If omitted or empty matrices are passed, then the
%           homotopy 0=f(x)+(lambda-1)*f(x0) is used. This works well with
%           a wide variety of functions, but is not necessarily of the form
%           to guarantee a solution for arbitrary inputs to any particular
%           function. If provided, these have the input formats
%           homF(y,fx,x0,f0) and dHomF(y,fx,dfx,x0,f0), where y=[lambda;x]
%           is a vector of the homotopy parameter lambda stacked on the
%           vector to be estimated x. fx and dfx are the xDimX1 and
%           xDimXxDim values of f(x) and df(x) with f and df as defined
%           above. x0 is the xDImX1 vector x0 as poaased above and
%           f0=f(x0).
% errorTols The algorithm has a number of error tolerances that can be
%           adjusted to adjust the performance. If one wishes to use
%           something other than the default parameters, then a structure
%           whose members are the parameters that should be changed can be
%           passed. The possible members and their default values are:
%           AbsTolArc, RelTolArc these are absolute and relative tolerances
%                      for how well the homotopy path is tracked. These
%                      adjust the accuracy of the Newton refinement step of
%                      the Euler prediction. The scale used for the
%                      relative error is scaleX=max(abs(x)) for the
%                      cufrrent estimate of x. These are tolerances of the
%                      change in the gradient. The default values are both
%                      1e-10.
%           AbsTolAns,RelTolAns These are tolerances related to the
%                      gradients in the stage of the estimator. The endgame
%                      is after the solution has been bracketed (a point
%                      with lambda<=1 and a point with lambda>=1 have been
%                      obtained, or a turnaround in the curve close to
%                      lambda=1 is detected). The default values are both
%                      1e-10.
%              endTolF After the bracketed endgame finished, it often
%                      occurs that despite convergence, f(xEst) is not all
%                      that close to 0. If any element are larger than
%                      endTolF, then xTolEst will be passed to Newtons's
%                      algorithm to iterate on using df and f without line
%                      search to refine the estimates. The default is 1e-7.
% stepSizeParams These are a number of parameters that affect the maximum
%           and minimum stepsizes. If one wishes to use something other
%           than the default parameters, then a structure whose members are
%           the parameters that should be changed can be passed. The
%           possible members and their default values are:
%           deltaLambdaMin, deltaFMin These are the minimum stepsizes
%                      allowed in terms of the change in the homotopy
%                      parameter lambda and in terms of the approximate
%                      fraction of the distance traveled to zero the
%                      function. This values are explained in more detail
%                      below. Default values are eps() and 1e-6.
%           deltaLambdaMax,deltaFMax These are the minimum stepsizes
%                      allowed in terms of the change in the homotopy
%                      parameter lambda and in terms of the approximate
%                      fraction of the distance traveled to zero the
%                      function. Note that deltaFMax should be >=1. and
%                      that lambda only does to 1, so deltaLambdaMax<=1.
%                      This values are explained in more detail below.
%                      Default values are 0.08 and 1.1.
%           LBar, RBar, DBar, q These parameter relate to the optimal
%                      stepsize prediction based on curbature criteria and
%                      are described in [1]. They are respectively the
%                      ideal values of the contraction factor, the residual
%                      factor, the distance factor and also the assumed
%                      operating order. Default values are respectively
%                      0.5, 0.01, 0.5, and 2.
%           BMin, BMax These are the minimum stepsize increase factor and
%                      the maximum stepsize increase factor. Default values
%                      are 1 and 10.
% lambdaTurnAround The value of the homotopy parameter lambda from which it
%           is assumed that a turnaround in the path means that the path
%           grazed the lambda=1 point without crossing. If a tournaround in
%           the path is detected after this point, then an endgame routine
%           for turnarounds is initiated, as described below. Such a
%           routine is not in the algorithm in [1] and [2]. The default
%           value if omitted or an empty matrix is passed is 0.9.
% iterParams There are a number of parameters related to iterations of
%           different parts of the algorithm. If one wishes to use
%           something other than the default parameters, then a structure
%           whose members are the parameters that should be changed can be
%           passed. The possible members and their default values are:
%           maxIter The maxinum number of steps to take along the homotpy
%                   path. The dault is 1024.
%           maxNewtonIter The maximum number of iterations to perform on
%                   the Newton corrector step of the algorithm without
%                   observaing convergence before declaring a stepsize too
%                   large and reducting it. The default if omitted or an
%                   empty matrix is passed is 4. This must be >=2.
%           maxRefineIter If a solution has been bracketed, steps with
%                   lambda<1 and lambda>=1 have been encountered, then this
%                   is the number of iterations of the bracketed endgame
%                   refinement to perform before declaring a failure to
%                   converge. The default is 
%                   2*(ceil(abs(log10(AbsTolAns+RelTolAns)))+1);
%                   maxIterNewton=25;
%           maxIterNewton When in the endgame when a turnaround has been
%                   detected, or if in a bracketed endgame, the precision
%                   requested by endTolF has not been attained, then a few
%                   iterations of Newtons method on f will be executed.
%                   This is the maximum number of iterations to perform
%                   before declaring a failure to converge. The default is
%                   25.
%
%OUTPUTS: xEst The xDimX1 estimate of the value making f(x)=0, or an empty
%              matrix if the algorithm failed. The algorithm only finds one
%              solution; there could be more than 1 solution.
%     exitCode A value indicating how the algorithm terminated. Possible
%              values are:
%              0 The algorithm finished successfully.
%              1 Convergence not attained after maxIter steps along the
%                homotopy path.
%              2 Unable to find an acceptable stepsize. The stepsize become
%                too small.
%              3 A non-finite value was encountered.
%              4 Unable to interpolate the turnaround point in the extremem
%                endgame.
%              5 Unable to interpolate the solution in the bracketed
%                endgame.
%              6 The bracketed endgame did not converge within
%                maxRefineIter iterations.
%         path This is the homotopy path that the algorithm followed until
%              termination.
%
%The basic idea being homotopy algorithms is that one has a homotopy
%h(x,lambda)=0 at x0 and lambda=0 and x0 and is formulated so that it
%equals a desired value such that f(x)=0 at lambda=1. One then writes a
%differential equation in terms of lambda and integrates from 0 to 1. In
%probability-1 homotopy algorithms, one uses an augmented state
%y=[lambda;x] and writes a differential equation in terms of another
%parameter s, on which lambda and x are assumed to depend. One then starts
%at y=[0;x0] and integrates however long along s it takes until y=[1;x],
%where x will be the solution. This function implements the probability-1
%homotopy algorithm described in [1] and [2] with a number of
%modifications.
%
%First, the algorithm in [1] and [2] is not capable of solving trivial
%problems such as f(x)=x^2, because the homotopy path in s only touches
%lambda=1 but never crosses it. To fix this, an extra endgame solver was
%added. If lambda at the current or previous step is >=lambdaTurnAround,
%and by observing the sign of the homotopy path direction in lambda at the
%current and previous steps, it is determined that the path has turned
%around, then it is assumed that the path grazed lambda=1 in between the
%steps. Two point Hermite interpolation is used to determine the extremum
%point and then the estimate of x is refined using Newtons method on f, the
%assumption being that the estimate is sufficiently close to f(x)=0 to
%converge.
%
%Another change of the algorithm from [1] and [2] is in the determination
%of the minimum and maximum stepsizes alogn the path s. The problem is that
%the scale of problem is unknown a prioi and changes based on how bad the
%estimates are. An initial/maximum step size of 1, as used in [1] might be
%far too small for the algorithm too converge, or it could be far too
%large. One measure of scale if how far one might have to travel in s to
%make a change in lambda. Since lambda can only go from 0 to 1, limiting
%the change in lambda is more scale invariant. If yp is the direction along
%the homotopy that would be traveled (the tangent unit vector), then
%lambdaMax/yp(1) is a linear approximation to how much s can change to take
%a maximum sized step. Of course, the path can turn around, so this is not
%a meaningful step size if yp(1)=0. Thus, we also se a limit to how far the
%algorithm can step in terms of zeroing f(x). A linear approximation to
%zero f can be obtained by solving  f-dfds*DeltaS=0 for DeltaS, where
%DeltaS is the stepsize and dfds is the derivative of f with respect to s.
%dfds=df*yp(2:end). The issue is that not all elements of f will be zeroed
%with the same stepsize when performing a linear approximation, so we use 
%DeltaS=pinv(df*yp(2:end))*f. By limiting the stepsize to DeltaS being some
%value over 1, we allow the algorithm to converge byt prevent overly large
%stepsizes.
%
%The final changes from [1] and [2] are that the definitions of absolute
%and relative tolerance for declaring convergence were slightly change from
%[1] and that the default homotopy is one that is suggested in neither [1]
%nor [2]. Better results on generic functions were observed with the
%homotopy 0=f(x)+(lambda-1)*f(x0).
%
%EXAMPLE 1:
%Here we solve a trivial problem where the algorithm as written in [1[
%would fail. We give it a terrible initial estimate:
% f=@(x)(x-1)^2;
% df=@(x)2*(x-1);
% x0=8000;
% [xEst,exitCode]=homotopyPathFollow(x0,f,df)
%One will get xEst very close to 1, as expected.
%
%EXAMPLE 2:
%Here, we show how we can solve a 8-dimensional second order set of
%equations, but since not all conditions of [3] are satisfied to guarantee
%global convergence, some initial estimatesresult in periodic homotopies
%that do not converge. Randomly trying a few initial estimates generally
%gets a solution. This example is taken from [4]
% f=@(x)[x(1)^2+x(2)^2-1;
%        x(3)^2+x(4)^2-1;
%        x(5)^2+x(6)^2-1;
%        x(7)^2+x(8)^2-1;
%        0.004731*x(1)*x(3)-0.3578*x(2)*x(3)-0.1238*x(1)-0.001637*x(2)-0.9338*x(4)+x(7)-0.3571;
%        0.2238*x(1)*x(3)+0.7623*x(2)*x(3)+0.2638*x(1)-0.07745*x(2)-0.6734*x(4)-0.6022;
%        x(6)*x(8)+0.3578*x(1)+0.004731*x(2);
%        -0.7623*x(1)+0.2238*x(2)+0.3461];
% df=@(x)[        2*x(1),                 2*x(2),                         0,       0,      0,      0, 0,      0;
%                      0,                      0,                    2*x(3),  2*x(4),      0,      0, 0,      0;
%                      0,                      0,                         0,       0, 2*x(5), 2*x(6), 0,      0;
%                      0,                      0,                         0,       0,      0,      0, 2*x(7), 2*x(8);
%  -0.1238+0.004731*x(3),  -0.001637-0.3578*x(3), 0.004731*x(1)-0.3578*x(2), -0.9338,      0,      0, 1,      0;
%     0.2638+0.2238*x(3),   -0.07745+0.7623*x(3),   0.2238*x(1)+0.7623*x(2), -0.6734,      0,      0, 0,      0;
%                 0.3578,               0.004731,                         0,      0,       0,   x(8), 0,    x(6);
%                -0.7623,                 0.2238,                         0,      0,       0,      0, 0,      0];
% x0=ones(8,1);
% [xEst,exitCode,path]=homotopyPathFollow(x0,f,df);
% xEst
% f(xEst)%One can see this is zero within precision bounds.
% exitCode
% figure(1)
% clf
% hold on
% for idx=2:9
%   plot(path(1,:),path(idx,:))
% end
% %This converges to a solution and the path goes there rather directly.
% %The horizontal axis of the plot is lambda and the vertical axis is all
% %of the other values. However, the following initialization provies a
% %period result, as one can see plotting the path.
% x0=[1.5442;0.0859;-1.4916;-0.7423;-1.0616;2.3505;-0.6156;0.7481];
% [xEst,exitCode,path]=homotopyPathFollow(x0,f,df);
% xEst
% exitCode
% figure(2)
% clf
% hold on
% for idx=2:9
% plot(path(1,:),path(idx,:))
% end
%
%REFERENCES:
%[1] L. T. Watson, S. C. Billups, and A. P. Morgan, "ALGORITHM 652 HOMPACK:
%    A suite of codes for globally convergent homotopy algorithms," ACM
%    Transactions on Mathematical Software, vol. 13, no. 3, Sep. 1987, pp.
%    281-310.
%[2] L. T. Watson,  M. Sosonkina, R. C. Melville, A. P. Morgan, and H. F.
%    Walker, "Algorithm 777: HOMPACK90: A Suite of Fortran 90 Codes for
%    Globally Convergent Homotopy Algorithms," ACM transactions on
%    Mathematical Software, vol. 23, no. 4, Dec. 1997, pp. 514-449.
%[3] L. T. Watson, "Probability-One Homotopies in Computational Science,"
%    Journal of Computational and Applied Mathematics, vol. 140, no. 1-2.
%    Mar. 2002, pp. 785-807.
%[4] A. Morgan and V. Shapiro. "Box-bisection for solving second-degree
%    systems and the problem of clustering," ACM Transactions on
%    Mathematical Software, vol. 13 no. 2, Jun. 1987, pp, 152-167.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

f0=f(x0);
if(nargin<4||isempty(homF))
    if(nargin>4&&isempty(dHomF))
        error('homF and dHomF must be given together if provided.')
    end
    
    %The default homotopy
    H=@(y,fy)homotopyF(y,fy,x0,f0);
    DH=@(y,fy,dfy)homotopyFJacob(y,fy,dfy,x0,f0);
elseif(nargin>=4&&isempty(dHomF))
    error('homF and dHomF must be given together if provided.')
else
    H=@(y,fy)homF(y,fy,x0,f0);
    DH=@(y,fy,dfy)dHomF(y,fy,dfy,x0,f0);
end

%Default error tolerances.
%Absolute an relative tolerances for tracking the homotopy.
AbsTolArc=1e-10;
RelTolArc=1e-10;
%Absolute and relative tolerances in the "endgame" for obtaining the final
%answer.
AbsTolAns=1e-10;
RelTolAns=1e-10;
%Tolerance on f when it is not close enough to zero in the end.
endTolF=1e-7;

if(nargin>5&&~isempty(errorTols))
    if(isfield(errorTols,'AbsTolArc'))
       AbsTolArc=errorTols.AbsTolArc;
    end
    if(isfield(errorTols,'RelTolArc'))
        RelTolArc=errorTols.RelTolArc;
    end
    if(isfield(errorTols,'AbsTolAns'))
        AbsTolAns=errorTols.AbsTolAns;
    end
    if(isfield(errorTols,'RelTolAns'))
        RelTolAns=errorTols.RelTolAns;
    end
    if(isfield(errorTols,'endTolF'))
        endTolF=errorTols.endTolF;
    end
end

%Default stepsize control parameters.
%This is used to determine hMin.
deltaLambdaMin=eps();
%This is used to determine hMax.
deltaLambdaMax=0.08;

deltaFMin=1e-6;
deltaFMax=1.1;

LBar=0.5;%Ideal contraction factor.
RBar=0.01;%Ideal residual factor.
DBar=0.5;%Ideal distance factor.
q=2;%Assumed operating order.

BMin=1;%Minimum step size increase factor.
BMax=10;%Maximum step size increase factor.

if(nargin>6&&~isempty(stepSizeParams))
    if(isfield(stepSizeParams,'deltaLambdaMin'))
        deltaLambdaMin=stepSizeParams.deltaLambdaMin;
    end
    if(isfield(stepSizeParams,'deltaLambdaMax'))
        deltaLambdaMax=stepSizeParams.deltaLambdaMax;
    end
    if(isfield(stepSizeParams,'deltaFMin'))
        deltaFMin=stepSizeParams.deltaFMin;
    end
    if(isfield(stepSizeParams,'deltaFMax'))
        deltaFMax=stepSizeParams.deltaFMax;
    end
    if(isfield(stepSizeParams,'LBar'))
        LBar=stepSizeParams.LBar;
    end
    if(isfield(stepSizeParams,'RBar'))
        RBar=stepSizeParams.RBar;
    end
    if(isfield(stepSizeParams,'DBar'))
        DBar=stepSizeParams.DBar;
    end
    if(isfield(stepSizeParams,'q'))
        q=stepSizeParams.q;
    end
    if(isfield(stepSizeParams,'BMin'))
        BMin=stepSizeParams.BMin;
    end
    if(isfield(stepSizeParams,'BMax'))
        BMax=stepSizeParams.BMax;
    end
end

%Above this value of lambda, a turn-around in the path is assumed to be
%indicative of a solution that touches but does not cross lambda=1.
if(nargin<8||isempty(lambdaTurnAround))
    lambdaTurnAround=0.9;
end

maxIter=1024;%Maximum overall iterations.
maxNewtonIter=4;%This must be >=2.
maxRefineIter=2*(ceil(abs(log10(AbsTolAns+RelTolAns)))+1);
maxIterNewton=25;

%Parameters related to iterations.
if(nargin>8&&~isempty(iterParams))
    if(isfield(iterParams,'maxIter'))
        maxIter=iterParams.maxIter;
    end
    
    if(isfield(iterParams,'maxNewtonIter'))
        maxNewtonIter=iterParams.maxNewtonIter;
    end

    if(isfield(iterParams,'maxRefineIter'))
        maxRefineIter=iterParams.maxRefineIter;
    end

    if(isfield(iterParams,'maxIterNewton'))
        maxIterNewton=iterParams.maxIterNewton;
    end
end

xDim=size(x0,1);

y=[0;x0];%The time-augmented state.
%This previous tangent ensures that the first step goes in the direction of
%increasing lambda.
ypPrev=[1;zeros(xDim,1)];
yPrev=zeros(xDim+1,1);
df0=df(x0);
yp=getTangent(y,ypPrev,DH,f0,df0);

if(yp(1)==0)
    %We do not know the direction of increasing lambda (since the
    %derivative with respect to lambda here is 0). There is thus +/-
    %ambiguity regarding the vector returned by getTangent. Thus, we will
    %determine the direction that is most in the direction decreasing the
    %values of the components. These will correspond to positive values of
    %deltaS needed to minimize the components.

    DeltaSVals=-f0./(df0*yp(2:end));
    if(sum(DeltaSVals)<0)
        yp=-yp;
    end
end

%Initial stepsize heuristic.
DeltaS=abs(pinv(df0*yp(2:end))*f0);
hMax=min(0.1*DeltaS,deltaLambdaMax/abs(yp(1)));
h=hMax;
path=zeros(xDim+1,maxIter+1);
path(:,1)=y;
for k=1:maxIter
    if((y(1)>=lambdaTurnAround||yPrev(1)>=lambdaTurnAround)&&sign(yp(1))*sign(ypPrev(1))<=0)
        [xEst,exitCode]=extremumEndgame(y,yp,yPrev,ypPrev,f,df,maxIterNewton,AbsTolAns,RelTolAns);
        path=path(:,1:k);%Size to fit.
        return;
    end

    firstTry=true;
    while(1)
        %An explicit Euler step:
        z0=y+h*yp;

        NewtonConverged=false;
        zN=z0;
        fzN=f(zN(2:end));%Function value
        dfzN=df(zN(2:end));%Function Jacobian
        HzN=H(zN,fzN);%Homotopy value
        normHz0=norm(HzN);%For use in Equation 22
        DHzN=DH(zN,fzN,dfzN);%Homotopy Jacobian.
        NewtonIter=1;
        while(1)
            if(any(~isfinite(DHzN(:)))||any(~isfinite(HzN)))
                %A non-finite value was encountered.
                xEst=[];
                exitCode=3;
                path=path(:,1:(k+1));%Size to fit.
                return
            end
            
            DeltaZ=-pinv(DHzN)*HzN;
            zN=zN+DeltaZ;

            %The values saved here are used in stepsize selection.
            if(NewtonIter==1)
                z1=zN;
                z2=z1;%So L=0 if convergence occurs after 1 iteration.
                fzN=f(zN(2:end));%Function value
                HzN=H(zN,fzN);%Homotopy value
                normHz1=norm(HzN);
            elseif(NewtonIter==2)
                z2=zN;
            end
            
            diffX=max(abs(DeltaZ(2:end)));
            scaleX=max(abs(zN(2:end)));
            if(abs(DeltaZ(1))<=AbsTolArc&&(diffX<=AbsTolArc||diffX<=RelTolArc*scaleX))
                NewtonConverged=true;
                break;
            end
            if(NewtonIter==maxNewtonIter)
                break;
            end
            if(NewtonIter>1)
                %It was already computed for the first iteration.
                fzN=f(zN(2:end));%Function value
                HzN=H(zN,fzN);%Homotopy Valu
            end
            dfzN=df(zN(2:end));%Function Jacobian
            DHzN=DH(zN,fzN,dfzN);
            NewtonIter=NewtonIter+1;
        end

        if(NewtonConverged)
            break
        end

        hOld=h;
        h=h/2;
        DeltaS=abs(pinv(dfzN*yp(2:end))*fzN);
        hMin=max(eps(h),min(deltaFMin*DeltaS,deltaLambdaMin/abs(yp(1))));
        if(h<hMin)
            xEst=[];
            exitCode=2;
            path=path(:,1:(k+1));%Size to fit.
            return;
        end
        firstTry=false;
    end

    yPrev=y;
    ypPrev=yp;
    y=zN;
    x=y(2:end);
    fy=f(x);%Function value
    dfy=df(x);%Function Jacobian
    yp=getTangent(y,ypPrev,DH,fy,dfy);

    path(:,k+1)=y;
    if(y(1)<1)
        L=norm(z2-z1)/norm(z1-z0);%Equation 21.
        R=normHz1/normHz0;%Equation 22.
        D=norm(z1-y)/norm(z0-y);%Equation 23.

        hHat=min([LBar/L,RBar/R,DBar/D])^(1/q)*h;
        
        DeltaS=abs(pinv(dfy*yp(2:end))*fy);
        
        %Compute for the given stepsize.
        hMin=max(h,min(1e-3*DeltaS,deltaLambdaMin/abs(yp(1))));
        %Here, the maximum allowed value of h depends on how large the
        %current step would have to be to exceed the maximum change in
        %lambda.
        hMax=min(deltaFMax*DeltaS,deltaLambdaMax/abs(yp(1)));
        
        hBar=min([max([hMin,BMin*h,hHat]),BMax*h,hMax]);

        if(firstTry==false)
            hBar=min([hOld,hBar]);%Equation 28
        end
        if(NewtonIter==1)
            h=max(h,hBar);%Equation 27
        elseif(NewtonIter==maxNewtonIter)
            h=min(h,hBar);%Equation 29
        else
            h=hBar;
        end
        continue;
    else
        [xEst,exitCode]=bracketedEndgame(y,yp,yPrev,ypPrev,H,DH,f,df,endTolF,maxRefineIter,maxIterNewton,AbsTolAns,RelTolAns);
        path=path(:,1:(k+1));%Size to fit.
        return;
    end
end
%If we get here, then it did not converge.
xEst=[];
exitCode=1;
path=path(:,1:(k+1));%Size to fit.
end

function yp=getTangent(y,ypPrev,DH,f,df)
%%GETTANGENT Get the tangent vector to the homotopy. DH returns the
%       gradient of the homotopy (a rectangular matrix). The tangent is any
%       vector orthogonal to the gradient. however, we specifically choose
%       one so as to maintain the most continuity in the direction traveled
%       along the homotopy path (the inner product of this vector and the
%       previously found tangent is >=0, meaning that they are no more than
%       90 degrees apart).
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %Get the tangent direction.
    yp=null(DH(y,f,df));

    %In the event of a bifurcation, just take the first one.
    yp=yp(:,1);

    %yp should be unit length; we force it to be so.
    yp=yp/norm(yp);
    %Ensure continuity along the curve. --i.e. get the sign correct so that
    %it does not backtrack on itself.
    if(yp'*ypPrev<0)
        yp=-yp;
    end
end

function val=homotopyF(y,fy,x0,f0)
%%HOMOTOPYF The function for the homotopy 0=f(x)+(lambda-1)*f(x0).
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

lambda=y(1);

val=fy+(lambda-1)*f0;
end

function val=homotopyFJacob(y,fy,dfy,x0,f0)
%%HOMOTOPYFJACOB The Jacobian for the homotopy 0=f(x)+(lambda-1)*f(x0).
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=[f0,dfy];
end

function [xEst,exitCode]=extremumEndgame(y,yp,yPrev,ypPrev,f,df,maxIterNewton,AbsTolAns,RelTolAns)
%%EXTREMUMENDGAME When a turnaround sufficiently close to lambda=1 has been
%          detected, it can often be assumed that the homotopy path has
%          touched lambda=1. This function starts on that assumption and
%          searches for the point where it touches the find the solution.
%          The algorithm uses interpolation and then iterations of Newtons
%          method.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    yDim=size(y,1);
    %If the homotopy path turned around.
    %Interpolate the extemum.
    x=[yPrev,y;
       ypPrev,yp];
    [a,c]=HermiteInterpMultiPoly([0,1],x,yDim);
    aPower=convertPolynomialForm('Newton','PowerSeries',a(1,:),c(1,:));
    aPowerDer=polyder(aPower);
    %Find the zero of the derivative.
    sRoots=real(roots(aPowerDer));
    sRoots=sRoots(sRoots>=0&sRoots<=1);
    if(isempty(sRoots))
        %Failure.
        xEst=[];
        exitCode=4;
        return;
    end

    if(size(sRoots,1)>1)
        z=polyval(aPower,sRoots);
        [~,idx]=min(abs(1-z));
        sRoots=sRoots(idx);
    end
    x=polyValNewton(sRoots,a(2:end,:),c(2:end,:));
    %Now, we assume that the algorithm is sufficiently close to the
    %solution such the Newton iterations directly on f will converge.
    NewtFunc=@(x)deal(0,f(x));
    %No delta-based testing, no line search, no Cholesky
    %decomposition.
    [xEst,~,exitCode]=NewtonsMethod(NewtFunc,df,x,[AbsTolAns,RelTolAns],0,0,-1,-1,maxIterNewton);
    if(exitCode>=0)
        exitCode=0;
        return;
    else
        xEst=[];%exitCode was set by NewtonsMethod.
        return
    end
end

function [xEst,exitCode]=bracketedEndgame(y,yp,yPrev,ypPrev,H,DH,f,df,endTolF,maxRefineIter,maxIterNewton,AbsTolAns,RelTolAns)
%%BRACKETEDENDGAME After having passed the lambda=1 threshold, the estimate
%           needs to be refined. This implements the engame algorithm
%           described in [1] for the dense Jacobian matrix normal flow
%           algorithm with the modifications described in [2].
%           Additionally, as it was observed that with some very poorly
%           conditioned problems, the final xEst can be sufficiently bad
%           that f(xEst) is not all that close to 0, after the algorithm
%           finishes, if f(xEst) is not very close to 0, then a number of
%           iterations of Newtons algorithm are performed using f and df to
%           refine the estimate. The assumption is that the result from the
%           homotopy algorithm will be close enough for Newton's algorithm
%           ( without step size adjustment) to converge.
%
%REFERENCES:
%[1] L. T. Watson, S. C. Billups, and A. P. Morgan, "ALGORITHM 652 HOMPACK:
%    A suite of codes for globally convergent homotopy algorithms," ACM
%    Transactions on Mathematical Software, vol. 13, no. 3, Sep. 1987, pp.
%    281-310.
%[2] L. T. Watson,  M. Sosonkina, R. C. Melville, A. P. Morgan, and H. F.
%    Walker, "Algorithm 777: HOMPACK90: A Suite of Fortran 90 Codes for
%    Globally Convergent Homotopy Algorithms," ACM transactions on
%    Mathematical Software, vol. 23, no. 4, Dec. 1997, pp. 514-449.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

yDim=size(y,1);
for refineIter=1:maxRefineIter
    if(refineIter==1)
        %Get the interpolation polynomial.
        x=[yPrev,y;
           ypPrev,yp];
        [a,c]=HermiteInterpMultiPoly([0,1],x,yDim);
        aPower=convertPolynomialForm('Newton','PowerSeries',a(1,:),c(1,:));
        %Now, find the value of s such that lambda=1.
        aPower(end)=aPower(end)-1;
        sRoots=real(roots(aPower));
        %Take the first root that is in the range of sPrev and s.
        idx=find(sRoots>=0&sRoots<=1,1);
        if(isempty(idx))
            %Failure.
            xEst=[];
            exitCode=5;
            return;
        end
        sEst=sRoots(idx);
        z=polyValNewton(sEst,a,c);

        p1=y;
        p2=yPrev;
        pOpp=p2;
    else
        %Secant method.
        z=p1+(p2-p1)*(1-p1(1))/(p2(1)-p1(1));

        if(norm(z-p1)>norm(p1-pOpp))
            %Chord method.
            z=p1+(pOpp-p1)*(1-p1(1))/(pOpp(1)-p1(1));
        end
    end

    x=z(2:end);
    fz=f(x);%Function value
    dfz=df(x);%Function Jacobian
    
    Hz=H(z,fz);
    DHz=DH(z,fz,dfz);

    if(any(~isfinite(DHz(:)))||any(~isfinite(Hz)))
        %A non-finite value was encountered.
        xEst=[];
        exitCode=3;
        return
    end
    
    %One iteration of Newton's method.
    DeltaZ=-pinv(DHz)*Hz;
    z=z+DeltaZ;

    if(abs(z(1)-1)<=max(AbsTolAns,RelTolAns))
        scale=max(abs(z));
        diff=max(abs(DeltaZ));
        if(diff<AbsTolAns||diff<scale*RelTolAns)
            xEst=z(2:end);
            %If f has not been sufficienlty zeroed.
            if(max(abs(f(xEst)))>endTolF)
                NewtFunc=@(x)deal(0,f(x));
                %No delta-based testing, no line search, no Cholesky
                %decomposition.
                [xEst,~,exitCode]=NewtonsMethod(NewtFunc,df,xEst,[AbsTolAns,RelTolAns],0,0,-1,-1,maxIterNewton);
            else
                exitCode=0;
            end
            return
        end
    end

    if(z(1)>1)
        if(p1(1)<1)
           pOpp=p1;
        end
    else%z(1)<1
        if(p1(1)>1)
            pOpp=p1;
        end
    end

    p2=p1;
    p1=z;
end

%It did not converge.
xEst=[];
exitCode=6;

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
