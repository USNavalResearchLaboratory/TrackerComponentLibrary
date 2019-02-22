function [xVals,exitCodes]=homotopyPathFollowODE(x0,f,FJacob,c,systemType,startForward,maxSol,initStepSize,order,solutionChoice,RelTol,AbsTol,maxSteps,lambdaTol,numStepsBeforeRestart)
%%HOMOTOPYPATHFOLLOWODE Use a probability-1 homotopy method to either find
%                a fixed point (f(x)=x)of a real, scalar or vector function
%                or to find when the function equals a certain value. Such
%                homotopy methods are useful for finding initial estimates
%                given a minimal number of measurements and essentially no
%                useful initial estimate as the convergence region is
%                significantly larger than Newton's method. Convergence
%                conditions are given in [9]. However, this algorithm is
%                much less reliable than homotopyPathFollow.
%
%INPUTS: x0 An xDimX1 initial estimate to the homotopy problem f(x)=x or
%           f(x)=c. The initial estimate can often be  fairely bad.
%         f A function handle that takes x as a parameter and returns a
%           value having the same dimensionality as x.
%    FJacob A function handle that provides the Jacobian matrix of F. This
%           is the matrix of partial derivatives such that FJacob(i,j) hold
%           the derivative of the ith component of f with respect to the
%           jth component of x. If omitted, the default value using the
%           function.
%         c The value c to use when solving the problem f(x)=c. If one
%           wishes to solve the problem f(x)=x, then this parameter is
%           ignored. If omitted, a default value of zero is used.
% systemType Specified the homotopy that is being solved. Possible values
%           are:
%           0) (The default if omitted). The problem being solved is
%              f(x)=x. The corresponding homotopy is:
%              0=lambda*(x-f(x))+(1-lambda)*(x-a)
%           1) The problem being solved is f(x)=c. The corresponding
%              homotopy is:
%              0=lambda*(c-f(x))+(1-lambda)*(a-f(x))
%              which is equivalent to
%              0=a+lambda*(c-a)-f(x)
%           2) The problem being solved is f(x)=c. The corresponding
%              homotopy is:
%              0=lambda*(c-f(x))+(1-lambda)*(x-a)
%              This homotopy generally seems to not work as well as
%              homotopy 1.
%           3) The problem being solved is f(x)=x. The corresponding
%              homotopy is:
%              0=lambda*(x-f(x))+(1-lambda)*(a-f(x))
%              which is equivalent to 
%              0=a+lambda*(x-a)-f(x)
%              This homotopy seems to not work as well as homotopy 0.
%           4) The problem being solved is f(x)=0. The corresponding
%              homotopy is
%              0=lambda*f(x)+(1-lambda)*((x-a)+(f(x)-f(a)))
%              which is equivalent to
%              0=lambda*((f(a)+a)-x)+x+f(x)-(f(a)+a)
%              or if the f(a)-a terms are lumped together, is equivalent to
%              0=lambda*(a-x)+x+f(x)-a
% startForward An optional parameter indicating whether the homotopy
%           solver should start going forward (in the direction of
%           increasing lambda) or backward (in the direction of decreasing
%           lambda). The default if omitted or an empty matrix is passed is
%           true. One should generally only consider starting backwards if
%           multiple solutions might exist.
%    maxSol Specify the maximum number of solutions that will be attempted
%           to be found. After one solution is found, the function tries to
%           continue on the path in the same direction to find another time
%           when the homotopy parameter lambda=1. Note that finite
%           precision errors can cause it to visit multiple solutions that
%           are essentially equal, or the homotopy path can loop around so
%           that it does revisit past solutions. Also, when it interpolates
%           between steps to find one solution, it does not check whether
%           the interpolating polynomial crosses lambda=1 multiple times,
%           so closely-spaced solutions can also be missed. The default if
%           this parameter is omitted or an empty matrix is passed is 1.
% initStepSize An optional value specifying the initial step size for the
%           adaptive Runge-Kutta method. If omitted or an empty matrix is
%           passed, an initial step size of 0.1 is used. Depending on the
%           scale of the problem and the inaccuracy of the initial
%           estimate, that might be too small. However, in general, it
%           usually works fairely well.
% order,solutionChoice  A pair of optional parameters that specify the
%           highest order of the embedded Runge-Kutta pair or Rosenbrock
%           method to use as well as the specific algorithm to use. Details
%           are given in the comments to the RungeKStep function for Runge
%           Kutta methods. To use a Rosenbrock method, set
%           solutionChoice=-1 and set the appropriate order for the method.
%           If omitted or empty matrices are passed, the default order of 5
%           is used and the default solutionChoice of 0 is used (A Runge-
%           Kutta method).
%    RelTol The maximum relative error tolerance allowed, a positive scalar
%           (its use is explained in more detail below). If omitted or an
%           empty matrix is passed, the default value of 1e-8 is used.
%    AbsTol The absolute error tolerance allowed, a positive scalar, of the
%           same for all components of x, or a positive NX1 vector. If
%           omitted or an empty matrix is passed, the default value of
%           1e-10 is used.
%  maxSteps The maximum allowable number of steps to perform the
%           integration. If omitted or an empty matrix is passed, the
%           default of 512 is used.
% lambdaTol The tolerance for declaring the homotopy parameter 1. If
%           omitted, the default value of 1e-8 is used.
%numStepsBeforeRestart As the function goes along the path, finite
%           precision errors will case the homotopy condition to be less
%           and less satisfied. THis parameter specified how many
%           iterations should go by before the homotopy is reset. The
%           default if omitted is 25. To reset the homotopy every
%           iteration, set it to 0.
%
%OUTPUTS: xVals The xDimXnumSol solutions to the chosen homotopy problem,
%               or an empty matrix if no solution can be found. Note that
%               this function only provides at most maxSol solution.
%               Depending on f(x), there might be multiple solutions or
%               zero solutions.
%     exitCodes A value indicating the exit status of the solver. When
%               asking for multiple solutions, this can be one longer than
%               xVal to indicate why the final solution method failed.
%               Possible values are
%               0: No error.
%               1: The step size decreased below the minimum allowed level
%                  while trying to reach the desired tolerances.
%               2: Maximum number of steps reached without completion.
%               3: Non-finite number encountered.
%               4: An error occurred interpolating the final solution.
%               5: A dependent function threw an exception. This generally
%                  means that a NaN was encountered in the null function
%                  when trying to solve for a descent direction.
%
%The algorithm is implemented as discussed in [1]. The algorithm is based
%on the probability 1 homotopy theory, which transforms the problem of
%solving difficult nonlinear equations into a problem of following a curve
%until a particular parameter becomes 1. Such techniques reach the desired
%solution with probability 1 given a completely random initialization and
%infinite precision. However, such proofs on probability-1 convergence are
%often done in the complex domain and this function operates in the real
%domain. Nonetheless, even though in reality the algorithm can fail with
%certain initializations and due to the use of finite precision, it is
%still robust with regard to a number of problems.
%
%The theory of probability 1 homotopy algorithms is discussed in [2] and
%[3]
%
%The algorithm implemented here is based somewhat on the description of the
%dense differential-equation based routine of [4], which is updated in [5].
%An understanding of some of the parts can also be obtained from [6].
%The dense differential equation based method is chosen despite being
%slower than alternatives using Newton's method, because the paper for
%algorithm 652 describes it as being less likely to fail.
%
%Basically, the homotopy chosen is differentiated with respect to a
%pathlength term s of which the state of interest x and the homotopy
%parameter lambda depend. Lambda starts at zero, and the path is followed
%(the differential equation is integrated) using an adaptive Runge-Kutta
%method. When value of lambda reaches 1, the problem is solved.
%
%However, as described in the above papers, finite precision errors can
%make the integrated path diverge from the homotopy condition that the
%chosen homotopy equation be zero. In such an instance, the algorithmic
%constant is reset periodically and near the end to force the homotopy
%conditions to remain satisfied. Also, at the end, as it is assumed that
%one will overstep the point of lambda=1, interpolation is performed to
%find the actual point where lambda=1. The interpolation is performed
%using Newton's method as described in Chapter 6.6 of [7].
%
%Homotopy methods are useful for many problem related to target tracking.
%For example, the problem of range-only track initiation is considered in 
%[8].
%
%REFERENCES:
%[1] D. F. Crouse, "Target Track Initiation in Difficult Scenarios Using
%    Probability-1 Homotopy Methods and Cubature Integration," in
%    Proceedings of the IEEE Aerospace Conference, Big Sky, MT, Mar. 2016.
%[2] L. T. Watson, "Globally convergent homotopy methods: A tutorial,"
%    Applied Mathematics and Computation, vol. 31, pp. 369-396, May 1989.
%[3] S. L. Richter and R. A. DeCarlo, "Continuation methods: Theory and
%    applications," IEEE Transactions on Systems, Man, and Cybernetics,
%    vol. SMC-13, no. 4, pp. 459-464, Jul./Aug. 1983.
%[4] L. T. Watson, S. C. Billups, and A. P. Morgan, "Algorithm 652:
%    HOMPACK: A suit of codes for globally convergent homotopy algorithms,"
%    ACM Transactions on Mathematical Software, vol. 13, no. 3, pp. 281-
%    310, Sep. 1987.
%[5] L. T. Watson, M. Sosonkina, R. C. Melville, A. P. Morgan, and H. F.
%    Walker, "Algorithm 777: HOMPACK: A suit of Fortran 90 codes for
%    globally convergent homotopy algorithms," ACM Transactions on
%    Mathematical Software, vol. 23, no. 4, pp. 514-549, Dec. 1997.
%[6] L. T. Watson, "A globally convergent algorithm for computing fixed
%    points of C2 maps," Applied Mathematics and Computation, vol. 5, no.
%    4, pp. 297-311, Oct. 1979.
%[7] J. R. Dormand, Numerical Methods for Differential Equations. Boca
%    Raton: CRC PRess, 1996.
%[8] R. L. Smith and C. Huang, "Study of a homotopy continuation method for
%    early orbit determination with the tracking and data relay satellite
%    system (TDRSS)," National Aeronautics and Space Administration/
%    Computer Sciences Corporation, Beltsville, MD, Tech. Rep. 86230,
%    Mar. 1986.
%[9] L. T. Watson, "Probability-One Homotopies in Computational Science,"
%    Journal of Computational and Applied Mathematics, vol. 140, no. 1-2.
%    Mar. 2002, pp. 785-807.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=length(x0);
    
    if(nargin<15)
        numStepsBeforeRestart=25;
    end
    
    if(nargin<14)
       lambdaTol=1e-9;
    end
    
    if(nargin<13||isempty(maxSteps))
        maxSteps=512;
    end
    
    if(nargin<12||isempty(AbsTol))
        AbsTol=1e-10;
    end
    
    if(nargin<11||isempty(RelTol))
        RelTol=1e-9;
    end
    
    if(nargin<=10||isempty(solutionChoice))
        solutionChoice=0;
    end
    
    if(nargin<9||isempty(order))
        order=5;
    end
    
    if(nargin<8||isempty(initStepSize))
        %The maximum distance that can be traveled along the curve in one
        %step.
        initStepSize=0.1;
    end
    
    if(nargin<7||isempty(maxSol))
        maxSol=1;
    end
    
    if(nargin<6||isempty(startForward))
       startForward=true; 
    end
    
    if(nargin<5||isempty(systemType))
        systemType=0;
    end
    
    if(nargin<4||isempty(c))
        c=zeros(xDim,1);
    end
    
    if(nargin<3||isempty(FJacob))
       FJacob=@(x)numDiff(x,f,xDim); 
    end
    
    if(systemType~=4)
        %The initial value of a.
        a=f(x0);
        %Change x0 as needed to make the initial homotopy equation true.
        if(systemType~=1&&systemType~=3)
            x0=a;
        end
    else
        a=f(x0)+x0;
    end
    
    %The initial step size.
    deltaS=initStepSize;
    
    %We have not yet travered any distance along the path.
    s=0;
    %The initial augmented state. We must integrate until the homotopy
    %parameter (the first component) becomes 1.
    u=[0;x0];
    
    %The point of setting dudsPrev to something even though there is no
    %first step, is because the derivative function needs the previous
    %derivative to ensure continuity.
    dudsPrev=zeros(xDim+1,1);
    if(startForward==false)
        dudsPrev(1)=-1;
    else
        dudsPrev(1)=1;
    end
    
    %Allocate the maximum amount of space for the return values.
    xVals=zeros(xDim,maxSol);
    exitCodes=zeros(maxSol,1);
    for curSol=1:maxSol
        [xVal,exitCode,s,deltaS,u,dudsPrev]=takeOneStep(a,s,deltaS,u,dudsPrev,f,FJacob,c,systemType,order,solutionChoice,RelTol,AbsTol,maxSteps,lambdaTol,numStepsBeforeRestart);
        
        exitCodes(curSol)=exitCode;
        %If the algorithm did not fail.
        if(exitCode==0)
            xVals(:,curSol)=xVal;
            a=getAForHomotopy(f,c,u,systemType);
        else
            %Shrink to fit the solutions.
            xVals=xVals(:,1:(curSol-1));
            exitCodes=exitCodes(1:curSol);
            return;
        end
    end
end


function [xVal,exitCode,s,deltaS,u,dudsPrev]=takeOneStep(a,s,deltaS,u,dudsPrev,f,FJacob,c,systemType,order,solutionChoice,RelTol,AbsTol,maxSteps,lambdaTol,numStepsBeforeRestart)
    xDim=length(u)-1;
    
    %The maximum number of iterations of Newton's method to perform when
    %trying to interpolate to the lambda=1 between two points stradding the
    %correct lambda.
    numNewtonIter=25;

    %When the algorithm resets, that means that it computes a new value for
    %a to force the chosen homotopy to be true.
    
    %Counts the number of times the algorithm has restarted with a new x0
    %value. At least one reset should occur when lambda is close to 1 to
    %ensure all constraints are satisfied. Also, if it makes a lot of steps
    %without making any progress, then it will also reset.
    nearEndResetOccurred=false;
    
    numStepsSinceReset=0;
    
    %The algorithm is put in a try-catch statement so that if it fails, the
    %failure can be easily reported. Failures that raise exceptions that
    %would need to be caught are typically going to come from the null
    %function encountering a NaN. NOTE that when debugging this, the
    %try-catch framework will catch all actual errors in the code that
    %would normally make Matlab stop. Thus, it should be commented out if
    %the code below needs to be modified/ debugged.
    try
        %Get the initial value of the derivative.
        derivVal=diffEq(u,a,dudsPrev,f,FJacob,c,systemType);
        for curStep=1:maxSteps        
            derivFunc=@(u,s)diffEq(u,a,dudsPrev,f,FJacob,c,systemType);

            %The previous values are saved in case they new to be set if a
            %reset occurs.
            uPrev=u;
            sPrev=s;
            derivValPrev=derivVal;
            dudsPrevPrev=dudsPrev;
            [deltaS,u,s,k,derivVal,exitCode]=performOneStep(u,s,derivFunc,deltaS,derivVal,order,solutionChoice,RelTol,AbsTol);

            %Check for an error taking the adaptive Runge-Kutta step.
            if(exitCode~=0)
                xVal=[];
                return;
            end
            %Extract the value of duds for uPrev.
            dudsPrev=k(:,1);

            numStepsSinceReset=numStepsSinceReset+1;

            %Make it reset when near completion so that the solution is
            %(hopefully) more accurate.
            if(u(1)>=0.99&&nearEndResetOccurred==false)
                shouldResetForEnd=true;
            else
                shouldResetForEnd=false;
            end

            %Determine whether the algorithm can terminate. This looks to
            %determine whether the sign changed. That is, whether it has
            %stepped over lambda=1.
            signChanged=(u(1)>=1&&uPrev(1)<1)||(u(1)<=1&&uPrev(1)>1);
            
            %Alternatively, the algorithm can terminate if the last step
            %was not within the tolerance bounds, and the current step is.
            %the condition that the last step not be within the bounds lets
            %the algorithm be restarted and continue along the path after
            %it finds a solution.
            enteredTolBound=(abs(u(1)-1)<=lambdaTol)&&(abs(uPrev(1)-1)>lambdaTol);
            
            if(nearEndResetOccurred==true&&(signChanged||enteredTolBound))
            %If here, we are at or have overstepped the end. We will
            %interpolate to the point where lambda=u(1)=1 if the value of
            %lambda is not within the tolerance set for lambda.
                if(enteredTolBound==false)
                    %Get the interpolating polynomials
                    fInterp=@(u,s)derivFunc(u,s);
                    [interpPolyA,interpPolyC]=RKInterpPolys(uPrev,sPrev,u,s,fInterp,order,solutionChoice,k);

                    %Now, use Newton's method to find where lambda=u(1)=1.
                    %Note that the interpolation polynomial takes an
                    %argument from 0->1 indicating the fractional distance
                    %along the path from sPrev to s.
                    sigmaEst=(1-uPrev(1))/(u(1)-uPrev(1));%Initial guess.

                    lambdaCur=Inf;
                    curIter=1;
                    while(abs(lambdaCur-1)>lambdaTol)
                        %The Newton iteration should converge quickly. If
                        %it does not, then exit with an error.
                        if(curIter>numNewtonIter)
                            xVal=[];
                            exitCode=4;
                            return;
                        end

                        [lambdaDerCur,lambdaCur]=polyDerValNewton(sigmaEst,interpPolyA(1,:),interpPolyC(1,:),1);
                        sigmaEst=sigmaEst-(lambdaCur-1)/lambdaDerCur;

                        curIter=curIter+1;
                    end

                    %If a good interpolation value was found, then return
                    %it.
                    xVal=zeros(xDim,1);
                    for curDim=1:xDim
                        xVal(curDim)=polyValNewton(sigmaEst,interpPolyA(curDim+1,:),interpPolyC(curDim+1,:));
                    end
                else
                %The end value is already close enough, so no interpolation
                %is necessary.
                    xVal=u(2:end);
                end

                %Successful exit.
                exitCode=0;
                return;
            elseif(shouldResetForEnd||numStepsSinceReset>numStepsBeforeRestart)
            %The algorithm should be restarted. If it is the first restart, it
            %might have overshot the end. We do not want to approach lambda=1
            %by decreasing lambda, so in that specific case, we use uPrev and
            %sPrev instead of u and s. In any reset, the value dudsPrev is not
            %changed, in the hope that the sense of the curve is not flipped
            %after the reset.
                if(shouldResetForEnd)
                    nearEndResetOccurred=true;
                end

                if(u(1)>=1)
                    deltaS=s-sPrev;
                    s=sPrev;
                    u=uPrev;
                    derivVal=derivValPrev;
                    dudsPrev=dudsPrevPrev;
                end
                
                a=getAForHomotopy(f,c,u,systemType);
            
                numStepsSinceReset=0;
                continue;
            end
        end

        %If we get here, then the maximum number of iterations was reached 
        %without reaching the end.
        xVal=[];
        exitCode=2;
        return;
    catch
        %If some function raised an exception, catch it and record the
        %error.
        xVal=[];
        exitCode=5;
        return;
    end
end

function a=getAForHomotopy(f,c,u,systemType)
%%GETAFORHOMOTOPY Get the constant value a that forces the homotopy
%                 specified by systemType to be satisfied for a given u
%                 vector.

    lambdaCur=u(1);
    xCur=u(2:end);
    if(systemType==0)%Make 0=lambda*(x-f(x))+(1-lambda)*(x-a) hold.
        a=(lambdaCur/(1-lambdaCur))*(xCur-f(xCur))+xCur;
    elseif(systemType==1)
        %Make 0=a+lambda*(c-a)-f(x) hold.
        a=(f(xCur)-c*lambdaCur)/(1-lambdaCur);
    elseif(systemType==2)
        %Make 0=lambda*(c-f(x))+(1-lambda)*(x-a) hold.
        a=lambdaCur/(1-lambdaCur)*(c-f(xCur))+xCur;
    elseif(systemType==3)
        %Make 0=a+lambda*(x-a)-f(x) hold.
        a=(f(xCur)-xCur*lambdaCur)/(1-lambdaCur);
    else
        %Make 0=lambda*(a-x)+x+f(x)-a hold.
        a=(f(xCur)+(1-lambdaCur)*xCur)/(1-lambdaCur);
    end
end


function duds=diffEq(u,a,dudsPrev,f,FJacob,c,systemType)
%DIFFEQ The differential equation corresponding to the selected homotopy.
    lambda=u(1);
    x=u(2:end);
    xDim=length(x);

    %Using the nullspace method rather than the explicit solution for
    %the derivative.
    if(systemType==0)
        %For 0=lambda*(x-f(x))+(1-lambda)*(x-a)
        D=[a-f(x),eye(xDim,xDim)-lambda*FJacob(x)];
    elseif(systemType==1)
        %For 0=a+lambda*(c-a)-f(x)
        D=[c-a,-FJacob(x)];
    elseif(systemType==2)
        %For 0=lambda*(c-f(x))+(1-lambda)*(x-a)
        D=[c-f(x)-(x-a),(1-lambda)*eye(xDim,xDim)-lambda*FJacob(x)];
    elseif(systemType==3)
        %For 0=a+lambda*(x-a)-f(x)
        D=[x-a,lambda*eye(xDim,xDim)-FJacob(x)];
    else
        %For 0=lambda*(a-x)+x+f(x)-a
        D=[a-x,(1-lambda)*eye(xDim,xDim)+FJacob(x)];
    end
    
    %If NaN values are present, this might throw an exception.
    temp=null(D);
    
    %In the event of a bifurcation, just take the first one.
    v=temp(:,1);
    duds=v/norm(v);
    %Ensure continuity along the curve. --i.e. get the sign correct so
    %that it does not backtrack on itself.
    if(duds'*dudsPrev<0)
        duds=-duds;
    end
end

function [deltaS,u,s,k,derivVal,exitCode]=performOneStep(u,s,derivFunc,deltaS,derivVal,order,solutionChoice,RelTol,AbsTol)
%%PERFORMONESTEP Take one step in a Runge-Kutta method with an adaptive
%                  stepsize. Much of this function is similar to one step
%                  in the function RKAdaptiveOverRange.

    uDim=size(u,1);
    %Determine the orders of the main and subsidiary embedded formulae that
    %were chosen. The order of convergence for the step size testing is the
    %smallest order taken.
    if(solutionChoice==-1)
        [orders,isFSAL]=RosenbrockStep(order);
        %Numerically find the Jacobian. The numjac function is used.
        fRev=@(t,x)derivFunc(x,t);
        dxdtCur=derivFunc(u,s);
        [J] = numjac(fRev,s,u,dxdtCur,AbsTol,[],false);
    else
        [orders,isFSAL]=RungeKStep(order,solutionChoice);
    end
    RKOrder=min(orders);
    
    %The minimum acceptable stepsize
    deltaSMin=2^4*eps(s);

    %Since th minimum step size changes every step, this makes sure 
    %that deltaS does not go beneath it just because the step changed.
    deltaS=max(deltaS,deltaSMin);

    %The first time the choice in step size fails, it is adjusted the
    %"optimal" way in the Runge-Kutta-Fehlberg method. Additional times,
    %the step size is just halved in the hope that it will reach an
    %accepted value more quickly.
    failedReducingStepSize=false;
    moveOnToNextStep=false;
    while(moveOnToNextStep==false)
        %The returned value dxValdt is for xCur and tCur, which is
        %curStep-1.
        if(solutionChoice==-1)
            [uPredMain,uPredSubsid,k]=RosenbrockStep(u,s,derivFunc,deltaS,derivVal,J,zeros(uDim,1),order);
        else
            [uPredMain,uPredSubsid,k]=RungeKStep(u,s,derivFunc,deltaS,derivVal,order,solutionChoice);
        end
        %Integration can only be over finite functions.
        if(any(~isfinite(uPredMain))||any(~isfinite(uPredSubsid)))
            exitCode=3;
            return;
        end

        %The local error estimate. This must be transformed into a 
        %combination relative/ absolute error term to determine whether
        %the step should be rejected.
        normFactor=max(max(abs(uPredMain),abs(u)),AbsTol/RelTol);
        theError=max(abs((uPredMain-uPredSubsid)./normFactor));

        if(theError>RelTol)
            %The step should be rejected.
            if(deltaS<deltaSMin)
                %If the step size got too small, then return.
                exitCode=1;
                return;
            end

            if(failedReducingStepSize==false)
                failedReducingStepSize=true;

                %The Fehlberg step reduction (using the relative error).
                deltaS=max(deltaSMin, deltaS * max(0.1, 0.8*(RelTol/theError)^(1/RKOrder)));
            else
                %Just halve the step size.
                deltaS=deltaS/2;
            end
        else
            %Save the results from the step.
            s=s+deltaS;

            %If a step is successful, then increase the step size for the 
            %next step in the standard manner used with Runge-Kutta-
            %Fehlberg methods, but limit the maximum size of the increase
            %to a scale factor of 10. This avoid huge step sizes when the
            %predicted error is very small. On the other hand, if the scale
            %of the problem is very huge, it might take a lot of steps to
            %get the step size up to a large enough number.
            deltaS=deltaS*min(10,0.8*(RelTol/theError)^(1/RKOrder));
            u=uPredMain;
            
            moveOnToNextStep=true;
        end
    end
    
    %Get the derivative at the current time to return.
    if(isFSAL)
        derivVal=k(:,end);
    else
        derivVal=derivFunc(u,s);
    end
    
    %No error on exit.
    exitCode=0;
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
