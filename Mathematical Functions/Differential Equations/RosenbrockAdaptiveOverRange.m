function [xVals,tVals,dxdtVals,exitCode,numRejections]=RosenbrockAdaptiveOverRange(xStart,tSpan,f,order,JacobianFun,dfdtFun,initStepSize,RelTol,AbsTol,maxSteps)
%%ROSENBROCKADAPTIVEOVERRANGE Integrate an ordinary differential equation
%                   using a modified Rosenbrock method with an adaptive
%                   step size. That is, integrating dx/dt=f(x,t) given
%                   initial conditions (xStart,tSpan(1)). Rosenbrock
%                   methods are better then Runge-Kutta methods when the
%                   differential equations are stiff. However, Rosenbrock
%                   methods require a derivative. More information
%                   on the algorithms used with different orders is given
%                   in the comments to the function RosenbrockStep.
%
%INPUTS: xStart The NX1 state vector at time tSpan(1).
%         tSpan  A 2X1 or 1X2 vector where tSpan(1) is the starting time
%                and tSpan2 is the desired stopping time for the
%                integration.
%             f  The function handle for f(x,t)=dxdt over which integration
%                is to be performed.
%          order The order of the main Rosenbrock routine to use. See the
%                function RosenbrockStep for possible values and algorithms
%                used as this function calls RosenbrockStep. If omitted or
%                an empty matrix is passed, the default value of 2 is used.
%   initStepSize An optional initial step size (in t) to use for the
%                integration. If omitted or an empty matrix is passed, an
%                ad-hoc method described below is used to find an initial
%                step size.
%    JacobianFun This is either a function handle having the form
%                JacobianFun(x,t) that provides the Jacobian of the
%                function f --the derivatives with respect to the parameter
%                x-- or this is a matrix that is the constant Jacobian
%                matrix. If this parameter is omitted or an empty matrix is
%                passed, then the NXN Jacobian matrix is computed using the
%                numjac function with the tolerance set to AbsTol.
%        dfdtFun This is either a function handle having the form
%                dfdtFun(x,t) that provides the derivative of f with
%                respect to the scalar parameter t, or this is a constant
%                NX1 vector for the derivative of f with respect to t. In
%                autonomous problems, a zero matrix should be passed. If
%                omitted, a central difference formula is used to
%                numerically approximate the derivative vector.
%   initStepSize An optional initial step size (in t) to use for the
%                integration. If omitted or an empty matrix is passed, an
%                ad-hoc method described below is used to find an initial
%                step size.
%         RelTol The maximum relative error tolerance allowed, a
%                positive scalar (its use is explained in more detail
%                below). If omitted or an empty matrix is passed, the
%                default value of 1e-3 is used.
%         AbsTol The absolute error tolerance allowed, a positive scalar,
%                of the same for all components of x, or a positive NX1
%                vector. If omitted or an empty matrix is passed, the
%                default value of 1e-6 is used.
%       maxSteps The maximum allowable number of steps to perform the
%                integration. If omitted, the default of 4096 is used.
%
%OUTPUTS:xVals The NXnumSteps set of values of x along the path.
%                  xVals(:,1) is xStart and xVals(:,end) is the value at
%                  final time tSpan(2). If the integration failed, e.g. due
%                  to encountering a NaN or being unable to get a sufficient
%                  step size, an empty matrix is returned.
%            tVals A numStepsX1 vector of values of t corresponding to
%                  the values of x in xVals. tVals(1) is equal to tSpan(1)
%                  and tVals(end) is equal to tSpan(2). If integration
%                  fails, an empty matrix is returned.
%         dxdtVals A numStepsX1 array containing deivatives of x evaluated
%                  at the times in tVals. The derivatives are just the
%                  result of evaluating f(x,t) at each point in
%                  (xVals,tVals). This combined with xVals and tVals could
%                  be used in functions to perform Hermite interpolation.
%                  If integration fails, an empty matrix is returned.
%         exitCode A code indicating how the algorithm terminated. Possible
%                  values are
%                  0: Integration was successful.
%                  1: Unable to get a small enough step size.
%                  2: Maximum number of steps reaced without completion.
%                  3: Non-finite number encountered.
%    numRejections The number of times a stepsize hypothesis was rejected.
%
%The same type of adaptive stepsize control is used with the Rosenbrock
%algorithm (a diagonal implicit Runge-Kutta routine) as with standard
%Runge-Kutta methods. Thus, the comments to the stepsize adaptation routine
%in RKAdaptiveOverRange explain how the stepsize is adjusted here. The onyl
%notable change is that a certain matrix has to be inverted for each of the
%Rosenbrock methods. The matrix depends on the stepsize. Thus, if the
%matrix is singular, this function has an additional rejection of the step
%and the stepsize is halved.
%
%The central difference formula used for the numeric derivative of f with
%respect to t if dfdtFun is not provided is taken from Chapter 4.1 of [1],
%whereby the three-point midpoint formula is used.
%
%REFERENCES:
%[1] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/ Cole, 2011.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10)
    maxSteps=4096;
end

if(nargin<9||isempty(AbsTol))
    AbsTol=1e-6;
end

if(nargin<8||isempty(RelTol))
   RelTol=1e-3;
end

xDim=size(xStart,1);
tStart=tSpan(1);
tEnd=tSpan(2);
tDiff=diff(tSpan);
tDiffMag=abs(tDiff);
deltaTSign=sign(tDiff);
%The sign is separated out so that one can use the algorithm running
%backwards as well as forwards in time.

%The maximum step size is arbitraily set to 1/3 the total distance.
maxStepSize=tDiffMag/3;

if(nargin<7||isempty(initStepSize))
    %If no initial step size is given, then just use the smallest uniform
    %step size for the given value of maxSteps.
    initStepSize=tDiffMag/maxSteps;
elseif(initStepSize<eps(tStart))
   error('The initial step size provided is too small') 
end

if(nargin<6)
    %Tell it to perform numerical differentiation for dfdt.
    dfdtFun=[];
end

if(nargin<4||isempty(order))
    order=2; 
end

%The lowest-order formula is used for the convergence rate. The isFSAL flag
%indicates whether k(:,end) is the value of f evaluated at the next step.
[orders,isFSAL]=RungeKStep(order);
RKOrder=min(orders);

%Allocate the maximum possible space needed for the return parameters.
xVals=zeros(xDim,maxSteps);
dxdtVals=zeros(xDim,maxSteps);
tVals=zeros(maxSteps,1);

xVals(:,1)=xStart;
tVals(:,1)=tStart;
dxdtVals(:,1)=f(xStart,tStart);

if(nargin<5||isempty(JacobianFun))
    %Indicate that the jacobian should be determined numerically
    JacobianType=0;
    FAC=[];%A pointer to store results for the jacobian function between
           %calls to make it faster.
    
    %The numjac function requires a minimum tolerance that is a vector.
    if(isscalar(AbsTol))
       AbsTol=repmat(AbsTol,[xDim,1]); 
    end
elseif(isa(JacobianFun,'function_handle'))
    %A function handle is provided for the Jacobian.
    JacobianType=1;
else
    JacobianType=2;
    %A constant Jacobian matrix is provided.
	J=JacobianFun;
end

deltaTMag=initStepSize;
numRejections=0;
for curStep=2:maxSteps
    xCur=xVals(:,curStep-1);
    tCur=tVals(curStep-1);
    dxdtCur=dxdtVals(:,curStep-1);
        
    %The absolute value of the minimum allowable step size. It is set so
    %that the step must makes something of a difference compared to the
    %numerical precision.
    deltaTMinMag=2^4*eps(tCur);
    %Since th minimum step size changes every loop, this makes sure that
    %deltaT does not go beneath it just because the loop changed.
    deltaTMag=max(deltaTMag,deltaTMinMag);
    deltaT=deltaTMag*deltaTSign;
    
    %If we would have overstepped the end, change deltaT to be the end.
    %Then, we can terminate at the end.
    if(deltaTSign*(tCur+deltaT)>tEnd*deltaTSign)
        deltaTMag=abs(tEnd-tCur);
        deltaT=deltaTMag*deltaTSign;
        %Allow the last step to be very small.
        deltaTMinMag=min(deltaTMinMag,deltaTMag);
    end
    
    if(JacobianType==0)
        %Numerically find the Jacobian. The numjac function is used.
        fRev=@(t,x)f(x,t);
        [J,FAC] = numjac(fRev,tCur,xCur,dxdtCur,AbsTol,FAC,false);
    elseif(JacobianType==1)
        %User-provided function
        J=JacobianFun(xCur,tCur);
    end

    %The first time the choice in step size fails, it is adjusted the
    %"optimal" way in the Runge-Kutta-Fehlberg method. Additional times,
    %the step size is just halved in the hope that it will reach an
    %accepted value more quickly.
    failedReducingStepSize=false;
    moveOnToNextStep=false;
    while(moveOnToNextStep==false)
        %If an explicit function (or matrix) for the derivative with
        %respect to time is provided.
        if(isempty(dfdtFun))
            %Just use a central finite-difference formula with a set being
            %a quarter of the current step size.
            dfdt=(f(xCur,tCur+deltaT/8)-f(xCur,tCur-deltaT/8))/(2*8);
        elseif(isa(dfdtFun,'function_handle'))
            dfdt=dfdtFun(xCur,tCur);
        else%It is a constant matrix.
            dfdt=dfdtFun;
        end

        [xPredMain,xPredSubsid,k]=RosenbrockStep(xCur,tCur,f,deltaT,dxdtCur,J,dfdt,order);    
        %If the step failed because the W matrix was singular, reduce the
        %step size and the matrix should not stay singular.
        if(isempty(xPredMain))
            numRejections=numRejections+1;
            
            if(deltaTMag<deltaTMinMag)
                %If the step size got too small, then return.
                xVals=[];
                dxdtVals=[];
                tVals=[];
                exitCode=1;
                return;
            end
            
            %Just halve the step size to make W non-singular.
            deltaTMag=deltaTMag/2;
            deltaT=deltaTSign*deltaTMag;
            continue;
        end
        
        %Integration can only be over finite functions.
        if(any(~isfinite(xPredMain))||any(~isfinite(xPredSubsid)))
            xVals=[];
            tVals=[];
            dxdtVals=[];
            exitCode=3;
            return;
        end
        
        %The step-size adaptation is the same as in RKAdaptiveOverRange
        %The local error estimate. This must be transformed into a 
        %combination relative/ absolute error term to determine whether
        %the step should be rejected.
        normFactor=max(max(abs(xPredMain),abs(xCur)),AbsTol/RelTol);
        theError=max(abs((xPredMain-xPredSubsid)./normFactor));
        if(theError>RelTol)
            %The step should be rejected.
            numRejections=numRejections+1;

            if(deltaTMag<deltaTMinMag)
                %If the step size got too small, then return.
                xVals=[];
                dxdtVals=[];
                tVals=[];
                exitCode=1;
                return;
            end

            if(failedReducingStepSize==false)
                failedReducingStepSize=true;

                %The Fehlberg step reduction (using the relative error).
                deltaTMag=max(deltaTMinMag, deltaTMag * max(0.1, 0.8*(RelTol/theError)^(1/RKOrder)));
            else
                %Just halve the step size.
                deltaTMag=deltaTMag/2;
            end
            deltaT=deltaTSign*deltaTMag;
        else
            %If a step is successful, then increase the step size for the 
            %next step in the standard manner used with Runge-Kutta-
            %Fehlberg methods, but limit the maximum size of the increase
            %to a scale factor of 4. This avoid huge step sizes when the
            %predicted error is very small.
            deltaTMag=min(maxStepSize,deltaTMag*min(4,0.8*(RelTol/theError)^(1/RKOrder)));

            %Save the results from the step.
            xVals(:,curStep)=xPredMain;
            tVals(curStep)=tCur+deltaT;
            if(isFSAL)
                dxdtVals(:,curStep)=k(:,end);
            else
                dxdtVals(:,curStep)=f(xVals(:,curStep),tVals(curStep));
            end
            
            moveOnToNextStep=true;
        end
    end

    if(tVals(curStep)==tEnd)
        %If integration completed, reshape the outputs to match the actual
        %number of steps.
        xVals=xVals(:,1:curStep);
        tVals=tVals(1:curStep);
        dxdtVals=dxdtVals(:,1:curStep);
        exitCode=0;
        return
    end
end

%If we get here, then the maximum number of iterations was reached without
%reaching the end.
xVals=[];
dxdtVals=[];
tVals=[];
exitCode=2;
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
