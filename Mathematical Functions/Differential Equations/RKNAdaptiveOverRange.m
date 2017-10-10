function [xVecVals,tVals,d2xdt2Vals,exitCode,numRejections]=RKNAdaptiveOverRange(xVecStart,tSpan,df,probIsGeneral,numRefinedSteps,initStepSize,order,solutionChoice,RelTol,AbsTol,maxSteps)
%%RKNADAPTIVEOVERRANGE Perform explicit general or special Runge-Kutta-
%                      Nyström integration over a given range of values
%                      using an adaptive stepsize. Runge-Kutta Nyström
%                      methods are derivative-free techniques for
%                      integrating second order differential equations. The
%                      general formulation allows the integration of
%                      equations of the form d^2xdt^2=df(xVec,t), where
%                      xVec consists of the stacked values x and dxdt. In
%                      the special formulation d^2xdt^2=df(x,t), where no
%                      first derivative values are passed. More information
%                      on Runge-Kutta-Nyström integration is given in the
%                      comments to the functions RungeKNystroemGStep and
%                      RungeKNystroemSStep.
%
%INPUTS:xVecStart The (2N)X1 state vector (consisting of N components and
%                their N derivatives with respect to t) at time tSpan(1)
%                (xVecStart=[x;dxdt]).
%          tSpan A 2X1 or 1X2 vector where tSpan(1) is the starting time
%                and tSpan2 is the desired stopping time for the
%                integration.
%             df The function df returns the second derivatuve d2xdt2. If 
%                the general problem is being solved, then
%                probIsGeneral=true and df is called as df(xVec,t) where
%                the full vector of N values and their N first derivatives
%                is passed. If probIsGeneral=false, then df is called as
%                df(x,t), where x is just the values without their
%                derivatives.
%  probIsGeneral Indicates whether the problem being solved is the general
%                problem, where df depends on the frist derivatives, of the
%                special problem, where df does not depend on the first
%                derivatives. The default if omitted or am empty matrix is
%                passed is true.
%numRefinedSteps An optional parameter specifying the number of steps to
%                interpolate between each normal step. With high-order
%                formulas, the steps produced might not be very nice to
%                use for plotting, so setting this to a value above 0 adds
%                more intermediate values (using the polynomial returned by
%                the RKNInterpPolys function) so that plotted data can look
%                nicer. If omitted or an empty matrix is passed, the
%                default value of 0 (no extra/ interpolated steps) is used.
%                Interpolation is done at an order that ranges between the
%                main and subsidiary orders of the Runge-Kutta-Nyström
%                method.
%   initStepSize An optional initial step size (in t) to use for the
%                integration. If omitted or an empty matrix is passed, an
%                ad-hoc method described below is used to find an initial
%                step size.
% order,solutionChoice  A pair of optional parameters that specify the
%                highest order of the embedded Runge-Kutta pair to use as
%                well as the specific algorithm to use. Details are given
%                in the comments to the RungeKNystroemGStep function, for
%                the general problem, and in the RungeKNystroemSStep
%                function for the special problem. If omitted or empty
%                matrices are passed, the default order of 5 is used
%                and the default solutionChoice of 0 is used.
%         RelTol The maximum relative error tolerance allowed, a
%                positive scalar (its use is explained in more detail
%                below). If omitted or an empty matrix is passed, the
%                default value of 1e-3 is used.
%         AbsTol The absolute error tolerance allowed, a positive scalar,
%                of the same for all components of x, or a positive NX1
%                vector. If omitted or an empty matrix is passed, the
%                default value of 1e-6 is used.
%       maxSteps The maximum allowable number of steps to perform the
%                integration. If omitted, the default of 1024 is used.
%
%OUTPUTS: xVecVals The NXnumSteps set of values of x along the path.
%                  xVecVals(:,1) is xVecStart and xVecVals(:,end) is the
%                  value at the final time tSpan(2). If the Runge-Kutta-
%                  Nyström integration fails, e.g. due to encountering a
%                  NaN or being unable to get a sufficient step size, an
%                  empty matrix is returned.
%            tVals A numStepsX1 vector of values of t corresponding to
%                  the values of x in xVecVals. tVals(1) is equal to tSpan(1)
%                  and tVals(end) is equal to tSpan(2). If integration
%                  fails, an empty matrix is returned.
%       d2xdt2Vals An NXnumSteps array containing the second derivatives of
%                  x evaluated at the times in tVals. The derivatives are
%                  just the result of evaluating df(x,t) at various points.
%                  This combined with xVecVals and tVals could be used in
%                  functions to perform Hermite interpolation. If
%                  integration fails, an empty matrix is returned.
%         exitCode A code indicating how the algorithm terminated. Possible
%                  values are
%                  0: Integration was successful.
%                  1: Unable to get a small enough step size.
%                  2: Maximum number of steps reaced without completion.
%                  3: Non-finite number encountered.
%    numRejections The number of times a stepsize hypothesis was rejected.
%
%The basic idea behind adaptive stepsize control for general
%Runge-Kutta-Nyström algorithms is essentially the same as that behind the
%control for standard explicit Runge-Kutta methods. However, the error
%control value used can either be the entire (2N)X1 vector of values and
%their derivatives (position and velocity), or just the NX1 vector of the
%values (position) depending on whether the subsidiary formula chosen
%produces 2N or N values. More information on the stepsize selection
%routine is given in the comments to the function RKAdaptiveOverRange as
%the same method is used here.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<11)
    maxSteps=1024;
end

if(nargin<10||isempty(AbsTol))
    AbsTol=1e-6;
end

if(nargin<9||isempty(RelTol))
   RelTol=1e-3;
end

if(nargin<8||isempty(solutionChoice))
   solutionChoice=0; 
end

if(nargin<7||isempty(order))
    order=5;
end

if(nargin<5||isempty(numRefinedSteps))
    numRefinedSteps=0;
end

if(nargin<4||isempty(probIsGeneral))
    probIsGeneral=true;
end

xDim=size(xVecStart,1);
%The number of components of the state that are not first derivatives. The
%state should contain values and then the derivatives (e.g. position and
%velocity).
xPosDim=xDim/2;

tStart=tSpan(1);
tEnd=tSpan(2);
tDiff=diff(tSpan);
tDiffMag=abs(tDiff);
deltaTSign=sign(tDiff);
%The sign is separated out so that one can use the algorithm running
%backwards as well as forwards in time.

%The maximum step size is arbitraily set to 1/3 the total distance.
maxStepSize=tDiffMag/3;

%Determine the orders of the main and subsidiary embedded Runge-Kutta
%formulae that were chosen. The order of convergence for the step size
%testing is the smallest order taken. The isFSAL flag indicates whether
%g(:,end)' is the value of f evaluated at the next step. The subsidType
%return value indicates the type of subsidiary solution returned, and thus
%affects how the error value for the step size adaptation is computed.
if(probIsGeneral)
    [orders,isFSAL,subsidType]=RungeKNystroemGStep(order,solutionChoice);
else
    [orders,isFSAL,subsidType]=RungeKNystroemSStep(order,solutionChoice);
end
RKOrder=min(orders);

if(nargin<6||isempty(initStepSize))
    %If no initial step size is given, then just use the smallest uniform
    %step size for the given value of maxSteps.
    initStepSize=tDiffMag/maxSteps;
elseif(initStepSize<eps(tStart))
   error('The initial step size provided is too small') 
end

%Allocate the maximum possible space needed for the return parameters.
numSpaces=1+maxSteps*(numRefinedSteps+1);
xVecVals=zeros(xDim,numSpaces);
d2xdt2Vals=zeros(xPosDim,numSpaces);
tVals=zeros(numSpaces,1);

xVecVals(:,1)=xVecStart;
tVals(:,1)=tStart;
d2xdt2Vals(:,1)=df(xVecStart,tStart);
lastStep=1;

deltaTMag=initStepSize;
numRejections=0;
for curStep=2:maxSteps
    xCur=xVecVals(:,lastStep);
    tCur=tVals(lastStep);
    d2xdt2Cur=d2xdt2Vals(:,lastStep);
    
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
    
    %The first time the choice in step size fails, it is adjusted the
    %"optimal" way in the Runge-Kutta-Fehlberg method. Additional times,
    %the step size is just halved in the hope that it will reach an
    %accepted value more quickly.
    failedReducingStepSize=false;
    moveOnToNextStep=false;
    while(moveOnToNextStep==false)
        %Take one step.
        if(probIsGeneral)
            [xPredMain,xPredSubsid,g]=RungeKNystroemGStep(xCur,tCur,df,deltaT,d2xdt2Cur,order,solutionChoice);
        else
            [xPredMain,xPredSubsid,g]=RungeKNystroemSStep(xCur,tCur,df,deltaT,d2xdt2Cur,order,solutionChoice);
        end

        %Integration can only be over finite functions.
        if(any(~isfinite(xPredMain))||any(~isfinite(xPredSubsid)))
            xVecVals=[];
            tVals=[];
            d2xdt2Vals=[];
            exitCode=3;
            return;
        end
        
        %Use the appropriate components when computing the error estimate
        %value.
        if(subsidType==0)%Full state returned as subsidiary
            xPredMainErrEst=xPredMain;
            xCurErrEst=xCur;
        elseif(subsidType==1)%No derivatives returned as sibsidiary
            xPredMainErrEst=xPredMain(1:xPosDim);
            xCurErrEst=xCur(1:xPosDim);
        else%Only derivatives returned as subsidiary
            xPredMainErrEst=xPredMain((xPosDim+1):end);
            xCurErrEst=xCur((xPosDim+1):end);
        end
        
        %The local error estimate. This must be transformed into a 
        %combination relative/ absolute error term to determine whether
        %the step should be rejected.
        normFactor=max(max(abs(xPredMainErrEst),abs(xCurErrEst)),AbsTol/RelTol);
        theError=max(abs((xPredMainErrEst-xPredSubsid)./normFactor));

        if(theError>RelTol)
            %The step should be rejected.
            numRejections=numRejections+1;
            
            if(deltaTMag<deltaTMinMag)
                %If the step size got too small, then return.
                xVecVals=[];
                d2xdt2Vals=[];
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
            %The step is successful, save the results from the step,
            %interpolating as necessary. If refined steps are desired, get
            %the appropriate interpolating polynomial.
            
            if(numRefinedSteps~=0)
                %Get a Hermite interpolating polynomial over the step. 
                [interpPolyA,interpPolyC]=RKNInterpPolys(xCur,tCur,xPredMain,tCur+deltaT,df,probIsGeneral,order,solutionChoice,g);
                %The interpolating polynomials take stapes parameterized by
                %a value between zero and 1 representing the distance from
                %tCur to tCur+deltaT.
                
                deltaFraction=1/(numRefinedSteps+1);
                curStepFrac=deltaFraction;
                for curRStep=1:numRefinedSteps
                    lastStep=lastStep+1;
                    tVals(lastStep)=tVals(lastStep-1)+deltaFraction*deltaT;

                    xVecVals(:,lastStep)=polyValNewton(curStepFrac,interpPolyA,interpPolyC);
                    if(probIsGeneral)
                        d2xdt2Vals(:,lastStep)=df(xVecVals(:,lastStep),tVals(lastStep));
                    else
                        d2xdt2Vals(:,lastStep)=df(xVecVals(1:xPosDim,lastStep),tVals(lastStep));
                    end
                    curStepFrac=curStepFrac+deltaFraction;
                end
            end
            
            lastStep=lastStep+1;
            xVecVals(:,lastStep)=xPredMain;
            tVals(lastStep)=tCur+deltaT;
            
            %Save the current value to be reused on the next step, if the method
            %is an FSAL function, so that an evaluation of f can be avoided.
            if(isFSAL)
                d2xdt2Vals(:,lastStep)=g(:,end);
            else
                if(probIsGeneral)
                    d2xdt2Vals(:,lastStep)=df(xVecVals(:,lastStep),tVals(lastStep));
                else
                    d2xdt2Vals(:,lastStep)=df(xVecVals(1:xPosDim,lastStep),tVals(lastStep));
                end
            end

            %If a step is successful, then increase the step size for the 
            %next step in the standard manner used with Runge-Kutta-
            %Fehlberg methods, but limit the maximum size of the increase
            %to a scale factor of 4. This avoid huge step sizes when the
            %predicted error is very small.
            deltaTMag=min(maxStepSize,deltaTMag*min(4,0.8*(RelTol/theError)^(1/RKOrder)));

            moveOnToNextStep=true;
        end
    end

    if(tVals(curStep)==tEnd)
        %If integration completed, reshape the outputs to match the actual
        %number of steps.
        xVecVals=xVecVals(:,1:curStep);
        tVals=tVals(1:curStep);
        d2xdt2Vals=d2xdt2Vals(:,1:curStep);
        exitCode=0;
        return
    end
end

%If we get here, then the maximum number of iterations was reached without
%reaching the end.
xVecVals=[];
d2xdt2Vals=[];
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
