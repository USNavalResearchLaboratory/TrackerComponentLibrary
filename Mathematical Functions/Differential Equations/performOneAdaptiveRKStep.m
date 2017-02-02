function [deltaTNew,xNew,tNew,k,dxdtCur,exitCode]=performOneAdaptiveRKStep(xCur,tCur,f,deltaT,deltaTMinMag,deltaTMaxMag,dxdtCur,order,solutionChoice,AbsTol,RelTol)
%%PERFORMONEADAPTIVERKSTEP Perform a single adaptive Runge-Kutta step,
%                   returning the new adjusted stepsize. This function is a
%                   subroutine of the function RKAdaptiveOverRange. Though
%                   most folks will just directly use the function
%                   RKAdaptiveOverRange to perform adaptive Runge-Kutta
%                   integration over a particular timespan, this subroutine
%                   has been broken out separately as one might want to
%                   program Runge-Kutta integration routines that
%                   adaptively integrate until a certain criterion is
%                   satisfied. As this function is meant as an efficient
%                   subroutine, none of the inputs provide default
%                   parameters if omitted.
%
%INPUTS:   xCur The NX1 state vector at time tCur.
%          tCur The scalar current time of integration.
%             f The function handle for f(x,t)=dxdt over which integration
%               is to be performed. The output is NX1-dimensional.
%        deltaT The current stepsize to be taken in t. This can be positive
%               or negative.
%  deltaTMinMag The minimum allowable magnitude of the step size.
%  deltaTMaxMag The maximum allowable magnitude of the step size.
%      dxdtCur  The value f(xCur,tCur). This is requested so that
%               methods that are FSAL (See comments to the function
%               RungeKStep) can avoid additional computations.
% order,solutionChoice  A pair of optional parameters that specify the
%                highest order of the embedded Runge-Kutta pair to use as
%                well as the specific algorithm to use. Details are given
%                in the comments to the RungeKStep function. If omitted or
%                empty matrices are passed, the default order of 5 is used
%                and the default solutionChoice of 0 is used.
%         RelTol The maximum relative error tolerance allowed, a
%                positive scalar (its use is explained in more detail
%                below).
%         AbsTol The absolute error tolerance allowed, a positive scalar,
%                of the same for all components of x, or a positive NX1
%                vector.
%
%OUTPUTS: deltaTNew The value of deltaT that can be used for the next
%               adaptive step (i.e. pass it to this function). If the
%               function fails (exitCode~=0), then this will be set to the
%               last step size used.
%         xNew  The updated NX1 state vector.
%         tNew  The time of the updated state vector.
%            k  The values of the derivatives f evaluated at various
%               points as determined by the algorithm used for the
%               step. This can be passed to functions like RKInterpPolys to
%               perform interpolation.
%      exitCode A code indicating whether the step could be successfully
%               performed. Possible values are
%               0: Integration was successful.
%               1: Unable to get a small enough step size.
%               3: Non-finite number encountered.
%
%The algorithm used for the step is described in the comments to the
%function RKAdaptiveOverRange.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Determine the orders of the main and subsidiary embedded Runge-Kutta
    %formulae that were chosen. The order of convergence for the step size
    %testing is the smallest order taken. The isFSAL flag indicates whether
    %k(:,end) is the value of f evaluated at the next step.
    [orders,isFSAL]=RungeKStep(order,solutionChoice);
    RKOrder=min(orders);

    deltaTMag=abs(deltaT);
    deltaTSign=sign(deltaT);
    
    %The first time the choice in step size fails, it is adjusted the
    %"optimal" way in the Runge-Kutta-Fehlberg method. Additional times,
    %the step size is just halved in the hope that it will reach an
    %accepted value more quickly.
    failedReducingStepSize=false;
    moveOnToNextStep=false;
    while(moveOnToNextStep==false)
        %The returned value dxValdt in k(:,1) is for xCur and tCur, which
        %is curStep-1.
        [xPredMain,xPredSubsid,k]=RungeKStep(xCur,tCur,f,deltaT,dxdtCur,order,solutionChoice);

        %Integration can only be over finite functions.
        if(any(~isfinite(xPredMain))||any(~isfinite(xPredSubsid)))
            xNew=[];
            tNew=[];
            exitCode=3;
            deltaTNew=deltaTMag;
            return;
        end
        
        %The local error estimate. This must be transformed into a 
        %combination relative/ absolute error term to determine whether
        %the step should be rejected.
        normFactor=max(max(abs(xPredMain),abs(xCur)),AbsTol/RelTol);
        theError=max(abs((xPredMain-xPredSubsid)./normFactor));

        if(theError>RelTol)
            if(deltaTMag<deltaTMinMag)
                %If the step size got too small, then return.
                xNew=[];
                tNew=[];
                exitCode=1;
                deltaTNew=[];
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
            %The step was successful; leave the loop.
            moveOnToNextStep=true;
        end
    end
    
    %The step is successful, save the results from the step, including
    %information so that interpolation can be perfromed, if necessary.
    xNew=xPredMain;
    tNew=tCur+deltaT;

    %Save the current value to be reused on the next step, if the method
    %is an FSAL function, so that an evaluation of f can be avoided.
    if(isFSAL)
        dxdtCur=k(:,end);
    else
        dxdtCur=f(xNew,tNew);
    end

    %If a step is successful, then increase the step size for the 
    %next step in the standard manner used with Runge-Kutta-
    %Fehlberg methods, but limit the maximum size of the increase
    %to a scale factor of 4. This avoid huge step sizes when the
    %predicted error is very small.
    deltaTMag=min(deltaTMaxMag,deltaTMag*min(4,0.8*(RelTol/theError)^(1/RKOrder)));

    %Since the minimum step size changes every loop, this makes sure that
    %deltaT does not go beneath it just because the loop changed.
    deltaTMag=max(deltaTMag,deltaTMinMag);
    deltaTNew=deltaTMag*deltaTSign;
            
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
