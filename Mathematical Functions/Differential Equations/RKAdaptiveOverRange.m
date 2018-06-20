function [xVals,tVals,dxdtVals,exitCode,nextStepSize]=RKAdaptiveOverRange(xStart,tSpan,f,initStepSize,numRefinedSteps,order,solutionChoice,RelTol,AbsTol,maxSteps)
%%RKADAPTIVEOVERRANGE Perform explicit Runge-Kutta integration over a
%                     given range of values using an adaptive step size.
%                     Runge-Kutta methods are derivative-free techniques
%                     for solving ordinary differential equations. That is,
%                     integrating dx/dt=f(x,t) given initial conditions
%                     (xStart,tSpan(1)). More information on available
%                     algorithms for the steps is given in the comments to
%                     the function RungeKStep.
%
%INPUTS: xStart The NX1 state vector at time tSpan(1).
%         tSpan A 2X1 or 1X2 vector where tSpan(1) is the starting time and
%               tSpan2 is the desired stopping time for the integration.
%             f The function handle for f(x,t)=dxdt over which integration
%               is to be performed. The output is NX1-dimensional.
%  initStepSize An optional initial step size (in t) to use for the
%               integration. If omitted or an empty matrix is passed, an
%               ad-hoc method described below is used to find an initial
%               step size.
% numRefinedSteps An optional parameter specifying the number of steps to
%               interpolate between each normal step. With high-order
%               formulas, the steps produced might not be very nice to
%               use for plotting, so setting this to a value above 0 adds
%               more intermediate values (using the polynomial returned by
%               the RKInterpPolys function) so that plotted data can look
%               nicer. If omitted or an empty matrix is passed, the default
%               value of 0 (no extra/ interpolated steps) is used.
%               Interpolation is done at the main order of the Runge-Kutta
%               method.
% order,solutionChoice  A pair of optional parameters that specify the
%               highest order of the embedded Runge-Kutta pair to use as
%               well as the specific algorithm to use. Details are given in
%               the comments to the RungeKStep function. If omitted or
%               empty matrices are passed, the default order of 5 is used
%               and the default solutionChoice of 0 is used.
%        RelTol The maximum relative error tolerance allowed, a
%               positive scalar (its use is explained in more detail
%               below). If omitted or an empty matrix is passed, the
%               default value of 1e-3 is used.
%        AbsTol The absolute error tolerance allowed, a positive scalar, of
%               the same for all components of x, or a positive NX1 vector.
%               If omitted or an empty matrix is passed, the default value
%               of 1e-6 is used.
%      maxSteps The maximum allowable number of steps to perform the
%               integration. If omitted, the default of 1024 is used.
%
%OUTPUTS: xVals The NXnumSteps set of values of x along the path.
%               xVals(:,1) is xStart and xVals(:,end) is the value at final
%               time tSpan(2). If the Runge-Kutta integration failed, e.g.
%               due to encountering a NaN or being unable to get a
%               sufficient step size, an empty matrix is returned.
%         tVals A numStepsX1 vector of values of t corresponding to the
%               values of x in xVals. tVals(1) is equal to tSpan(1) and
%               tVals(end) is equal to tSpan(2). If integration fails, an
%               empty matrix is returned.
%      dxdtVals A numStepsXN array containing derivatives of x evaluated at
%               the times in tVals. The derivatives are just the result of
%               evaluating f(x,t) at each point in (xVals,tVals). This
%               combined with xVals and tVals could be used in functions to
%               perform Hermite interpolation. If integration fails, an
%               empty matrix is returned.
%      exitCode A code indicating how the algorithm terminated. Possible
%               values are
%               0: Integration was successful.
%               1: Unable to get a small enough step size.
%               2: Maximum number of steps reached without completion.
%               3: Non-finite number encountered.
%  nextStepSize The step size (in t) that would be used for the next
%               integration step. If empty, then the integration ended
%               because the step size was too small.
%
%The basic idea behind adaptive stepsize control can be ascertained from
%Chapter 5.2 of [3]. However, the implementation here has a number of
%changes. For example, unlike the basic algorithm of [3], depending on the
%user's choice of order and solutionChoice, the main Runge-Kutta formula,
%which is the one whose result is saved, might have a higher or lower order
%than the subsidiary formula, which is only used to check the error to
%determine how the stepsize should change. The smallest order of the two is
%the order that is used in formulae for adaptively changing the step size.
%
%The algorithm needs an initial stepsize to start. If the user does not
%provide an intial stepsize, then the initial stepsize is just set to the
%width of the integration region in t divided by the maximum number of
%steps. That is, initStepSize = (tSpan(2)-tSpan(1))/maxSteps;
%
%Once an initial stepsize has been determined, the function RungeKStep is
%used to perform a single step of the Runge-Kutta algorithm. Then one must
%determine whether the stepsize was small enough. If the stepsize is small
%enough, then the stepsize should still be adaptively changed for the nest
%step. The procedure for updating the stepsize is essentially the procedure
%for non-stiff equations from [2], Equation 2.3). The method has been
%slightly modified in that the error being adjusted is not just the
%absolute error, as explained below. Also, the maximum stepsize is limited
%to 1/5 the integration region and increases in the stepsize are limited
%so as not to be too large.
%
%We would like the error in each step to be less than the absolute and the
%relative error tolerances. For an absolute error tolerance, this means
%that
%abs(xPredMain-xPredSubsid)<AbsTol
%for all components of the predicted value xPred in a particular step of
%size deltaT. For a relative error tolerance, this means that
%abs(xPredMain-xPredSubsid)./max(max(abs(xPred),abs(xCur)),RelTol)<RelTol
%for all components. The extra RelTol in the denominator on the left helps
%deal with issues when the true value is close to zero.
%
%The requirements for absolute and relative tolerances can be combined such
%that the absolute tolerance only comes into play when the numbers are
%sufficiently small that the relative tolerance becomes problematic. The
%two lines
% normFactor=max(max(abs(xPred),abs(xHatPred)),AbsTol/RelTol);
% theError=max(abs((xPred-xHatPred)./normFactor));
%Produce a single quantity, theError, which, if less than RelTol, means
%that the stepsize was too large. If the denominator of the normFactor is
%too small, then it is just AbsTol/RelTol, in which case the RelTol values
%on both sides of the equal sign cancel. And the comparison is just the
%comparison for the absolute tolerance (looking only at the largest
%error component). Thus theError term uses a threshold factor of
%AbsTol/RelTol. This is used as a similar threshold ratio when creating an
%ad-hoc initial estimate.
%
%REFERENCES:
%[1] I. Gladwell, L. F. Shampine, and R. W. Brankin, "Automatic selection
%    of the initial step size for an ODE solver," Journal of Computational
%    and Applied Mathematics, vol. 18, no. 2, pp. 175-192, May 1987.
%[2] G. Hall, "A new stepsize strategy for explicit Runge-Kutta codes,"
%    Advances in Computational Mathematics, vol. 3, no. 4, pp. 343-352, May
%    1995.
%[3] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/ Cole, 2011.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10)
    maxSteps=1024;
end

if(nargin<9||isempty(AbsTol))
    AbsTol=1e-6;
end

if(nargin<8||isempty(RelTol))
   RelTol=1e-3;
end

if(nargin<7||isempty(solutionChoice))
   solutionChoice=0; 
end

if(nargin<6||isempty(order))
    order=5;
end

if(nargin<5||isempty(numRefinedSteps))
    numRefinedSteps=0;
end

xDim=size(xStart,1);
tStart=tSpan(1);
tEnd=tSpan(2);
tDiff=diff(tSpan);
tDiffMag=abs(tDiff);
deltaTSign=sign(tDiff);
%The sign is separated out so that one can use the algorithm running
%backwards as well as forwards in time.

%The maximum step size is arbitrarily set to 1/3 the total distance.
deltaTMaxMag=tDiffMag/3;

if(nargin<4||isempty(initStepSize))
    %If no initial step size is given, then just use the smallest uniform
    %step size for the given value of maxSteps.
    initStepSize=tDiffMag/maxSteps;
elseif(initStepSize<eps(tStart))
   error('The initial step size provided is too small.') 
end

%Allocate the maximum possible space needed for the return parameters.
numSpaces=1+maxSteps*(numRefinedSteps+1);
xVals=zeros(xDim,numSpaces);
dxdtVals=zeros(xDim,numSpaces);
tVals=zeros(numSpaces,1);

xVals(:,1)=xStart;
tVals(1)=tStart;
dxdtVals(:,1)=f(xStart,tStart);
lastStep=1;

deltaT=deltaTSign*initStepSize;
for curStep=2:(maxSteps+1)
    xCur=xVals(:,lastStep);
    tCur=tVals(lastStep);
    dxdtCur=dxdtVals(:,lastStep);
    
    %The absolute value of the minimum allowable step size. It is set so
    %that the step must makes something of a difference compared to the
    %numerical precision.
    deltaTMinMag=2^4*eps(tCur);

    %If we would have overstepped the end, change deltaT to be the end.
    %Then, we can terminate at the end.
    if(deltaTSign*(tCur+deltaT)>tEnd*deltaTSign)
        deltaTMag=abs(tEnd-tCur);
        deltaT=deltaTMag*deltaTSign;
        %Allow the last step to be very small.
        deltaTMinMag=min(deltaTMinMag,deltaTMag);
    end
    
    [deltaTNew,xNew,tNew,k,dxdtCur,exitCode]=performOneAdaptiveRKStep(xCur,tCur,f,deltaT,deltaTMinMag,deltaTMaxMag,dxdtCur,order,solutionChoice,AbsTol,RelTol);

    %If an error occurred, then just return.
    if(exitCode~=0)
        xVals=[];
        tVals=[];
        dxdtVals=[];
        nextStepSize=abs(deltaTNew);
        return
    end
    
    dxdtVals(:,lastStep)=dxdtCur;

    if(numRefinedSteps~=0)
        %If extra interpolation between the steps is to be performed.
        %Get a Hermite interpolating polynomial over the step. 
        [interpPolyA,interpPolyC]=RKInterpPolys(xCur,tCur,xNew,tNew,f,order,solutionChoice,k);
        %The interpolating polynomials take stapes parameterized by
        %a value between zero and 1 representing the distance from
        %tCur to tNew.
        deltaTTaken=tNew-tCur;
        
        deltaFraction=1/(numRefinedSteps+1);
        curStepFrac=deltaFraction;
        for curRStep=1:numRefinedSteps
            lastStep=lastStep+1;
            tVals(lastStep)=tVals(lastStep-1)+deltaFraction*deltaTTaken;

            xVals(:,lastStep)=polyValNewton(curStepFrac,interpPolyA,interpPolyC);
            dxdtVals(:,lastStep)=f(xVals(:,lastStep),tVals(lastStep));
            curStepFrac=curStepFrac+deltaFraction;
        end
    end
    
    lastStep=lastStep+1;
    xVals(:,lastStep)=xNew;
    tVals(lastStep)=tNew;
    dxdtVals(:,lastStep)=dxdtCur;
    deltaT=deltaTNew;

    if(tVals(lastStep)==tEnd)
        %If integration completed, reshape the outputs to match the actual
        %number of steps.
        xVals=xVals(:,1:lastStep);
        tVals=tVals(1:lastStep);
        dxdtVals=dxdtVals(:,1:lastStep);
        exitCode=0;
        nextStepSize=deltaTMag;
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
