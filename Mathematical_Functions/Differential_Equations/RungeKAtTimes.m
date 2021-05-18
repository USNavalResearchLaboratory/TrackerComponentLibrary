function xList=RungeKAtTimes(xInit,theTimes,f,deltaTMax,order,solutionChoice)
%%RUNGEKATTIMES Perform multiple steps of Runge-Kutta propagation with a
%               possibly variable time-interval between steps. This 
%               accounts for the possibility that the maximum step size
%               allowed in the Runge-Kutta method might not line up with
%               the times at which the state is desired.
%
%INPUTS: xInit The initial value of the state (scalar or vector) over which
%              integration is being performed.
%     theTimes The times at which state estimates are desired. theTimes(1)
%              is the time of xInit.
%            f f(xVal,curT) returns the derivative of xVal taken at time
%              curT.
%    deltaTMax The maximum allowable step size in the Runge-Kutta
%              integration between times in theTimes. If this parameter is
%              omitted, it is just set to Inf, meaning that the steps will
%              be at the actual times.
%        order The order of the Runge-Kutta method. If this parameter is
%              omitted, then the default order of 4 is used. Order can
%              range from 1 to 7.
% solutionChoice When multiple formulae are implemented, this selects which
%              one to use. Otherwise, this parameter is not used.
%
%OUTPUTS: xList The state at the times given in times xList(:,1) is the
%               same as xInit.
%
%This technique tries to propagate the state to the desired times taking
%the largest uniform stepsize between samples possible that does not exceed
%deltaTMax. The function RungeKSteps is used to perform the propagation
%for the required number of steps between samples.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<6)
       solutionChoice=0; 
    end

    if(nargin<5)
        order=5;
    end

    if(nargin<4)
        deltaTMax=Inf;
    end

    xDim=length(xInit);
    numTimes=length(theTimes);

    xList=zeros(xDim,numTimes);
    xList(:,1)=xInit;
    for curTime=2:numTimes
        deltaTStep=theTimes(curTime)-theTimes(curTime-1);
        if(deltaTStep<deltaTMax)
            deltaT=deltaTStep;
            numSteps=1;
        else
            numSteps=ceil(deltaTStep/deltaTMax);
            deltaT=deltaTStep/numSteps;
        end
        
        xList(:,curTime)=RungeKSteps(xList(:,curTime-1),theTimes(curTime-1),f,deltaT,numSteps,order,solutionChoice,true);        
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
