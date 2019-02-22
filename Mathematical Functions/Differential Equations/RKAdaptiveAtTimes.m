function xList=RKAdaptiveAtTimes(xInit,theTimes,f,initStepSize,order,solutionChoice,RelTol,AbsTol,maxSteps)
%%RKADAPTIVEATTIMES Perform multiple steps of Runge-Kutta propagation
%                   using an adaptive step size. Runge-Kutta methods are
%                   derivative-free techniques for solving ordinary
%                   differential equations. That is, integrating
%                   dx/dt=f(x,t) given initial conditions
%                   (xStart,tSpan(1)). More information on available
%                   algorithms for the steps is given in the comments to
%                   the function RungeKStep.
%
%INPUTS: xInit The initial value of the state (scalar or vector) over which
%              integration is being performed.
%     theTimes The times at which state estimates are desired. theTimes(1)
%              is the time of xInit.
%            f f(xVal,curT) returns the derivative of xVal taken at time
%              curT.
% initStepSize An optional initial step size (in t) to use for the
%              integration. If omitted or an empty matrix is passed, an
%              ad-hoc method is used to find an initial step size.
%        order The order of the Runge-Kutta method. If this parameter is
%              omitted, then the default order of 4 is used. Order can
%              range from 1 to 7.
%  solutionChoice When multiple formulae are implemented, this selects
%              which one to use. Otherwise, this parameter is not used.
%       RelTol The maximum relative error tolerance allowed, a positive
%              scalar. If omitted or an empty matrix is passed, the default
%              value of 1e-3 is used.
%       AbsTol The absolute error tolerance allowed, a positive scalar, of
%              the same for all components of x, or a positive NX1 vector.
%              If omitted or an empty matrix is passed, the default value
%              of 1e-6 is used.
%     maxSteps The maximum allowable number of steps to perform the
%              integration. If omitted, the default of 1024 is used.
%
%OUTPUTS: xList The state at the times given in times xList(:,1) is the
%               same as xInit.
%
%A detailed description of the adaptive step size algorithm can be found in
%the comments of RKAdaptiveOverRange.
%
%May 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(maxSteps))
    maxSteps=1024;
end

if(nargin<8||isempty(AbsTol))
    AbsTol=1e-6;
end

if(nargin<7||isempty(RelTol))
    RelTol=1e-3;
end
if(nargin<6||isempty(solutionChoice))
    solutionChoice=0;
end

if(nargin<5||isempty(order))
    order=5;
end

if(nargin<4||isempty(initStepSize))
    initStepSize=[];
end

xDim=length(xInit);
numTimes=length(theTimes);

xList=zeros(xDim,numTimes);
xList(:,1)=xInit;
for curTime=2:numTimes
    tSpan=theTimes(curTime-1:curTime);
    [xVals,~,~,~,~,initStepSize]=RKAdaptiveOverRange(xList(:,curTime-1),tSpan,f,initStepSize,0,order,solutionChoice,RelTol,AbsTol,maxSteps);
    xList(:,curTime)=xVals(:,end);
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
