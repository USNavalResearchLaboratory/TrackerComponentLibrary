function xList=RungeKSteps(xInit,t0,f,deltaT,numSteps,order,solutionChoice)
%%RUNGEKSTEPS Perform multiple steps of Runge-Kutta propagation with a
%             fixed time interval between steps.
%
%INPUTS:    xInit    The initial value of the state (scalar or vector) over
%                    which integration is being performed.
%           t0       The time at which xInit is taken.
%           f        f(xVal,curT) returns the derivative of xVal taken at
%                    time curT.
%           deltaT   The size of the steps in the Runge-Kutta integration.
%           numSteps The number of steps over which Runge-Kutta integration
%                    is performed.
%           order   The order of the Runge-Kutta method. If this parameter
%                   is omitted, then the default order of 4 is used. Order
%                   can range from 1 to 7.
%   solutionChoice  When multiple formulae are implemented, this selects
%                   which one to use. Otherwise, this parameter is not
%                   used. Currently, only the order=4 method has multiple
%                   solutions implemented in which case omitting this
%                   parameter or setting it to zero used the Dormand and
%                   Prince Algorithm, and setting it to 1 uses the Fehlberg
%                   algorithm.
%
%OUTPUTS:   xList   The first column of xList is xInit. The subsequent
%                   columns are the state propagated forward at intervals
%                   of deltaT.
%
%This function calls the RungeKStep function to propagate forward the
%target state each step of duration deltaT. See the comments in RungeKStep
%for more information.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<7)
        solutionChoice=0;
    end
    
    if(nargin<6)
        order=4;
    end

    xDim=size(xInit,1);
    xList=zeros(xDim,numSteps+1);
    xList(:,1)=xInit;

    curT=t0;
    for curStep=1:numSteps    
        xList(:,curStep+1)=RungeKStep(xList(:,curStep),curT,f,deltaT,[],order,solutionChoice);
        curT=curT+deltaT;
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
