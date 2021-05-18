function [xPred,SPred,exitCode]=momentMatchPred(xPrev,SPrev,a,D,dColDim,tPrev,tPred,xi,w,RKOptions,stateDiffFunc)
%%MOMENTMATCHPRED Predict forward a Gaussian state estimate through time
%                 when the evolution of the state is described by a
%                 continuous-time stochastic diferential equation. This
%                 uses a Gaussian approximation throughout the entire
%                 prediction, allowing for a solution based on
%                 deterministic differential equations. The differential
%                 equations are integrated forward using either the
%                 RKAdaptiveOverRange function, an explicit Runge-Kutta
%                 method, which can take a number of options, or using a
%                 fixed number of Runge-Kutta steps via the RungeKAtTimes
%                 function.
%
%INPUTS: xPrev The xDim X 1 state estimate at time tPrior.
%        SPrev The xDim X xDim lower-triangular square root of the state
%              covariance matrix at time tPrior.
%            a The drift function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable as its
%              arguments.
%            D The diffusion function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable
%              as its arguments.
%      dColDim The number of columns that a matrix obtained by executing
%              D(x,t) will have. The number of rows must be the
%              dimensionality of the state.
%        tPrev The time of xPrev and SPrev.
%        tPred The time to which xPrev and SPrev should be predicted.
%           xi An xDim X numCubPoints matrix of cubature points.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points.
%    RKOptions An optional structure whose components have the same name
%              and meaning as the corresponding components in the
%              RKAdaptiveOverRange function with the exception of the
%              fixedStepNum option, in which case a Runge-Kutta method
%              using a fixed step size is used in place of the adaptive
%              method. The possible components of the RKOptions function
%              (with default values if omitted) are
%              RKOptions.initStepSize: (tPred-tPrev)/RKOptions.maxSteps
%              RKOptions.order:           5
%              RKOptions.solutionChoice:  0
%              RKOptions.RelTol:          1e-3.
%              RKOptions.AbsTol:          1e-6
%              RKOptions.maxSteps:        1024
%              RKOptions.fixedStepNum:   (none) If this parameter is
%                                        present, then that number of steps
%                                        of a fixed length (non-adaptive)
%                                        will be taken and the RelTol, 
%                                        AbsTol, maxSteps, and initStepSize
%                                        options are ignored.
%              If an empty matrix is passed in place of RKOptions, then
%              only the default values will be used.
% stateDiffFunc An optional function handle that transforms the value of
%              the difference between two state estimates. This must be
%              able to handle sets of differences between state estimates.
%              For an xDim state, this must handle an xDimXN  matrix of N
%              differences. This only needs to be supplied if the
%              difference has to be restricted to a certain range. For
%              example, if elements of the state are angles, then the
%              difference between angles in two states should be limited to
%              +/-pi using the wrapRange function for that element.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate. If the 
%               RKAdaptiveOverRange function failed, which might occur, for
%               example, if the maximum number of allowed steps is too
%               small, then an empty matrix is returned.
%         SPred The xDim X xDim predicted state covariance matrix, or an
%               empty matrix if the RKAdaptiveOverRange failed.
%      exitCode The exit code from the RKAdaptiveOverRange function. This
%               is zero if the integration was a success. Otherwise, the
%               value identifies what the problem was. If the
%               RKOptions.fixedStepNum option was specified, then this
%               will always be zero. See the RKAdaptiveOverRange function
%               for more details when using an adaptive method.
%
%The algorithm is the moment matching algorithm given in Section VI of
%[1]. The optional parameter stateDiffFunc is not described in the above
%reference, but has been added to allow the filter to be used
%with states containing angular quantities. For example, if the state
%consists of 2D position, heading (in radians) and speed, then 
%stateDiffFunc=@(xDiff)[xDiff(1:2,:);
%                       wrapRange(xDiff(3,:),-pi,pi);
%                       xDiff(4,:)];
%should be used.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, Part II, pp. 4-41, Feb. 2015.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numCubPoints=length(w);
    xDim=length(xPrev);
    sqrtW=sqrt(w);
    
    %This is only set if a non-adaptive Runge-Kutta algorithm is to be
    %used.
    fixedNumSteps=[];
    
    if(nargin<10)
        RKOptions=[];
    end
    
    %The default options for the RKAdaptiveOverRange function. The
    %initStepSize default option is not set until the other options are
    %read, since the default depends on the value of maxSteps.
    initStepSize=[];
    order=5;
    solutionChoice=0;
    RelTol=1e-3;
    AbsTol=1e-6;
    maxSteps=1024;
    if(~isempty(RKOptions))  
        if(isfield(RKOptions,'initStepSize'))
            initStepSize=RKOptions.initStepSize;
        end
        if(isfield(RKOptions,'order'))
            order=RKOptions.order;
        end
        if(isfield(RKOptions,'solutionChoice'))
            solutionChoice=RKOptions.solutionChoice;
        end
        if(isfield(RKOptions,'RelTol'))
            RelTol=RKOptions.RelTol;
        end
        if(isfield(RKOptions,'AbsTol'))
            AbsTol=RKOptions.AbsTol;
        end
        if(isfield(RKOptions,'maxSteps'))
            maxSteps=RKOptions.maxSteps;
        end
        if(isfield(RKOptions,'fixedStepNum'))
            fixedNumSteps=RKOptions.fixedStepNum;
        end
    end
    
    if(isempty(initStepSize))
        initStepSize=(tPred-tPrev)/maxSteps;
    end
    
    %If no function is given to wrap the state to a particular range of
    %values.
    if(nargin<11)
        stateDiffFunc=@(x)x;
    end
    
    %The zeros in S are not given to the differential equation solver.
    %lowerTriEls holds the indices of the lower-triangular elements.
    lowerTriEls=tril(true(xDim,xDim));
    
    if(~isempty(fixedNumSteps))
        deltaTMax=(tPred-tPrev)/fixedNumSteps;
        xVals=RungeKAtTimes([xPrev;SPrev(lowerTriEls)],[tPrev tPred],@f,deltaTMax,order,solutionChoice);
        exitCode=0;
    else%If an adaptive routine should be used.
        [xVals,~,~,exitCode]=RKAdaptiveOverRange([xPrev;SPrev(lowerTriEls)],[tPrev tPred],@f,initStepSize,0,order,solutionChoice,RelTol,AbsTol,maxSteps);
    end
    if(isempty(xVals))%If the RKAdaptiveOverRange function failed.
        xPred=[];
        SPred=[];
        return;
    end

    xSVec=xVals(:,end);
    xPred=xSVec(1:xDim,end);
    SEls=xSVec((xDim+1):end,end);
    SPred=zeros(xDim,xDim);
    SPred(lowerTriEls)=SEls;
    function dVal=f(x,t)
        %This function provides the derivatives of the state and the
        %lower-triangular square root covariance matrix at a particular
        %time.
        xCur=x(1:xDim);
        SElsCur=x((xDim+1):end);
        %Only the nonzero elements of S were passed.
        SCur=zeros(xDim,xDim);
        SCur(lowerTriEls)=SElsCur;
        
        %Cubature State Points
        XpMat=transformCubPoints(xi,xCur,SCur);
        
        %Fractional modified cubature points.
        XPMatTilde=zeros(xDim,numCubPoints);
        for curPoint=1:numCubPoints
            XPMatTilde(:,curPoint)=a(XpMat(:,curPoint),t);
        end
        
        %The function XPMatTilde*w is the mean.
        dx=XPMatTilde*w;
        
        dPoints=D(XpMat,t);
        DCat=reshape(bsxfun(@times,dPoints,reshape(sqrtW,1,1,numCubPoints)),xDim,dColDim*numCubPoints);

        dP=dPSum(XpMat,XPMatTilde,DCat,xCur,dx,sqrtW,stateDiffFunc);
        %If a numerical problem arose, then flag it here so that the svd
        %function inside of pinv does not have an error.
        if(any(~isfinite(SCur)))
           dVal=NaN*ones(size(x));
           return;
        end
        
        SCurInv=pinv(SCur);
        dS=SCur*triLower(SCurInv*dP*SCurInv');
        dVal=[dx;dS(lowerTriEls)];
    end
end

function dP=dPSum(XpMat,XPMatTilde,DCat,xCur,dx,sqrtW,stateDiffFunc)
%%DPSUM Compute dP using sums rather than the matrix form given in Equation
%       (75) of the paper. The sum form follows from the expected values in
%       Equation 73c. This allows functions to transform differences in
%       state and differences in the derivative of the state to be used.
%       For example, if the state has angular components, then angular
%       differences should be wrapped using wrapRange.
%
%INPUTS: xPMat The matrix of state cubature points.
%   XPMatTilde The matrix of state derivative cubature points.
%         DCat The matrix of diffusion matrix cubature points.
%         xCur The mean state estimate.
%           dx The mean state derivative estimate.
%        sqrtW The square root cubature point vector.
% stateDiffFunc The user-provided function to transform differences of
%              state estimates.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDiff=bsxfun(@times,stateDiffFunc(bsxfun(@minus,XpMat,xCur)),sqrtW');
    dxDiff=bsxfun(@times,bsxfun(@minus,XPMatTilde,dx),sqrtW');
    
    dP=xDiff*dxDiff'+dxDiff*xDiff'+DCat*DCat';
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
