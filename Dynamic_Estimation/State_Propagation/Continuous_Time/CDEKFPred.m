function [xPred,PPred,exitCode]=CDEKFPred(xPrev,PPrev,a,D,includesJacob,tPrev,tPred,RKOptions)
%%CDEKFPRED Predict forward a Gaussian state estimate through time when the
%           evolution of the state is described by a continuous-time
%           stochastic diferential equation using a continuous-discrete
%           extended Kalman filter. This uses a Gaussian approximation
%           throughout the entire prediction, allowing for a solution based
%           on deterministic differential equations. The differential
%           equations are integrated forward using either the
%           RKAdaptiveOverRange function, an explicit Runge-Kutta method,
%           which can take a number of options, or using a fixed number of
%           Runge-Kutta steps via the RungeKAtTimes function.
%
%INPUTS: xPrev The xDim X 1 state estimate at time tPrior.
%        PPrev The xDim X xDim state covariance matrix at time tPrior.
%            a The drift function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable as its
%              arguments. To avoid the need for numerical differentiation,
%              it is desirable if the second output of the function is the
%              derivative with respect to x of the drift function.
%              Whether or not the function returns this can be indicated
%              using the includesJacob input.
%            D The diffusion function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable as its
%              arguments.
% includesJacob A boolean value indicating whether the drift function a
%              also returns a Jacobian matrix. If this value is false or an
%              empty matrix is passed, it is assumed that a does not return
%              a Jacobian matrix and one will be found using numerical
%              differentiation via the numDiff function with default
%              parameters.
%        tPrev The time of xPrev and SPrev.
%        tPred The time to which xPrev and SPrev should be predicted.
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
%              RKOptions.fixedNumSteps:  (none) If this parameter is
%                                        present, then that number of steps
%                                        of a fixed length (non-adaptive)
%                                        will be taken and the RelTol, 
%                                        AbsTol, maxSteps, and initStepSize
%                                        options are ignored.
%              If an empty matrix is passed in place of RKOptions, then
%              only the default values will be used.
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
%               RKOptions.fixedStepNum option was specified, then this will
%               always be zero. See the RKAdaptiveOverRange function for
%               more details when using an adaptive method.
%
%The algorithm is the mean-covariance Runge-Kutta (MC-RK4) method described
%in [1]. Other methods in this paper are criticized in Section III of [2].
%
%REFERENCES:
%[1] P. Frogerais, J. Bellanger, and L. Senhadji. "Various ways to compute
%    the continuous-discrete extended Kalman filter," IEEE Transactions on
%    Automatic Control, vol. 57, no. 4, pp. 1000-1004, Apr. 2012.
%[2] G. Y. Kulikov and M. V. Kulikova. "Accurate numerical implementation
%    of the continuous-discrete extended Kalman filter," IEEE Transactions
%    on Automatic Control, vol. 59, no. 1, pp. 273-279, Jan. 2014.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=length(xPrev);
    
    if(nargin<8)
        RKOptions=[];
    end
    if(isempty(includesJacob)||includesJacob==false)
        aValJacob=@(x,t)deal(a(x,t),numDiff(x,@(y)a(y,t),xDim));
    else
        aValJacob=a;
    end
    
    %This is only set if a non-adaptive Runge-Kutta algorithm is to be
    %used.
    fixedNumSteps=[];
    
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
        if(isfield(RKOptions,'fixedNumSteps'))
            fixedNumSteps=RKOptions.fixedNumSteps;
        end
    end
    
    if(isempty(initStepSize))
        initStepSize=(tPred-tPrev)/maxSteps;
    end
    
    if(~isempty(fixedNumSteps))
        deltaTMax=(tPred-tPrev)/fixedNumSteps;
        xVals=RungeKAtTimes([xPrev;PPrev(:)],[tPrev tPred],@f,deltaTMax,order,solutionChoice);
        exitCode=0;
    else%If an adaptive routine should be used.
        [xVals,~,~,exitCode]=RKAdaptiveOverRange([xPrev;PPrev(:)],[tPrev tPred],@f,initStepSize,0,order,solutionChoice,RelTol,AbsTol,maxSteps);
    end
    
    if(isempty(xVals))%If the RKAdaptiveOverRange function failed.
        xPred=[];
        PPred=[];
        return;
    end

    xSVec=xVals(:,end);
    xPred=xSVec(1:xDim,end);
    PPred=xSVec((xDim+1):end,end);
    PPred=reshape(PPred,xDim,xDim);
    %Make sure that the matrix remains symmetric.
    PPred=(PPred+PPred')/2;
    
    %This is to deal with the propagation making it possibly not positive
    %definite.
    SPred=cholSemiDef(PPred,'lower');
    PPred=SPred*SPred';
    
    function dVal=f(x,t)
        %This function provides the derivatives of the state and the
        %covariance matrix at a particular time.
        xCur=x(1:xDim);
        PCur=x((xDim+1):end);
        PCur=reshape(PCur,xDim,xDim);
        
        [dx,F]=aValJacob(xCur,t);
        
        dPoints=D(xCur,t);
        dP=PCur*F'+F*PCur'+dPoints*dPoints';
        
        dVal=[dx;dP(:)];
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
