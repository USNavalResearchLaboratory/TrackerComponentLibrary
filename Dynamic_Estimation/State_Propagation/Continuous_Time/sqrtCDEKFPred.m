function [xPred,SPred,exitCode]=sqrtCDEKFPred(xPrev,SPrev,a,D,AJacob,tPrev,tPred,RKOptions)
%%SQRTCDEKFPRED Predict forward a Gaussian state estimate through time when
%               the evolution of the state is described by a
%               continuous-time stochastic diferential equation using a
%               square root continuous-discrete extended Kalman filter.
%               This uses a Gaussian approximation throughout the entire
%               prediction, allowing for a solution based on deterministic
%               differential equations. The differential equations are
%               integrated forward using the RKAdaptiveOverRange function,
%               an explicit Runge-Kutta method, which can take a number of
%               options.
%
%INPUTS: xPrev The xDim X 1 state estimate at time tPrior.
%        SPrev The xDim X xDim square root state covariance matrix at time
%              tPrior.
%            a The drift function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable as its
%              arguments.
%            D The diffusion function in the continuous-time stochastic
%              dynamic model. It takes the state and a time variable as its
%              arguments.
%       AJacob The derivative wrt x of the drift function in the
%              continuous-time stochastic dynamic model. It takes the
%              state and a time variable as its arguments. If an empty
%              matrix is passed, then AJacob will be found using
%              numerical differentiation via the numDiff function with
%              default parameters.
%        tPrev The time of xPrev and SPrev.
%        tPred The time to which xPrev and SPrev should be predicted.
%    RKOptions An optional structure whose components have the same name
%              and meaning as the corresponding components in the
%              RKAdaptiveOverRange function. The possible compnents of
%              the RKOptions function (with default values if omitted)
%              are
%              RKOptions.initStepSize: (tPred-tPrev)/RKOptions.maxSteps
%              RKOptions.order:           5
%              RKOptions.solutionChoice:  0
%              RKOptions.RelTol:          1e-3.
%              RKOptions.AbsTol:          1e-6
%              RKOptions.maxSteps:        1024
%              If an empty matrix is passed in place of RKOptions, then
%              only the default values will be used.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate. If the
%               RKAdaptiveOverRange function failed, which might occur,
%               for example, if the maximum number of allowed steps is
%               too small, then an empty matrix is returned.
%         SPred The xDim X xDim predicted square root state covariance
%               matrix, or an empty matrix if the RKAdaptiveOverRange
%               failed.
%      exitCode The exit code from the RKAdaptiveOverRange function. This
%               is zero if the integration was a success. Otherwise, the
%               value identifies what the problem was. See the
%               RKAdaptiveOverRange function for more details.
%
%The algorithm is modified from the mean-covariance Runge-Kutta (MC-RK4)
%method described in [1]. Other methods in this paper are criticized in
%Section III of [2].
%
%Parts of the calculation for the derivative of the square root state
%covariance matrix are taken from Section VI of [3].
%
%REFERENCES:
%[1] P. Frogerais., J. Bellanger, and L. Senhadji. "Various ways to compute
%    the continuous-discrete extended Kalman filter," IEEE
%    Transactions on Automatic Control, vol. 57, no. 4, pp. 1000-1004,
%    2012.
%[2] G. Y. Kulikov and M. V. Kulikova. "Accurate numerical implementation
%    of the continuous-discrete extended Kalman filter," IEEE Transactions
%    on Automatic Control, vol. 59, no. 1, pp. 273-279, Jan. 2014.
%[3] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, Part II, pp. 4-41, Feb. 2015.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=length(xPrev);
    
    if(nargin<8)
        RKOptions=[];
    end
    if(isempty(AJacob))
        AJacob=@(x,t)numDiff(x,@(y) a(y,t),xDim);
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
    end
    
    if(isempty(initStepSize))
        initStepSize=(tPred-tPrev)/maxSteps;
    end
    
    %The zeros in S are not given to the differential equation solver.
    %lowerTriEls holds the indices of the lower-triangular elements.
    lowerTriEls=tril(true(xDim,xDim));
    [xVals,~,~,exitCode]=RKAdaptiveOverRange([xPrev;SPrev(lowerTriEls)],[tPrev tPred],@f,initStepSize,0,order,solutionChoice,RelTol,AbsTol,maxSteps);
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
        PCur=SCur*SCur';
        
        dx=a(xCur,t);
        
        dPoints=D(xCur,t);
        F=AJacob(xCur,t);
        dP=PCur*F'+F*PCur'+dPoints*dPoints';
        SInv=pinv(SCur);
        dS=SCur*triLower(SInv*dP*SInv');
        
        dVal=[dx;dS(lowerTriEls)];
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
