function [xPred, SPred, xPropCenPoints]=sqrtDiscCubKalPred(xPrev,SPrev,f,SQ,xi,w,stateDiffTrans,stateAvgFun,stateTrans)
%%SQRTDISCCUBKALPRED Perform the discrete-time prediction step that comes  
%                    with the square root implementation of the cubature
%                    Kalman filter with additive process noise.
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        SPrev The xDim X xDim lower-triangular square root of the state
%              covariance matrix at the previous time-step.
%            f A function handle for the state transition function that
%              takes the state as its parameter.
%           SQ The xDimXxDim lower-triangular square root of the process
%              noise covariance matrix.
%           xi An xDimXnumCubPoints matrix of cubature points. If this and
%              the next parameter are omitted or empty matrices are passed,
%              then fifthOrderCubPoints(xDim+cDim) is used. It is suggested
%              that xi and w be provided to avoid needless recomputation of
%              the cubature points.
%            w A numCubPointsX1 vector of the weights associated with the
%              cubature points. These must be all positive.
% stateDiffTrans An optional function handle that takes an xDimXN matrix of
%              N differences between states estimates and transforms them
%              however might be necessary. For example, a state continaing
%              angular components will generally need differences between
%              angular components wrapped to the range +/-pi.
%  stateAvgFun An optional function that given an xDimXN matrix of N state
%              estimates and an NX1 vector of weights, provides the
%              weighted average of the state estimates. This is necessary
%              if, for example, states with angular components are
%              averaged.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one should generally
%              want to bind it to the region +/-pi. This is not applied to
%              the output of f.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate.
%         SPred The xDim X xDim lower-triangular square root of the
%               predicted state covariance estimate.
% xPropCenPoints The centered propagated cubature state points. This matrix
%               is needed if the backwards smoothing step is being applied.
%
%The mathematics behind the function sqrtDiscCubKalPred are described in 
%more detail in Section IX of [1] and in [2]. The use of stateDiffTrans,
%stateAvgFun, and stateTrans arises when dealing with state space models
%involving, for example, angles. This is analogous to similar issues that
%arise when dealing with angular measurements as described in [3].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems Magazine,
%    vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] I. Arasaratnam and S. Haykin, "Cubature Kalman filters," IEEE
%    Transactions on Automatic Control, vol. 54, no. 6, pp. 1254-1269,
%    Jun. 2009.
%[3] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=length(xPrev);

    if(nargin<5||isempty(xi))
        [xi,w]=fifthOrderCubPoints(xDim);
    end

    if(nargin<7||isempty(stateDiffTrans))
        stateDiffTrans=@(x)x;
    end
    
    if(nargin<8||isempty(stateAvgFun))
        stateAvgFun=@(x,w)calcMixtureMoments(x,w);
    end
    
    if(nargin<9||isempty(stateTrans))
        stateTrans=@(x)x;
    end
    
    xDim=length(xPrev);
    numCubPoints=length(w);
    
    xPropPoints=zeros(xDim,numCubPoints);%Allocate space

    %Calculate the cubature state points
    xPoints=stateTrans(transformCubPoints(xi,xPrev,SPrev));

    %Propagate the cubature state points
    for curP=1:numCubPoints
        xPropPoints(:,curP)=f(xPoints(:,curP));
    end

    %Calculate the predicted state
    xPred=stateAvgFun(xPropPoints,w);

    %Centered, propagated cubature state points
    xDiff=stateDiffTrans(bsxfun(@minus,xPropPoints,xPred));
    xPropCenPoints=bsxfun(@times,xDiff,sqrt(w'));

    %The root prediction covariance.
    SPred=tria([xPropCenPoints, SQ]);
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
