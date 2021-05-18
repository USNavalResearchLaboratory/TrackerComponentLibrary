function [xPred, SPred]=sqrtDiscCubKalNonAdditivePred(xPrev,SPrev,f,SQ,xi,w,stateDiffTrans,stateAvgFun,stateTrans)
%%SQRTDISCCUBKALNONADDITIVEPRED Perform the discrete-time prediction step   
%                 that comes with the square root implementation of the
%                 cubature Kalman filter with non-additive process noise.
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        SPrev The xDim X xDim lower-triangular square root of the state
%              covariance matrix at the previous time-step.
%            f A function handle for the state transition function. This
%              can either be of the form f(xStacked), where xStacked is
%              [x;w], where x is the state and w is the process noise
%              value, or this can be of the form f(x,v), where the inputs
%              are separately specified. The function returns the
%              propagated state given the realization of the process noise.
%           SQ The xDimX xDim lower-triangular square root of the process
%              noise covariance matrix.
%           xi An (xDim+noiseDim) X numCubPoints matrix of cubature points
%              for the joint state and process noise.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points. These must be all postive.
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
%
%The mathematics behind the function sqrtDiscCubKalNonAdditivePred are 
%described in more detail in Section IX of [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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
    noiseDim=size(SQ,1);
    numCubPoints=length(w);
    
    %This let's allows us to parameterize the function either with the
    %noise stacked on the state or as two separate inputs.
    if(nargin(f)>1)
        fStacked=@(x)f(x(1:xDim),x((xDim+1):end));
    else
        fStacked=f;
    end

    xPropPoints=zeros(xDim,numCubPoints);%Allocate space

    %Calculate the augmented prior cubature state points
    SStacked=[SPrev,               zeros(xDim,noiseDim);
              zeros(noiseDim,xDim),SQ];
    xStackedPoints=transformCubPoints(xi,[xPrev;zeros(noiseDim,1)],SStacked);
    
    xStackedPoints(1:xDim,:)=stateTrans(xStackedPoints(1:xDim,:));
    
    %Propagate the cubature state points
    for curP=1:numCubPoints
        xPropPoints(:,curP)=fStacked(xStackedPoints(:,curP));
    end

    %Calculate the predicted state
    xPred=stateAvgFun(xPropPoints,w);

    %Centered, propagated cubature state points
    xPropCenPoints=bsxfun(@times,stateDiffTrans(bsxfun(@minus,xPropPoints,xPred)),sqrt(w'));

    %The root prediction covariance.
    SPred=tria(xPropCenPoints);
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
