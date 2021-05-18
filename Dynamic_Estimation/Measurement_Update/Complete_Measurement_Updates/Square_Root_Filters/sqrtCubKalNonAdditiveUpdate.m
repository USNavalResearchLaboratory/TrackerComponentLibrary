function [xUpdate, SUpdate,innov,Szz,W]=sqrtCubKalNonAdditiveUpdate(xPred,SPred,z,SR,h,xi,w,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%%SQRTCUBKALNONADDITIVEUPDATE Perform the measurement update step in the   
%                      square root cubature Kalman filter with non-additive
%                      measurement noise.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular predicted state covariance
%              matrix.
%            z The zDim X 1 vector measurement.
%           SR The zDim X zDim lower-triangular square root of the
%              measurement covariance matrix in the native coordinate
%              system of the measurement.
%            h A function handle for the measurement function. This can
%              either be of the form h(xStacked), where xStacked is [x;w],
%              where x is the state and w is the measurement noise value,
%              or this can be of the form h(x,w), where the inputs are
%              separately specified.
%           xi An (xDim+zDim) X numCubPoints matrix of cubature points. If
%              this and the next parameter are omitted or empty matrices
%              are passed, then fifthOrderCubPoints(xDim) is used. It is
%              suggested that xi and w be provided to avoid needless
%              recomputation of the cubature points.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points.
%   innovTrans An optional function handle that computes and optionally
%              transforms the value of the difference between the
%              observation and any predicted points. This is called as
%              innovTrans(a,b) and the default if omitted or an empty
%              matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%              able to handle sets of values. For a zDimX1 measurement,
%              either of the inputs could be zDimXN in size while one of
%              the inputs could be zDimX1 in size.  This only needs to be
%              supplied when a measurement difference must be restricted
%              to a certain range. For example, the innovation between two
%              angles will be 2*pi if one angle is zero and the other
%              2*pi, even though they are the same direction. In such an
%              instance, a function handle to the
%              wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%              appropriate parameters should be passed for innovTrans.
%   measAvgFun An optional function handle that, when given N measurement
%              values with weights, produces the weighted average. This
%              function only has to be provided if the domain of the
%              measurement is not linear. For example, when averaging
%              angular values, then the function meanAng should be used.
% stateDiffTrans An optional function handle that, like innovTrans does for
%              the measurements, takes an xDimXN matrix of N differences
%              between states and transforms them however might be
%              necessary. For example, a state containing angular
%              components will generally need to be transformed so that the
%              difference between the angles is wrapped to -pi/pi.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one might generally
%              want to bind it to the region +/-pi.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         SUpdate The updated xDim X xDim lower-triangular square root
%                 state covariance matrix.
%      innov, Szz The zDimX1 innovation and the zDimXzDim square root
%                 innovation covariance matrix are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The gain used in the update. This can be useful when
%                 gating and using the function calcMissedGateCov.
%
%%The mathematics behind the function sqrtCubKalNonAdditiveUpdate are
%described in more detail in Section IX of [1]
%
%The optional parameters innovTrans and measAvgFun are not described in the
%reference [1], but allow for possible modifications to the filter as
%described in [2] The parameters have been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%measAvgFun=@(z,w)[calcMixtureMoments(z(1,:),w);
%                  meanAng(z(2,:),w')];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    zDim=size(z,1);
    xDim=size(xPred,1);
    
    if(nargin<6||isempty(xi))
        [xi,w]=fifthOrderCubPoints(xDim+zDim);
    end

    if(nargin<8||isempty(innovTrans))
        %The function just returns the input.
        innovTrans=@(a,b)bsxfun(@minus,a,b);
    end
    
    if(nargin<9||isempty(measAvgFun))
       measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
    end
    
    if(nargin<10||isempty(stateDiffTrans))
       stateDiffTrans=@(x)x; 
    end
    
    if(nargin<11||isempty(stateTrans))
       stateTrans=@(x)x; 
    end
    
    numCubPoints=size(xi,2);
    sqrtW=sqrt(w);
    
    %This let's allows us to parameterize the function either with the
    %noise stacked on the state or as two separate inputs.
    if(nargin(h)>1)
        hStacked=@(x)h(x(1:xDim),x((xDim+1):end));
    else
        hStacked=h;
    end
    
    %Predicted cubature state points
    S=[SPred,               zeros(xDim,zDim);
       zeros(zDim,xDim),SR];
    xPredStackedPoints=transformCubPoints(xi,[xPred;zeros(zDim,1)],S);
    xPredStackedPoints(1:xDim,:)=stateTrans(xPredStackedPoints(1:xDim,:));
    
    %Predicted, centered cubature state points
    xPredCenPoints=bsxfun(@times,bsxfun(@minus,xPredStackedPoints(1:xDim,:),xPred),sqrtW');

    %Predicted cubature measurement points
    zPredPoints=zeros(zDim,numCubPoints);
    for curP=1:numCubPoints
        zPredPoints(:,curP)=hStacked(xPredStackedPoints(:,curP));
    end
    
    %Measurement prediction
    zPred=measAvgFun(zPredPoints,w);
    
    %The innovation
    innov=innovTrans(z,zPred);

    %Centered, predicted cubature measurement points
    zPredCenPoints=bsxfun(@times,innovTrans(zPredPoints,zPred),sqrtW');

    %Root innovation covariance
    Szz=tria(zPredCenPoints);

    %The cross covariance.
    Pxz=xPredCenPoints*zPredCenPoints';

    %The filter gain
    W=(Pxz/Szz')/Szz;

    %Updated state estimate
    xUpdate=stateTrans(xPred+W*innov);

    %Updated state root covariance
    SUpdate=tria([stateDiffTrans(xPredCenPoints-W*zPredCenPoints),W*SR]);
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
