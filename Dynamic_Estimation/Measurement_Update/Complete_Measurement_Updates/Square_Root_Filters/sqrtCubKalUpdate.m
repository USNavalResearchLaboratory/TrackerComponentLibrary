function [xUpdate,SUpdate,innov,Szz,W]=sqrtCubKalUpdate(xPred,SPred,z,SR,h,xi,w,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%SQRTCUBKALUPDATE Perform the measurement update step in the square root 
%                 cubature Kalman filter with additive measurement noise.
%                 Unlike the non-square root version, the covariance update
%                 is such that finite precision problems cannot lead to a
%                 non-positive (semi-)definite covariance. However, unlike
%                 the non-square root version, all of the cubature weights
%                 must be positive.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular square root predicted state
%              covariance matrix.
%            z The zDim X 1  measurement vector.
%           SR The zDim X zDim lower-triangular square root of the
%              measurement covariance matrix in the native coordinate
%              system of the measurement.
%            h A function handle for the measurement function that the
%              state as its argument.
%           xi An xDim X numCubPoints matrix of cubature points. If this
%              and the next parameter are omitted or empty matrices are
%              passed, then fifthOrderCubPoints(xDim). It is suggested that
%              xi and w be provided to avoid needless recomputation of the
%              cubature points.
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
% stateDiffTrans An optional function handle that takes an xDimXN matrix of
%              N differences between states and transforms them however
%              might be necessary. If not transformation is necessary, this
%              parameter can be omitted or an empty matrix passed.
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
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%If the function h needs additional parameters beyond the state, then the
%parameters can be passed by using an anonymous function as the function
%handle. For example, suppose that the measurement function is measFunc and
%it needs the additional parameters param1 and param2. In this instance,
%rather than using
%h=@measFunc
%one should use
%h=@(x)measFunc(x,param1,param2)
%This way, every time sqrtCubKalUpdate calls measFunc (via h) with a
%different x, those two parameters are always passed.
%
%The mathematics behind the function sqrtCubKalUpdate are described in more
%detail in Section IX of [1] and in [2]. Note that this is essentially one
%type of square root "unscented Kalman filter" with additive noise. One
%simply has to provide the filter with the appropriate cubature points and
%weights.
%
%The optional parameters innovTrans and measAvgFun are not described in
%references [1] and [2], but allow for possible modifications to the filter
%as described in [3] The parameters have been added to allow the filter to
%be used with angular quantities. For example, if the measurement consisted
%of range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%measAvgFun=@(z,w)[calcMixtureMoments(z(1,:),w);
%                  meanAng(z(2,:),w')];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] D. F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems Magazine,
%    vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] I. Arasaratnam and S. Haykin, "Cubature Kalman filters," IEEE
%    Transactions on Automatic Control, vol. 54, no. 6, pp. 1254-1269,
%    Jun. 2009.
%[3] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(xPred,1);
    
    if(nargin<6||isempty(xi))
        [xi,w]=fifthOrderCubPoints(xDim);
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
    
    zDim=size(z,1);
    numCubPoints=size(xi,2);
    sqrtW=sqrt(w);
    %Calculate the updated estimates

    %Predicted cubature state points
    xPredPoints=stateTrans(transformCubPoints(xi,xPred,SPred));

    %Predicted, centered cubature state points
    xPredCenPoints=bsxfun(@times,stateDiffTrans(bsxfun(@minus,xPredPoints,xPred)),sqrtW');

    %Predicted cubature measurement points
    zPredPoints=zeros(zDim,numCubPoints);
    for curP=1:numCubPoints
        zPredPoints(:,curP)=h(xPredPoints(:,curP));
    end
    
    %Measurement prediction.
    zPred=measAvgFun(zPredPoints,w);
    
    %Centered, predicted cubature measurement points, transformed as
    %necessary to keep the values within a desired range.
    zPredCenPoints=bsxfun(@times,innovTrans(zPredPoints,zPred),sqrtW');

    %Root innovation covariance
    Szz=tria([zPredCenPoints,SR]);

    %The cross covariance.
    Pxz=xPredCenPoints*zPredCenPoints';

    %The filter gain
    W=(Pxz/Szz')/Szz;

    %The innovation, transformed as necessary to keep values in a desired
    %range.
    innov=innovTrans(z,zPred);
    
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
