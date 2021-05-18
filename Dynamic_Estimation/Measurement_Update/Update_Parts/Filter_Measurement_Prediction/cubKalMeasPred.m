function [zPred, PzPred,otherInfo]=cubKalMeasPred(xPred,PPred,zDim,h,xi,w,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%%CUBKALMEASPRED Perform the measurement prediction part of the measurement
%           update step of the cubature Kalman filter with additive
%           measurement noise. The function cubKalUpdateWithPred can be
%           used to complete the measurement update. Separating the
%           measurement prediction step from the rest of the update step
%           can make the creation of multiple measurement association
%           hypotheses from a single target prediction more efficient.
%           The full measurement update function is cubKalUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states.
%        PPred The xDimXxDimXnumComp predicted state covariance matrices. 
%         zDim The dimensionality of the output of h.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%           xi An xDimXnumCubPoints matrix of cubature points. If this
%              and the next parameter are omitted or empty matrices are
%              passed, then fifthOrderCubPoints(xDim) is used. It is 
%              suggested that xi and w be provided to avoid needless
%              recomputation of the cubature points.
%            w A numCubPointsX1 vector of the weights associated with the
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
%OUTPUTS: zPred The zDimXnumComp measurement predictions from the filter.
%        PzPred The zDimXzDimXnumComp covariance matrices associated with
%               zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to cubKalUpdateWithPred
%               when updating with a measurement.
%
%The mathematics behind the cubature Kalman filter are described in more
%detail in Section IX of [1] and in [2]. See the comments to cubKalUpdate
%for more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%cubKalUpdate in one step as with using cubKalMeasPred followed by
%cubKalUpdateWithPred.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% PPred=[28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13];
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% R=eye(zDim,zDim);%Measurement covariance matrix.
% %The update in one step.
% [xUpdate,PUpdate,innov,Pzz,W]=cubKalUpdate(xPred,PPred,z,R,h);
% %The update in two steps.
% [zPred, PzPred,otherInfo]=cubKalMeasPred(xPred,PPred,zDim,h);
% [xUpdate1,PUpdate1,innov1,Pzz1,W1]=cubKalUpdateWithPred(z,R,zPred,PzPred,otherInfo);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;PUpdate1(:)-PUpdate(:);innov1(:)-innov;Pzz1(:)-Pzz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] D. F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems Magazine,
%    vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] I. Arasaratnam and S. Haykin, "Cubature Kalman filters," IEEE
%    Transactions on Automatic Control, vol. 54, no. 6, pp. 1254-1269,
%    Jun. 2009.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(xPred,1);
    numComp=size(xPred,2);

    if(nargin<5||isempty(xi))
        [xi,w]=fifthOrderCubPoints(xDim);
    end

    if(nargin<7||isempty(innovTrans))
        %The function just returns the input.
        innovTrans=@(a,b)bsxfun(@minus,a,b);
    end
    
    if(nargin<8||isempty(measAvgFun))
        measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
    end

    if(nargin<9||isempty(stateDiffTrans))
        stateDiffTrans=@(x)x; 
    end
    
    if(nargin<10||isempty(stateTrans))
        stateTrans=@(x)x; 
    end

    numCubPoints=size(xi,2);
    
    zPred=zeros(zDim,numComp);
    PzPred=zeros(zDim,zDim,numComp);
    xPredCenPoints=zeros(xDim,numCubPoints,numComp);
    zPredCenPoints=zeros(zDim,numCubPoints,numComp);
    Pxz=zeros(xDim,zDim,numComp);
    
    for k=1:numComp
        %cholSemiDef is used instead of chol in case a positive semi-definite
        %covariance matrix is passed.
        SPred=cholSemiDef(PPred,'lower');
        %Predicted cubature state points
        xPredPoints=stateTrans(transformCubPoints(xi,xPred,SPred));

        %Predicted cubature measurement points
        zPredPoints=zeros(zDim,numCubPoints);
        for curP=1:numCubPoints
            zPredPoints(:,curP)=h(xPredPoints(:,curP));
        end

        %Measurement prediction.
        zPred=measAvgFun(zPredPoints,w);

        %Centered, predicted cubature measurement points, transformed as
        %necessary to keep the values within a desired range.
        zPredCenPoints(:,:,k)=innovTrans(zPredPoints,zPred);
        xPredCenPoints(:,:,k)=stateDiffTrans(bsxfun(@minus,xPredPoints,xPred));
        for curP=1:numCubPoints
            diff=zPredCenPoints(:,curP);
            PzPred(:,:,k)=PzPred(:,:,k)+w(curP)*(diff*diff');
            Pxz(:,:,k)=Pxz(:,:,k)+w(curP)*xPredCenPoints(:,curP)*diff';
        end
        %Pxz is not needed for the measurement prediction, but we compute it
        %here, so that it need not be recomputed again and again if
        %cubKalUpdateWithPred is called for multiple measurements.
    end
    
    otherInfo.innovTrans=innovTrans;
    otherInfo.stateDiffTrans=stateDiffTrans;
    otherInfo.stateTrans=stateTrans;
    otherInfo.xPredCenPoints=xPredCenPoints;
    otherInfo.zPredCenPoints=zPredCenPoints;
    otherInfo.xPred=xPred;
    otherInfo.w=w;
    otherInfo.Pxz=Pxz;
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
