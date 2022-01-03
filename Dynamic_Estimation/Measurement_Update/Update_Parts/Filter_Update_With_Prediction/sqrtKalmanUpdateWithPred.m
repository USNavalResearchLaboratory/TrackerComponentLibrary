function [xUpdate, SUpdate,innov,Szz,W]=sqrtKalmanUpdateWithPred(z,SR,zPred,otherInfo)
%%SQRTKALMANUPDATEWITHPRED Given the output of the measurement prediction
%           step from sqrtKalmanMeasPred and a measurement, complete the
%           measurement update step of the square-root Kalman filter.
%           Separating the measurement prediction step from the rest of the
%           update step can make the creation of multiple measurement
%           association hypotheses from a single target prediction more
%           efficient. The full measurement update function is
%           sqrtKalmanUpdate.
%
%INPUTS: z The zDimX1 measurement vector.
%       SR The zDimXzDim lower-triangular square root of the measurement
%          covariance matrix in the native coordinate system of the
%          measurement.
%    zPred The zDimXnumComp measurement predictions from the filter.
%   PzPred The zDimXzDimXnumComp covariance matrices associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the sqrtKalmanMeasPred function.
%
%OUTPUTS: xUpdate The xDimXnumComp updated state vectors.
%         SUpdate The updated xDimXxDimXnumComp lower-triangular square-
%                 root state covariance matrices.
%      innov, Szz The zDimXnumComp innovations and the zDimXzDimXnumComp
%                 square-root innovation covariance matrix are returned in
%                 case one wishes to analyze the consistency of the
%                 estimator or use those values in gating or likelihood
%                 evaluation.
%               W The xDimXzDimXnumComp gain used in the update. This can
%                 be useful when gating and using the function
%                 calcMissedGateCov.
%
%See the comments to the function sqrtKalmanMeasPred for an example of
%usage of this function. See the comments to sqrtKalmanUpdate for more
%information on the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xPred=otherInfo.xPred;
SPred=otherInfo.SPred;
Pxz=otherInfo.Pxz;
H=otherInfo.H;

xDim=size(xPred,1);
numComp=size(xPred,2);
zDim=size(z,1);

xUpdate=zeros(xDim,numComp);
SUpdate=zeros(xDim,xDim,numComp);
innov=zeros(zDim,numComp);
Szz=zeros(zDim,zDim,numComp);
W=zeros(xDim,zDim,numComp);

for k=1:numComp
    Szz(:,:,k)=tria([H*SPred(:,:,k),SR]);
    W=(Pxz(:,:,k)/Szz(:,:,k)')/Szz(:,:,k);
    innov(:,k)=z-zPred(:,k);
    xUpdate(:,k)=xPred(:,k)+W(:,:,k)*innov(:,k);
    temp=W(:,:,k)*H;
    SUpdate(:,:,k)=tria([(eye(size(temp))-temp)*SPred,W(:,:,k)*SR]);
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
