function [xUpdate,MUpdate,DUpdate,innov,Pzz,W]=reducedStateMeasUpdateWithPred(z,R,zPred,PzPred,otherInfo,MPred,DPred)
%%REDUCEDSTATEMEASUPDATEWITHPRED Given the output of the measurement
%           prediction step from reducedStateMeasPred and a measurement,
%           complete the measurement update step of the reduced state
%           estimator. Separating the measurement prediction step from the
%           rest of the update step can make the creation of multiple
%           measurement association hypotheses from a single target
%           prediction more efficient. The full measurement update function
%           is reducedStateUpdate.
%
%INPUTS: z The zDimX1 vector measurement.
%        R The zDimXzDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimXnumComp measurement predictions from the filter.
%   PzPred The zDimXzDimXnumComp covariance matrices associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the KalmanMeasPred function.
%    MPred The xDimXxDimXnumComp matrices contributing to the total predicted
%          state covariance matrix based solely on measurement errors.
%    DPred The xDimXzDimXnumComp matrices of bias coefficients that are
%          supposed to relate target state errors to dynamic model
%          parameter uncertainty.
%
%OUTPUTS: xUpdate The xDimXnumComp updated state estimates.
%         MUpdate The xDimXxDimXnumComp updated state covariance matrices
%                 contribution due to measurement errors.
%         DUpdate The xDimXzDimXnumComp matrices of updated bias
%                 coefficients contributing to the state covariance matrix.
%      innov, Pzz The zDimXnumComp innovations and a zDimXzDimXnumComp set
%                 of matrices matrix Pzz that are akin to an innovation
%                 covariance matrices are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDimXnumComp gains used in the updates. This can
%                 be useful when gating and using the function
%                 calcMissedGateCov.
%
%See the comments to the function reducedStateMeasPred for an example of
%usage of this function. See the comments to reducedStateUpdate for more
%information on the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xPred=otherInfo.xPred;
PPred=otherInfo.PPred;
H=otherInfo.H;
xDim=size(xPred,1);
numComp=size(xPred,2);
zDim=size(z,1);

xUpdate=zeros(xDim,numComp);
MUpdate=zeros(xDim,xDim,numComp);
DUpdate=zeros(xDim,zDim,numComp);
innov=zeros(zDim,numComp);
Pzz=zeros(zDim,zDim,numComp);
W=zeros(xDim,zDim,numComp);
for k=1:numComp
    Pzz(:,:,k)=PzPred(:,:,k)+R;
    W=PPred(:,:,k)*H'/Pzz(:,:,k);
    L=eye(xDim,xDim)-W(:,:,k)*H;

    %Equation 23 in [1].
    MUpdate=L*MPred(:,:,k)*L'+W(:,:,k)*R*W(:,:,k)';

    %Equation 24 in [1].
    DUpdate(:,:,k)=L*DPred(:,:,k);

    innov(:,k)=z-zPred(:,k);
    %Equation 26 in [1].
    xUpdate(:,k)=xPred(:,k)+W(:,:,k)*innov(:,k);
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
