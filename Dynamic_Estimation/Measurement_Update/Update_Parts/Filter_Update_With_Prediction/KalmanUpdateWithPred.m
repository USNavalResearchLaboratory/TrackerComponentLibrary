function [xUpdate,PUpdate,innov,Pzz,W]=KalmanUpdateWithPred(z,R,zPred,PzPred,otherInfo)
%%KALMANUPDATEWITHPRED Given the output of the measurement prediction step
%           from KalmanMeasPred and a measurement, complete the measurement
%           update step of the Kalman filter. Separating the measurement
%           prediction step from the rest of the update step can make the
%           creation of multiple measurement association hypotheses from a
%           single target prediction more efficient. The full measurement
%           update function is KalmanUpdate.
%
%INPUTS: z The zDimX1 measurement vector.
%        R The zDimXzDim measurement covariance matrix associated with z.
%    zPred The zDimXnumComp measurement predictions from the filter.
%   PzPred The zDimXzDimXnumComp covariance matrices associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the KalmanMeasPred function.
%
%OUTPUTS: xUpdate The xDimXnumComp updated target state vectors.
%         PUpdate The updated xDimXxDimXnumComp state covariance matrices
%                 associated with xUpdate.
%      innov, Pzz The zDimXnumComp innovations and the zDimXzDimXnumComp
%                 innovation covariance matrices are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The xDimXzDimXnumComp gains used in the update. This can
%                 be useful when gating and using the function
%                 calcMissedGateCov.
%
%See the comments to the function KalmanMeasPred for an example of usage of
%this function. See the comments to KalmanUpdate for more information on
%the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xPred=otherInfo.xPred;
PPred=otherInfo.PPred;
Pxz=otherInfo.Pxz;
H=otherInfo.H;

zDim=size(z,1);
xDim=size(xPred,1);
numComp=size(xPred,2);

xUpdate=zeros(xDim,numComp);
PUpdate=zeros(xDim,xDim,numComp);
innov=zeros(zDim,numComp);
Pzz=zeros(zDim,zDim,numComp);
W=zeros(xDim,zDim);
for k=1:numComp
    Pzz(:,:,k)=PzPred(:,:,k)+R;
    innov(:,k)=z-zPred(:,:,k);

    W(:,:,k)=Pxz(:,:,k)/Pzz(:,:,k);

    xUpdate(:,k)=xPred(:,k)+W(:,:,k)*innov(:,k);

    temp=W(:,:,k)*H;
    temp=eye(size(temp))-temp;
    PUpdate(:,:,k)=temp*PPred(:,:,k)*temp'+W(:,:,k)*R*W(:,:,k)';
    %Ensure symmetry
    PUpdate(:,:,k)=(PUpdate(:,:,k)+PUpdate(:,:,k)')/2;
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
