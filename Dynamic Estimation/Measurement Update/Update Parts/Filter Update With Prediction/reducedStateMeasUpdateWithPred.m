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
%INPUTS: z The zDim X 1 vector measurement.
%        R The zDim X zDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimX1 measurement prediction from the filter.
%   PzPred The zDimXzDim covariance matrix associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the KalmanMeasPred function.
%    MPred The xDimXxDim matrix contributing to the total predicted state
%          covariance matrix based solely on measurement errors.
%    DPred The xDimXzDim matrix of bias coefficients that are supposed to
%          relate target state errors to dynamic model parameter
%          uncertainty.
%
%OUTPUTS: xUpdate The xDimX1 updated state estimate.
%         MUpdate The xDimXxDim updated state covariance matrix
%                 contribution due to measurement errors.
%         DUpdate The xDimXzDim matrix of updated bias coefficients
%                 contributing to the state covariance matrix.
%      innov, Pzz The zDimX1 innovation and a zDimXzDim matrix Pzz that is
%                 akin to an innovation covariance matrix are returned in
%                 case one wishes to analyze the consistency of the
%                 estimator or use those values in gating or likelihood
%                 evaluation.
%               W The gain used in the update. This can be useful when
%                 gating and using the function calcMissedGateCov.
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

Pzz=PzPred+R;

W=PPred*H'/Pzz;
L=eye(xDim,xDim)-W*H;

%Equation 23 in [1].
MUpdate=L*MPred*L'+W*R*W';

%Equation 24 in [1].
DUpdate=L*DPred;

innov=z-zPred;
%Equation 26 in [1].
xUpdate=xPred+W*innov;

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
