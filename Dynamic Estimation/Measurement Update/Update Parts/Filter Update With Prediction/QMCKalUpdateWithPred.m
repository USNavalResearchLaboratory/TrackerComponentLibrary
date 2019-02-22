function [xUpdate,PUpdate,innov,Pzz,W]=QMCKalUpdateWithPred(z,R,zPred,PzPred,otherInfo)
%%QMCKALUPDATEWITHPRED Given the output of the measurement prediction step
%           from QMCKalMeasPred and a measurement, complete the measurement
%           update step of the quasi-Monte Carlo Kalman filter with
%           additive measurement noise. Separating the measurement
%           prediction step from the rest of the update step can make the
%           creation of multiple measurement association hypotheses from a
%           single target prediction more efficient. The full measurement
%           update function is QMCKalUpdate.
%
%INPUTS: z The zDim X 1 vector measurement.
%        R The zDim X zDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimX1 measurement prediction from the filter.
%   PzPred The zDimXzDim covariance matrix associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the QMCKalMeasPred function.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update.
%
%See the comments to the function QMCKalMeasPred for an example of usage of
%this function. See the comments to QMCKalUpdate for more information on
%the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

innovTrans=otherInfo.innovTrans;
stateTrans=otherInfo.stateTrans;
xPred=otherInfo.xPred;
PPred=otherInfo.PPred;
Pxz=otherInfo.Pxz;

Pzz=PzPred+R;

%The innovation, transformed as necessary to keep values in a desired
%range.
innov=innovTrans(z-zPred);

%The filter gain
W=Pxz/Pzz;

%Updated state estimate
xUpdate=stateTrans(xPred+W*innov);

%Updated state covariance matrix
PUpdate=PPred-W*Pzz*W';

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
