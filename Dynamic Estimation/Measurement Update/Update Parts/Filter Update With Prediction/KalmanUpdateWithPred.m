function [xUpdate,PUpdate,innov,Pzz,W]=KalmanUpdateWithPred(z,R,zPred,PzPred,otherInfo)
%%KALMANUPDATEWITHPRED Given the output of the measurement prediction step
%           from KalmanMeasPred and a measurement, complete the measurement
%           update step of the Kalman filter. Separating the measurement
%           prediction step from the rest of the update step can make the
%           creation of multiple measurement association hypotheses from a
%           single target prediction more efficient. The full measurement
%           update function is KalmanUpdate.
%
%INPUTS: z The zDim X 1 vector measurement.
%        R The zDim X zDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimX1 measurement prediction from the filter.
%   PzPred The zDimXzDim covariance matrix associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the KalmanMeasPred function.
%
%OUTPUTS: xUpdate The xDim X 1 updated target state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update.
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

Pzz=PzPred+R;

innov=z-zPred;

W=Pxz/Pzz;

xUpdate=xPred+W*innov;

temp=W*H;
temp=eye(size(temp))-temp;
PUpdate=temp*PPred*temp'+W*R*W';
%Ensure symmetry
PUpdate=(PUpdate+PUpdate')/2;
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
