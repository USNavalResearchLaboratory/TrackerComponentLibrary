function [xUpdate,PUpdate,innov,Pzz]=KalmanUpdate(xPred,PPred,z,R,H)
%KALMANUPDATE Perform the measurement update step in the standard linear
%             Kalman filter. This assumes that the predicted target state
%             and the measurement are both multivariate Gaussian
%             distributed.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        PPred The xDim X xDim predicted state covariance matrix.
%            z The zDim X 1 measurement vector.
%            R The zDim X zDim measurement covariance matrix. 
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%
%OUTPUTS: xUpdate The xDim X 1 updated target state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%
%Given a prediction of the target state xPred with an associated covariance
%matrix PPred, and assuming that they are the mean and covariance matrix of
%a Gaussian distribution, under the measurement model
%z=H*x+w
%where z is the measurement, x is the truye target state, H is a zDim*xDim
%matrix and w is additive Gaussian noise, this function updates the mstate
%estimate and its associated covariance estimate.
%
%The Joseph-form covariance update is used for improved numerical
%stability. The algorithm is derived in Chapter 5 of [1].
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zPred=H*xPred;
innov=z-zPred;

Pzz=R+H*PPred*H';
W=PPred*H'/Pzz;

xUpdate=xPred+W*innov;

temp=W*H;
temp=eye(size(temp))-temp;
PUpdate=temp*PPred*temp'+W*R*W';
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
