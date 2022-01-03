function [xUpdate, SUpdate,innov,Szz,W]=sqrtKalmanUpdate(xPred,SPred,z,SR,H)
%SQRTKALMANUPDATE Perform the measurement update step in a square root
%                 implementation of the standard linear Kalman filter.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular predicted square-root state
%              covariance matrix.
%            z The zDim X 1 measurement vector.
%           SR The zDim X zDim lower-triangular square root of the
%              measurement covariance matrix in the native coordinate
%              system of the measurement.
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix SR*SR'.
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
%The mathematics behind the specific square root implementation used here
%are described in Appendix G of [1]. This is algebraically equivalent to
%KalmanUpdate, but there are benefits of using a square root covarinace
%representation.
%
%REFERENCES
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zPred=H*xPred;

Pxz=SPred*SPred'*H';

Szz=tria([H*SPred,SR]);
W=(Pxz/Szz')/Szz;

innov=z-zPred;
xUpdate=xPred+W*innov;
temp=W*H;

SUpdate=tria([(eye(size(temp))-temp)*SPred,W*SR]);
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
