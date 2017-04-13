function [yUpdate, PInvUpdate]=infoFilterUpdate(yPred,PInvPred,z,R,H)
%%INFOFILTERUPDATE Perform the measurement update step in the standard 
%                  linear information filter.
%
%INPUTS: yPred The xDimX1 predicted information state. The information
%              state is the inverse covariance matrix times the target
%              state.
%     PInvPred The xDimXxDim inverse of the predicted state covariance
%              matrix.
%            z The zDim X 1 vector measurement.
%            R The zDim X zDim measurement covariance matrix. This must be
%              positive definite.
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%
%OUTPUTS: yUpdate The xDim X 1 updated (posterior) information state
%                 vector.
%      PInvUpdate The updated xDim X xDim inverse state covariance matrix.
%
%The information filter is algebraically equivalent to the standard linear
%Kalman filter, but allows for track propagation with very uncertain
%states, such as when starting tracks. The implementation of the update
%given here is from the flow chart given in Appendix H of [1].
%
%The standard information filter has the implied linear measurement
%equation
%z=H*x+w
%where z is the measurement, w is the zero-mean Gaussian measurement noise
%with covariance matrix R and x is the true state. However, instead of
%propagating the true state x with its covariance matrix P, the information
%filter propagates
%y=inv(P)*x
%and instead of propagating the true covariance matrix P, it propagates
%PInv=inv(P)
%which means that the filter can be used even when PInv is singular.
%
%More information on information filtering is given in Chapter 7.2 of [2]
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    yUpdate=yPred+H'/R*z;
    PInvUpdate=PInvPred+H'/R*H;
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
