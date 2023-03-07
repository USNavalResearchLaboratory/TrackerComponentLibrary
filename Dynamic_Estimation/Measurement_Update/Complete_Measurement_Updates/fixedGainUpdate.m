function [xUpdate,innov,PUpdate]=fixedGainUpdate(W,xPred,z,H,PPred,R)
%FIXEDGAINUPDATE Perform the measurement update step in a fixed-gain linear
%                filter. This can be an alpha-beta or similar filter. The
%                asymptotic gain for a particular dynamic model can be
%                found using the function findAsymptoticGain. If a prediction
%                covariance matrix is available, it can be updated using the
%                suboptimal gain and the measurement covariance matrix.
%                Otherwise, PPred and R can be omitted.
%
%INPUTS: W The gain in a constant gain filter.
%    xPred The xDim X 1 predicted target state.
%        z The zDim X 1 vector measurement.
%        H The zDim X xDim measurement matrix for a linear measurement
%          model. That is z=H*x.
% PPred, R If provided, these will be used to obtain PUpdate using the
%          Joseph's form of the Kalman filter covariance update equation,
%          which can handle suboptimal gains. If these are omitted, then
%          PUpdate will not be returned.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%           innov The zDimX1 innovation of the filter. This is sometimes
%                 used in gating and likelihood evaluation.
%         PUpdate This output is optional and requires that PPred and R be
%                 provided. This is the posterior covariance matrix using the
%                 supplied Kalman filter gain.
%
%The update in a fixed gain linear filter is just
%xUpdate=xPred+W*(z-H*xPred);
%This encompassed alpha-beta-gamma-style filters, since if the state
%is just position, W=alpha; if the state is position and velocity, then
%W=[alpha;beta/T], where T is the sampling period between measurements, if
%the state is position velocity and acceleration, then
%W=[alpha;beta/T;gamma/(2*T^2)], although some authors define gamma such
%that W=[alpha;beta/T;gamma/(T^2)].
%
%This function does not provide state covariance matrix estimates by default,
%because if the filter has reached its asymptotic performance level, the
%posterior state covariance is just PPostAsymp, that one can get from the
%function findAsymptoticGain when computing the gain, and the measurement
%prediction covariance is F*PPostAsymp*F'+Q. If one provides R and Ppred,
%however, the function can return an updated covariance matrix.
%
%The prediction step in a fixed gain filter is just F*xUpdate,
%where F is the state transition matrix.
%
%The Kalman filter update is described in Chapter 5 of [1]. This is the same,
%except the Kalman filter gain is fixed.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zPred=H*xPred;
innov=z-zPred;
xUpdate=xPred+W*innov;

if(nargout>2)
    temp=W*H;
    temp=eye(size(temp))-temp;
    PUpdate=temp*PPred*temp'+W*R*W';
    %Ensure symmetry
    PUpdate=(PUpdate+PUpdate')/2;
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
