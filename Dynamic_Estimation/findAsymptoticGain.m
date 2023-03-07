function [W,PPostAsymp]=findAsymptoticGain(H,F,R,Q)
%%FINDASYMPTOTICGAIN  Find the gain for a fixed-gain linear filter that
%                     uses the asymptotic discrete-time Kalman filter gain.
%                     Filters like this are al-pha-beta-gamma-style
%                     filters. However, rather than providing alpha, beta,
%                     and gamma, this function provides the gain directly.
%                     When given an appropriate  1D motion model, one can
%                     then extract alpha, beta, gamma, and similar
%                     higher-order parameters, if desired. This gain can be
%                     used in the function fixedGainUpdate.
%
%INPUTS: H The zDim X xDim  measurement matrix such that H*x+w is the
%          measurement, where x is the state and w is zero-mean Gaussian
%          noise with covariance matrix R.
%        F The xDim X xDim state transition matrix. The state at discrete-
%          time k+1 is modeled as F times the state at time k plus zero-
%          mean Gaussian process noise with covariance matrix Q
%        R The zDim X zDim measurement covariance matrix.
%        Q The xDim X xDim process noise covariance matrix.
%
%OUTPUTS: W The gain in a constant gain filter. That is, if xk is the state
%           estimate at time k, then the state estimate at time k+1 is
%           F*xk+W*(z-H*F*xk) where z is the observation at time k+1.
% PPostAsymp The asymptotic posterior (after a measurement update) 
%           covariance of the state estimate when using a fixed-gain
%           filter. This is just the solution to the Riccatti equation. If
%           the asymptotic prior (after prediction, before a measurement
%           update) covariance matrix is desired, then that is just
%           F*PPostAsymp*F'+Q.
%
%This assumes that the sampling rate is of a fixed duration T so that F and
%Q are constants. If one uses a 1D motion model, then the coefficients from
%an alpha-beta-gamma-style filter can be backed out of W since if the state
%is just position, W=alpha; if the state is position and velocity, then
%W=[alpha;beta/T]; if the state is position velocity and acceleration, then
%W=[alpha;beta/T;gamma/(2*T^2)], although some authors define gamma such
%that W=[alpha;beta/T;gamma/(T^2)]. Either way, this function can be used
%to solve for the values of alpha beta and gamma that are asymptotically
%optimal.
%
%The algorithm is described in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "A General Solution to Optimal Fixed-Gain (Alpha-Beta-
%    Gamma Etc.) Filters," IEEE Signal Processing Letters, vol. 22, no. 7,
%    pp. 901-904,  Jul. 2015.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, find the asymptotic posterior covariance as a solution to the
%Riccatti equation.
PPostAsymp=RiccatiPostNoClutter(H,F,R,Q);
%Next, find the corresponding Kalman filter gain.
W=PPostAsymp*H'/R;

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
