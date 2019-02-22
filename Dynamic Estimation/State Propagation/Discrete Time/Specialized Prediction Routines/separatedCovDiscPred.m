function [xPred,LPred,TPred,PPred]=separatedCovDiscPred(xPrev,PPrev,LPrev,F,Ba,c)
%%SEPARATEDCOVDISCPRED Perform the state prediction step in the separated
%                  covariance filter. This filter takes a maximum assumed
%                  acceleration (or other moment) for the target and
%                  provides the optimum estimate in terms of a cost
%                  function that trades off between estimation accuracy
%                  and estimator delay.
%
%INPUTS: xPrev The xDimX1 state estimate at the previous time-step.
%        PPrev The xDim X xDim state covariance matrix at the previous
%              time-step.
%        LPrev The xDimXzDim delay vector (defined before Equation 12 in
%              [1]) from the previous time step, where zDim is the
%              dimensionality of the position components in LPrev (and also
%              the dimensionality of the measurement). The use of multiple
%              columns represents a choice in how the algorithm was
%              generalized to multiple dimensions.
%           F  The xDim X xDim state transition matrix.
%           Ba An xDimXzDim vector that represents the effects of the
%              maximum possible acceleration (or other moment) on the
%              target state over the duration of the prediction interval.
%              This only matters for the computation of the lag term in the
%              filter and the sign of all of the elements is generally
%              positive. For example, for a 3D linear dynamic model with
%              additive Gaussian noise consisting of a state having 3
%              position and 3 velocity components, one might use
%              Ba=[T^2/2,   0,      0;
%                  0,       T^2/2,  0;
%                  0,       0,      T^2/2;
%                  T,       0,      0;
%                  0,       T,      0;
%                  0,       0,      T]*aMax;
%              where aMax is the maximum allowable acceleration. In 2D,
%              this would be
%              Ba=[T^2/2,   0;
%                  0,       T^2/2;
%                  T,       0;
%                  0,       T]*aMax;
%              and in 1D, as in the paper, this is
%              Ba=[T^2/2;
%                  T]*aMax;
%              However, analogous filters can work with a target state
%              consisting of position velocity and acceleration where Ba
%              adds in the effects of a maximum jerk term.
%            c The confidence region under consideration by the filter.
%              0<c<1. If this parameter is omitted or an empty matrix is
%              passed, the default value of c=0.99 is used.
%
%OUTPUTS: xPred The xDimX1 predicted target state.
%         LPred The xDimXzDim predicted delay vector.
%         TPred The xDimXxDim total error matrix. This is a combination of
%               errors due to measurement noise and filter lag.
%         PPred The xDimXxDim predicted state covariance matrix. This is
%               not needed for the measurement update step, but can be
%               useful for gating or other uses.
%
%The equations for the algorithm are given in Table 1 of [1]. The algorithm
%is presented in 1D for a state consisting of position and velocity.
%However, the equations are given in vector form and can thus be used in
%multiple dimensions, which is done here.
%
%The filter in 1 is 1D. To generalize the filter to 3D, L is redefined so
%that each column contains the lag for one particular dimension of motion.
%To keep the solution the same as in Table 1, L becomes a matrix where the
%elements in each column that do not correspond to components for that
%dimensions of motion are zero.
%
%In [1], no clear method of initializing this type of tracking filter is
%provided. A simple way to initialize the filter would be to use two
%Cartesian converted measurements to obtain a state estimate and covariance
%PInit as one would do with a normal Kalman filter (one could, for example,
%use the KalmanFIRSmoother function) and then LPrev=zeros(xDim,zDim); and
%PPrev=PInit.
%
%REFERENCES:
%[1] G. J. Portmann, J. R. Moore, and W. G. Bath, "Separated covariance
%    filtering," in Proceedings of the IEEE International Radar Conference,
%    Arlington, VA, 7-10 May 1990, pp. 456-460.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(c))
   c=0.99;%99% confidence interval.
end

xPred=F*xPrev;
LPred=F*LPrev+Ba;
TPred=c^2*F*PPrev*F'+LPred*LPred';

%Ensure symmetry
TPred=(TPred+TPred')/2;

if(nargout>3)
    PPred=(TPred-LPred*LPred')/c^2;
    
    %Ensure symmetry
    PPred=(PPred+PPred')/2;
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
