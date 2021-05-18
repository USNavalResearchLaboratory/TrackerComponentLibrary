function boolVal=linTargetIsUntrackable(targetParams,PD,lambda,AbsTol,RelTol,maxIter)
%%LINTARGETISUNTRACKABLE Assuming that an asymptotic inoovation covariance
%           is available or that a target follows a particular specified
%           linear dynamic model (e.g. position, velocity, and
%           acceleration, and it could have correlation terms like a Gauss-
%           Markov model) with known process noise covariance and Cartesian
%           measurement covariance, and assuming that a probability of
%           detection is known, determine whether a target is untrackable.
%           A target is deemed untrackable if the cost function used in a
%           multiple hypothesis tracker (MHT)  is always greater than or
%           equal to  the null hypothesis than for any track hypothesis
%           when given asymptotic track accuracy parameters. If a track is
%           untrackable, even an MHT of infinite length would not be able
%           to track the target. This assumes Cartesian measurements. The
%           use of range rate, complex target amplitude, micro-Doppler and
%           other information could potentially render untrackable targets
%           trackable.
%
%INPUTS: targetParams This can be one of two types of inputs. If an
%          asymptotic innovation covariance matrix is known, then this is
%          the innovation covariance matrix. Otherwise, this is a structure
%          holding parameters for a linear system so that an asymptotic
%          innovation covariance matrix can be computed. The members of the
%          structure are:
%           H The zDim X xDim measurement matrix such that H*x+w is the
%             measurement, where x is the state and w is zero-mean Gaussian
%             noise with covariance matrix R.
%           F The xDim X xDim state transition matrix The state at
%             discrete-time k+1 is modeled as F times the state at time k
%             plus zero- mean Gaussian process noise with covariance matrix
%             Q.
%           R The zDim X zDim measurement covariance matrix.
%           Q The xDim X xDim process noise covariance matrix. This can
%             not be a singular matrix.
%       PD The optional detection probability of the target at each scan.
%          PD should not be 1.
%   lambda The false alarm density in Cartesian coordinates. For example,
%          this might be 1e-8 false alarms per unit volume in cubic meters.
%   RelTol The maximum relative error tolerance allowed when computing the
%          asymptotic predictive covariance matrix. This is a positive
%          scalar. If omitted or an empty matrix is passed, the default
%          value of 1e-13 is used. This value is only used if targetParams
%          is a structure. 
%   AbsTol The absolute error tolerance allowed, a positive scalar. When
%          computing the asymptotic predictive covariance matrix. This is a
%          positive scalar. If omitted or an empty matrix is passed, the
%          default value of 1e-10 is used. This value is only used if
%          targetParams is a structure. 
%  maxIter An optional integer specifying the maximum number of iterations
%          to use when cmputing the asymptotic predictive covariance
%          matrix. By default, if omitted, this is 5000. This value is only
%          used if targetParams is a structure. 
%
%As discussed in Chapter 7.5 of [1], track oriented MHTs declare target
%existence based on finite sums of log-likelihood ratios. These ratios are
%typically made into dimensionless score functions as in Chapter 7.5 of [1]
%and in [2]. For the probability that a measurement z is from a particular
%target, the contribution to the score function is f(z)*PD/lambda, where
%f(z) is the likelihood of the measurement conditioned on the predicted
%state estimate. The null hypothesis receives a score of 1-PD. If the null
%hypothesis is always greater than or equal to the track association
%hypothesis then a target is untrackable.
%
%Given a Cartesian measurement model (z=H*x+w where x is the target state,
%w is measurement noise with covariance matrix R and H is the measurement
%matrix) and assuming a linear dynamic model (xpredicted=F*xPrior+v where v
%is noise with covariance matrix Q) where the state transition matrix F and
%the process noise covariance matrix Q are constant (implying a uniform
%revisit rate) as well as having a known detection probability of PD,  one
%determine an asymptotic lower bound on the estimation accuracy of a
%predicted state measurement using the RiccatiPredNoClutter function. Given
%PPred, a lower bound on the state precision covariance, the score function
%contribution for a target-originated measurement is maximum when the
%predicted measurement equals the measurement value. From Chapter 7.5.3 of
%[1], for a linear measurement model, this is
%PD/(lambda*sqrt(det(2*pi*SPred))) where SPred is the innovation covariance
%matrix, which as given in the derivation of the Kalman filter in Chapter 5
%of [3] is H*PPred*H'+R. Thus, this function determines whether this
%maximum score for a track-originated measurement is greater than the
%alternative (1-PD). If not, then the target is untrackable.
%
%The trackability bound here can be considered something of a "stealth"
%bound, because for a fixed false alarm level in a system, a target can
%become untrackable either by lowing its detection probability sufficiently
%(becoming "stealth") and/or by increasing its maneuverability, which,
%here, is reflected in increasign the process noise covariance.
%
%Note that if the innovation covariance matrix is passed for targetParams,
%this function just uses its size and determiniant.
%
%EXAMPLE:
% T=1;
% targetParams=[];
% targetParams.F=FPolyKal(T,6,1);
% q0=1;
% targetParams.Q=QPolyKal(T,6,1,q0);
% targetParams.R=diag([10;10;10]);
% targetParams.H=[1,0,0,0,0,0;
%    0,1,0,0,0,0;
%    0,0,1,0,0,0];
% PD=0.5;
% lambda=1e-8;
% boolVal=linTargetIsUntrackable(targetParams,PD,lambda)
% %The above example is trackable (boolVal=false). However, if the
% %detection probability decreases/ the false alarm rate increases, then
% %one gets an untrackable target. Consider:
% PD=0.1;
% lambda=1e-6;
% boolVal1=linTargetIsUntrackable(targetParams,PD,lambda)
%
%REFERENCES:
%[1] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%[2] Y. Bar-Shalom, S. S. Blackman, and R. J. Fitzgerald, "Dimensionless
%    score function for multiple hypothesis tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 43, no. 1, pp. 392-400, Jan.
%    2007.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(AbsTol))
    AbsTol=1e-13;
end

if(nargin<5||isempty(RelTol))
    RelTol=1e-10;
end

if(nargin<6||isempty(maxIter))
   maxIter=5000; 
end

if(isstruct(targetParams))
    H=targetParams.H;
    F=targetParams.F;
    R=targetParams.R;
    Q=targetParams.Q;
    m=size(H,1);

    %Asymptotic predicted state covariance matrix with given PD assuming
    %completely correct association hypotheses.
    PPred=RiccatiPredNoClutter(H,F,R,Q,PD,AbsTol,RelTol,maxIter);

    %Asymptotic innovation covariance matrix.
    detSPred=det(H*PPred*H'+R);
else
    S=targetParams;
    m=size(S,1);
    detSPred=det(S);
end

maxIsTargetScore=PD/(lambda*sqrt((2*pi)^m*detSPred));
missedDetectScore=1-PD;

boolVal=maxIsTargetScore<=missedDetectScore;
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
