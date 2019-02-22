function [xUpdate,MUpdate,DUpdate,innov,Pzz,W]=reducedStateUpdate(xPred,PPred,MPred,DPred,z,R,H)
%%REDSTATEFILTERUPDATE Perform the measurement update step in the reduced
%               state estimator. This is a filter that takes measurements
%               in Cartesian coordinates and assumes that the dynamic model
%               includes an unknown parameter that is somehow bounded, but
%               whose contribution is modeled with a mean and covariance
%               matrix. The filter separates contributions due to
%               measurement errors and dynamic model mismatch errors.
%
%INPUTS: xPred The xDimX1 predicted state estimate.
%        PPred The xDimXxDim predicted total state covariance estimate.
%              This is provided by the state prediction step.
%        MPred The xDimXxDim matrix contributing to the total predicted
%              state covariance matrix based solely on measurement errors.
%        DPred The xDimXzDim matrix of bias coefficients that are supposed
%              to relate target state errors to dynamic model parameter
%              uncertainty.
%            z The zDimX1 measurement vector.
%            R The zDimXzDim measurement covariance matrix.
%            H The optional zDimXxDim measurement matrix. The measurement
%              is modeled as z=H*x+noise. If this parameter is omitted or
%              an empty matrix is passed, then H will be taken as a
%              zDimXxDim identity matrix followed by columns of zeros
%              (Assuming that zDim<=xDim. Otherwise, H must be provided).
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
%The filter is taken from [1]. Other applications are described in [1].
%Note that the formulation of the dynamic model requires that the
%dimensionality of the measurements does not vary.
%
%In [1] and [2], no clear method of initializing this type of tracking
%filter is provided. A simple way to initialize the filter would be to use
%two Cartesian converted measurements to obtain a state estimate and
%covariance as one would do with a normal Kalman filter (one could, for
%example, use the KalmanFIRSmoother function) and then set MPrev to the
%covariance value obtained while setting DPrev to zero.
%
%REFERENCES:
%[1] P. Mookerjee and F. Reifler, "Reduced state estimator for systems with
%    parametric inputs," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 2, pp. 446-461, Apr. 2004.
%[2] P. Mookerjee and F. Reifler, "Reduced state estimators for consistent
%    tracking of maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 41, no. 2, pp. 608-619, Apr. 2005.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);

if(nargin<7||isempty(H))
	zDim=size(z,1); 
    H=[eye(zDim,zDim),zeros(zDim,xDim-zDim)];
end

%Equation 20 in [1].
Pzz=H*PPred*H'+R;

%Ensure symmetry
Pzz=(Pzz+Pzz')/2;

W=PPred*H'/Pzz;
L=eye(xDim,xDim)-W*H;

%Equation 23 in [1].
MUpdate=L*MPred*L'+W*R*W';

%Ensure symmetry
MUpdate=(MUpdate+MUpdate')/2;

%Equation 24 in [1].
DUpdate=L*DPred;

innov=z-H*xPred;
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
