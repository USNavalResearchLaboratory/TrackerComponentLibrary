function [A,xHyp,PHyp,GateMat]=makeStandardCartOnlyLRMatHyps(xPred,SPred,zCart,SRCart,PD,lambda,rPred,gammaVal,zNative,measJacob)
%%MAKESTANDARDCARTONLYLRMATHYPS Create the likelihood ratio matrix used for
%               2D assignment under standard Gaussian approximations in the
%               coordinate system of the target states using Gaussian-
%               approximated Cartesian position-only measurements and
%               provide all updated state estimates using the Cartesian
%               measurements. This function can also perform simultaneous
%               brute-force gating, returning a matrix indicating whether
%               each measurement gates with each target. This function is
%               appropriate for use in Cartesian-converted measurement
%               tracking without range rate information. This assumes that
%               the false alarms occur due to a Poisson point process with
%               a constant density lambda that is either in the coordinate
%               system of the states (Cartesian) or in the coordinate
%               system of a converted measurement when appropriate
%               parameters are passed. This implements the unitless ratio
%               of [1], except detection probabilities are just multiplied
%               by the target existence probabilities, if provided, and an
%               option for the false alarm density to be specified in a
%               different coordinate system is given.
%
%INPUTS: xPred A numDimXnumTar set of predicted target states. The
%              first posDim Components are position.
%        SPred A numDimXnumDimXnumTar set of lower-triangular square-root
%              state covariance matrices.
%        zCart The zDimXnumMeas set of measurements in Cartesian
%              coordinates.
%       SRCart A zDimXzDimXnumMeas set of lower-triangular square-root
%              measurement covariance matrices. If all of the covariance
%              matrices are the same, then SRCart can be a single zDimXzDim
%              lower triangular matrix.
%           PD The target detection probabilities. This can either be a
%              numTarX1 vector if the targets have different probabilities,
%              or this can be a scalar value if all of the target detection
%              probabilities are the same. This value does not include the
%              effects of gating (set by gammaVal) eliminating measurements
%              from consideration.
%       lambda The false alarm density. When given in Cartesian
%              coordinates, this has units of # false alarm/volume, for
%              example, false alarms/m^3. When given in the coordinate
%              system of the measurement (anticipating a conversion at a
%              single point using a Jacobian), this has units of inverse of
%              the "volume" in the measurement coordinate system. For
%              example, for polar measurements, this might be inverse
%              meters*radians. If this is not in Cartesian coordinates,
%              then the inputs measJacob and zNative must be provided.
%        rPred When using a single-scan tracking algorithm with integrated
%              target existence probabilities, this is the numTarX1 set of
%              target existence probabilities. If omitted or an empty
%              matrix is passed, then rPred=1 is used (everything exists).
%     gammaVal This is the threshold for gating. Gating is performed
%              by a brute-force evaluation of Mahalanobis distances between
%              predicted target locations and the measurements. gammaVal is
%              related to the probability region about the target as
%              described in Chatper 2.3.2 of [2]. For a zDim-dimensional
%              measurement, one can obtain the value of gammaVal for a
%              99.97% probability region using
%              ChiSquareD.invCDF(0.9997,zDim). The default if this
%              parameter is omitted or an empty matrix is passed is Inf,
%              meaning that everything gates.
%  zNative, measJacob, If these two inputs are omitted, then it is assumed
%              that lambda is given in Cartesian coordinates. Otherwise, it
%              is assumed that lambda is given in the native coordinate
%              system of the (Cartesian-converted) measurements and thus
%              has to be "converted" to Cartesian cordinates too. However,
%              such a conversion will typically yield a different lambda at
%              every single point. Thus, the value of lambda given at the
%              measurement location is used. For the ith measurement, this
%              is det(measJacob(zNative(:,i)))*lambda. Thus, zNative is a
%              zOrigDimXnumMeas matrix of measurements and measJacob
%              is a function handle to get the Jacobian of the measurement.
%              For example, if the original measurements were in spherical
%              coordinates, one could use 
%              measJacob=@(z)calcSpherConvJacob(z,0).
%              
%OUTPUTS: A A numTarX(numMeas+numTar) matrix of target-measurement and 
%           missed detection likelihood ratios computed as in [1] for
%           Gaussian states/(converted) measurements. Columns > numMeas
%           hold missed-detection likelihoods. Thus, off-diagonal terms
%           for columns > numMeas are set to 0 and the diagonal terms
%           set to the costs of a missed detection for each given target. 
%           When given rPred, PD is just multiplied by these values before
%           computing the assignment matrix. If a measurement does not
%           gate, then the corresponding entry in A is explicitly set to
%           0.
%      xHyp A numDimXnumTarX(numMeas+1) set of conditionally updated target
%           states (using the prior, the measurement, and the
%           sqrtKalmanUpdate function), one for each association hypothesis
%           in A. The last hypothesis (missed detection) is the same as the
%           prediction value.
%      PHyp A numDimXnumDimXnumTarX(numMeas+1) set of conditionally updated
%           target covariance matrices corresponding to the state estimate
%           in xHyp. Hypotheses are not provided in square root form so as
%           to simplify their usage in JPDAF-style filters.
%   GateMat A numTarXnumMeas boolean matrix indicating whether each target
%           gates with each measurement. 
%
%We are using the dimensionless score function from [1]. This formulation
%in terms of likelihood ratios eliminates the need to compute the volumes
%of detection gates. Of course, if these likelihood ratios are computed in
%Cartesian coordinates, but lambda is in a measurement coordinate system,
%this is just an approximation. The updated target hypotheses are obtained
%using the sqrtKalmanUpdate function.
%
%The assignment matrix returned by this function is commonly used in
%functions such as assign2D, singleScanUpdate, and
%singleScanUpdateWithExistence. This function is useful when performing
%tracking using Cartesian-converted measurements.
%
%REFERENCES:
%[1] Y. Bar-Shalom, S. S. Blackman, and R. J. Fitzgerald, "Dimensionless
%    score function for multiple hypothesis tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 43, no. 1, pp. 392-400, Jan.
%    2007.
%[2] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);

if(nargin<10)
    zNative=[]; 
end

if(nargin<9)
    measJacob=[]; 
end

if(nargin<8||isempty(gammaVal))
    gammaVal=Inf;
end

if(nargin<7||isempty(rPred))
    numTar=size(xPred,2);
    %Assume all targets exist with probability 1.
    rPred=ones(numTar,1);
end

posDim=size(zCart,1);

%The measurement matrix.
H=[eye(posDim,posDim),zeros(posDim,xDim-posDim)];

measUpdateFuns=@(x,S,z,SR)sqrtKalmanUpdate(x,S,z,SR,H);

[A,xHyp,PHyp,GateMat]=makeStandardLRMatHyps(xPred,SPred,measUpdateFuns,zCart,SRCart,PD,lambda,rPred,gammaVal,[],zNative,measJacob);

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
