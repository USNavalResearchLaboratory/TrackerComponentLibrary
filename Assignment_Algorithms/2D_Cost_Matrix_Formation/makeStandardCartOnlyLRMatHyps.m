function [A,xHyp,PHyp,GateMat]=makeStandardCartOnlyLRMatHyps(xPred,SPred,zCart,SRCart,SRCartMax,PD,lambda,rPred,gammaVal,measJacobDet)
%%MAKESTANDARDCARTONLYLRMATHYPS Create the likelihood ratio matrix used for
%         2D assignment under standard Gaussian approximations in the
%         coordinate system of the target states using Gaussian-
%         approximated Cartesian position-only measurements and provide
%         updated state estimates when tracks gate with the measurements. 
%         Brute-force gating is used. This function is appropriate for use
%         in Cartesian-converted measurement tracking without range rate
%         information. This assumes that the false alarms occur due to a
%         Poisson point process with a constant density lambda that is
%         either in the coordinate system of the states (Cartesian) or in
%         the coordinate system of a converted measurement when appropriate
%         parameters are passed. This implements the unitless ratio of [1],
%         except the missed detection probability with gating is 1-PD*PG
%         and not just 1-PD due to Section 3.5.3 of [2] (for multiple
%         gammaVal values, the larged PG is chosen). Additionally, the
%         covariance inflation due to gating applied to a missed detection
%         hypothesis (due to [3], Equation 28) can be used, if desired.
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
%    SRCartMax A zDimXzDimXnumTar set of maximum measurement covariance
%              matrices that can be taken into consideration for gating and
%              will be used to compute the covariance inflation of [3]
%              (Equation 28) for missed detection hypotheses. If only a
%              single zDimXzDim matrix is passed, it is assumed the same
%              for all targets. If an empty matrix is passed, then the
%              covariance inflation is not performed. The inflation is
%              insignificant if the gate probability PG is almost 1 and
%              PD is less than 1 by a decent amount. On the other hand, if
%              PD is high but PG is low, then this can be the difference
%              between track divergence or not.
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
%              then the inputs measJacobDet must be provided.
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
% measJacobDet If this input is omitted, then it is assumed that lambda is
%              given in the same coordinate system as the target state
%              (typically Cartesian). Otherwise, this is a length numMeas
%              array such that measJacobDet(i) is
%              det(measJacob(zNative(:,i))) where measJacob is a function
%              that computes the Jacobian matrix of the measurement and
%              zNative is the measurement in the original coordinate system
%              of the measurement, not the coordinate system of the state.
%              For example, if the original measurements were in spherical
%              coordinates, measJacob could be calcSpherConvJacob(z,0), and
%              thus measJacobDet would be determinates of the outputs of
%              calcSpherConvJacob evaluated at all of the measurements.
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
%           sqrtKalmanMeasPred and sqrtKalmanUpdateWithPred functions,
%           which are equivalent to the sqrtKalmanUpdate function), one for
%           each association hypothesis in A. For hypotheses that do not
%           gate, then the corresponding entry in xHyp is all zeros. The
%           last hypothesis (missed detection) is the same as the
%           prediction value, but the covariance might be inflated as in
%           [3] if the SRCartMax input is provided.
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
%this is just an approximation.
%
%The assignment matrix returned by this function is commonly used in
%functions such as assign2D, singleScanUpdate, and
%singleScanUpdateWithExistence. This function is useful when performing
%tracking using Cartesian-converted measurements.
%
%Developments of the score function are given in [1] and relevant parts of
%Section 3.5.3 and other sections of [2].
%
%REFERENCES:
%[1] Y. Bar-Shalom, S. S. Blackman, and R. J. Fitzgerald, "Dimensionless
%    score function for multiple hypothesis tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 43, no. 1, pp. 392-400, Jan.
%    2007.
%[2] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%[3] X. R. Li, "Tracking in clutter with strongest neighbor measurements -
%    Part i: Theoretical analysis," IEEE Transactions on Automatic Control,
%    vol. 43, no. 11, pp. 1560-1578, Nov. 1998.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);

if(nargin<10)
    measJacobDet=[];
else
    %They must all be positive.
    measJacobDet=abs(measJacobDet);
end

if(nargin<9||isempty(gammaVal))
    gammaVal=Inf;
end

numTar=size(xPred,2);
zDim=size(zCart,1);
numMeas=size(zCart,2);

if(nargin<8||isempty(rPred))
    %Assume all targets exist with probability 1.
    rPred=ones(numTar,1);
end

%If only one is given, assume they are all the same.
if(size(SRCart,3)==1)
    SRCart=repmat(SRCart,[1,1,numMeas]);
end

%If only one is given, assume they are all the same.
if(~isempty(SRCartMax)&&size(SRCartMax,3)==1)
    SRCartMax=repmat(SRCartMax,[1,1,numTar]);
end

if(~isempty(SRCartMax)&&isempty(zCart))
    %Get the correct measurement dimensionality to use to update the missed
    %detection covariance matrix with PG for the possibility that the
    %measurement is actually outside of the gate.
    zDim=size(SRCartMax,1);
end

if(isscalar(PD))
    PD=PD*ones(numTar,1);
end

%The gate probability (from an inverse chi-squared PDF). Choose PG to be
%the maximum gate when multiple gammaVal values are given.
PG=ChiSquareD.CDF(max(gammaVal(:)),zDim);

%Whether or not the target is detected is affected by whether it exists
%(the rPred).
PD=PD(:).*rPred(:);

%The measurement matrix.
H=[eye(zDim,zDim),zeros(zDim,xDim-zDim)];

numHyp=numMeas+1;%The extra one is the missed detection hypothesis.

%The likelihood ratio matrix.
A=zeros(numTar,numMeas+numTar);
%The gating matrix
GateMat=false(numTar,numMeas);

%Conditional hypothesis updates.
xHyp=zeros(xDim,numTar,numHyp);
PHyp=zeros(xDim,xDim,numTar,numHyp);
for curTar=1:numTar
    xPredCur=xPred(:,curTar);
    SPredCur=SPred(:,:,curTar);

    %The missed detection likelihood (includes the gating probability).
    A(curTar,numMeas+curTar)=(1-PD(curTar)*PG);
    xHyp(:,curTar,numHyp)=xPredCur;
    PHyp(:,:,curTar,numHyp)=SPredCur*SPredCur';
    
    if(numMeas>0||~isempty(SRCartMax))
        %Get the filter measurement prediction once (not for all
        %measurements). This is independent of the measurement covariance
        %matrix.
        [zPred,PzPred,otherInfo]=sqrtKalmanMeasPred(xPredCur,SPredCur,H);
    
        %If upper bounds on the measurement covariances to consider for
        %each target are provided.
        if(~isempty(SRCartMax))
            %The gain corresponding to a gate with the maximum considered
            %Cartesian covariance matrix for that target.
            W=calcSqrtKalmanGain(SRCartMax(:,:,curTar),otherInfo);
            Pzz=PzPred+(SRCartMax(:,:,curTar)*SRCartMax(:,:,curTar)');
            PMissed=calcMissedGateCov(PHyp(:,:,curTar,numHyp),Pzz,W,PD(curTar),gammaVal,PG);
            PHyp(:,:,curTar,numHyp)=PMissed;
        end
    end
    
    for curMeas=1:numMeas
        %First, we evaluate the Mahalanobis distance necessary for gating,
        %then we evaluate the measurement update and the likelihood ratio
        %if the target gates.
        innov=zCart(:,curMeas)-zPred;
        
        Pzz=PzPred+(SRCart(:,:,curMeas)*SRCart(:,:,curMeas)');
        mahabDist=innov'*inv(Pzz)*innov;
        if(mahabDist>gammaVal)
            %If it does not gate, then assign a zero likelihood ratio.
            A(curTar,curMeas)=0;
            continue;
        end
        
        %It gates.
        GateMat(curTar,curMeas)=true;
        
        %Perform the measurement update.
        [xUpdate, SUpdate,innov,Szz]=sqrtKalmanUpdateWithPred(zCart(:,curMeas),SRCart(:,:,curMeas),zPred,otherInfo);
 
        xHyp(:,curTar,curMeas)=xUpdate;
        PHyp(:,:,curTar,curMeas)=SUpdate*SUpdate';
        
        if(~isempty(measJacobDet))
            %If the clutter density is given in some local measurement
            %coordinate system, then use the Jacobian to do a single-
            %point transformation.
            JDet=measJacobDet(curMeas);
        else
            JDet=1;
        end

        %Evaluate the likelihood ratio.
        A(curTar,curMeas)=PD(curTar)/(JDet*lambda)*GaussianD.PDFS(innov,zeros(zDim,1),Szz);
    end
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
