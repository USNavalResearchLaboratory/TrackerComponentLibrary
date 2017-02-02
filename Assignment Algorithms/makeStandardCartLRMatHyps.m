function [A,xHyp,PHyp]=makeStandardCartLRMatHyps(xPred,SPred,H,zCart,SRCart,PD,lambda,rPred,zNative,measJacob)
%%MAKECARTLRMATHYPS Create the likelihood ratio matrix used for 2D
%               assignment under standard Gaussian approximations in the
%               coordinate system of the target states using Gaussian-
%               approximated position-only measurements and provide all
%               Kalman-filter updated state estimates using the Cartesian
%               measurements. This function is appropriate for use in
%               Cartesian-converted measurement tracking. This assumes that
%               the false alarms occur due to a Poisson point process with
%               a constant density lambda that is either in the coordinate
%               system of the states (Cartesian) or in the coordiante
%               system of a converted measurement when appropriate
%               parameters are passed. This implements the unitless ratio
%               of [1], except detection probabilities are just multiplied
%               by the target existence probabilities, if provided and an
%               option for the false alarm density to be specified in a
%               different coordinate system is given.
%
%INPUTS: xPred A numDimXnumTar set of predicted target states.
%        SPred A numDimXnumDimXnumTar set of lower-triangular square-root
%              state covariance matrices.
%            H A zDimXnumDim matrix such that H*xPred extracts the
%              position components of the state (equivalent to those in the
%              Cartesian [converted] measrement).
%        zCart The zDimXnumMeas set of measurements in Cartesian
%              coordinates. This can either be the actual measurement
%              directly in Cartesian coordinates, or this can be an
%              unbiased Cartesian-converted measurement.
%       SRCart A numDimXnumDimXnumTar set of lower-triangular square-root
%              measurement covariance matrices. It is assumed that the
%              innovation covariance matrix for the ith target and jth
%              measurement is 
%              H*SPred(:,:,i)*SPred(:,:,i)'*H'+SRCart*SRCart'. That value
%              is as used in the Kalman filter.
%           PD The target detection probabilities. This can either be a
%              numTarX1 vector if the targets have different probabilities,
%              or this can be a scalar value if all of the target detection
%              probabilities are the same.
%       lambda The false alarm density. When given in Cartesian
%              coordinates, this has units of # false alarm/volume, for
%              example, false alarms/m^3. When given in the coordinate
%              system of the measurement (anticipating a conversion at a
%              single point using a Jacobian), this as units of inverse of
%              the "volume" in the measurement coordinate system. For
%              example, for polar measurements, this might be inverse
%              meters*radians. If this is not in Cartesian coordinates,
%              then the inputs measJacob and zNative must be provided.
%        rPred When using a single-scan tracking algorithm with integrated
%              target existence probabilities, this is the numTarX1 set of
%              target existence probabilities. If omitted or an empty
%              matrix is passed, then rPred=1 is used (everything exists).
%  measJacob,zNative If these two inputs are omitted, then it is assumed
%              that lambda is given in Cartesian cooridnates. Otherwise, it
%              is assumed that lambda is given in the native coordinate system of
%              the (Cartesian-converted) measurements and thus has to be
%              "converted" to Cartesian cordinates too. However, such a
%              conversion will typically yield a different lambda at every
%              single point. Thus, the value of lambda given at the
%              measurement location is used. For the ith measurement, this
%              is det(measJacob(zNative(:,i)))*lambda. Thus, zNative is a
%              zOrigDimXnumMeas matrix of measurements and measJacob
%              is a function handle to get the Jacobian of the measurement.
%              For example, if the original measurements were in spherical
%              coordinates, one could use 
%              measJacob=@(z)calcSpherJacob(z,0).
%              
%OUTPUTS: A A numTarX(numMeas+numTar) matrix of target-measurement and 
%           missed detection likelihood ratios computed as in [1] for
%           Gaussian states/(converted) measurements. Columns > numMeas
%           hold missed-detection likelihoods. Thus, off-diagonal terms
%           for columns > numMeas should be set to 0 and the diagonal terms
%           set to the costs of a missed detection for each given target. 
%           When given rPred, PD is just multiplied by these values before
%           computing the assignment matrix.
%      xHyp A numDimXnumTarX(numMeas+1) set of conditionally updated target
%           states (using the prior, the measurement, and the
%           sqrtKalmanUpdate function), one for each association hypothesis
%           in A. The last hypothesis (missed detection) is the same as the
%           prediction value.
%      PHyp A numDimXnumTarX(numMeas+1) set of conditionally updated target
%           covariance matrices corresponding to the state estimate in
%           xHyp. Hypotheses are not provided in square root form so as to
%           simplify their usage in JPDAF-style filters.
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
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);
zDim=size(zCart,1);
numTar=size(xPred,2);
numMeas=size(zCart,2);
numHyp=numMeas+1;%The extra one is the missed detection hypothesis.

if(nargin<7||isempty(rPred))
    %Assume all targets exist with probability 1.
    rPred=ones(numTar,1);
end

if(isscalar(PD))
    PD=PD*ones(numTar,1);
end

PD=PD.*rPred;

%The likelihood ratio matrix.
A=zeros(numTar,numMeas+numTar);
%Conditional hypothesis updates.
xHyp=zeros(xDim,numTar,numHyp);
PHyp=zeros(xDim,xDim,numTar,numHyp);

for curTar=1:numTar
    xPredCur=xPred(:,curTar);
    SPredCur=SPred(:,:,curTar);
    
    for curMeas=1:numMeas
        zCartCur=zCart(:,curMeas);
        SRCartCur=SRCart(:,:,curMeas);
        
        [xUpdate, SUpdate,innov,Szz]=sqrtKalmanUpdate(xPredCur,SPredCur,zCartCur,SRCartCur,H);
        xHyp(:,curTar,curMeas)=xUpdate;
        PHyp(:,:,curTar,curMeas)=SUpdate*SUpdate';
        
        if(nargin>9&&~isempty(zNative))
            %If the clutter density is given in some local measurement
            %coordinate system, then compute the Jacobian to do a single-
            %point transformation.
            JDet=det(measJacob(zNative(:,curMeas)));
        else
            JDet=1;
        end

        A(curTar,curMeas)=PD(curTar)/(JDet*lambda)*GaussianD.PDFS(innov,zeros(zDim,1),Szz);
    end
    
    %The missed detection likelihood.
    A(curTar,numMeas+curTar)=(1-PD(curTar));
    xHyp(:,curTar,end)=xPredCur;
    PHyp(:,:,curTar,end)=SPredCur*SPredCur';
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
