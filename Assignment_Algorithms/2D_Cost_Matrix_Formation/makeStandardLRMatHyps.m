function [A,xHyp,PHyp,GateMat]=makeStandardLRMatHyps(xPred,SPred,measUpdateFuns,zMeas,SRMeas,PD,lambda,rPred,gammaVal,measTypes,measJacobDet)
%%MAKESTANDARDLRMATHYPS Create the likelihood ratio matrix used for 2D
%               assignment under standard Gaussian approximations that are
%               used for the prior and posterior of the state. This
%               function allows the measurements to be of different types,
%               but the measurement update functions must all return an
%               innovation and a lower-triangular square root innovation
%               covariance matrix to be used for gating. This function
%               creates all possible target-measurement association
%               hypotheses and performs simultaneous brute-force gating.
%               It is assumed that the false alarms occur due to a Poisson
%               point process with a constant density that is either in the
%               coordinate system of the states or in the coordinate system
%               of a converted measurement when appropriate parameters are
%               passed. This implements the unitless ratio of [1], except
%               detection probabilities are just multiplied by the target
%               existence probabilities, if provided and an option for the
%               false alarm density to be specified in a different
%               coordinate system is given. Additionally, the missed
%               detection cost is modified with the gate probability as
%               described below.
%
%INPUTS: xPred A numDimXnumTar set of predicted target states.
%        SPred A numDimXnumDimXnumTar set of lower-triangular square-root
%              state covariance matrices.
% measUpdateFuns A numTypesX1 or 1XnumTypes cell array containing
%              function handles that are called to perform the measurement
%              update step with each type of measurement. Each function
%              handle has the form
%              curMeasUpdateFun(xPredCur,SPredCur,zMeasCur,SRMeasCur)
%              That is, it takes a state estimate, its lower-triangular
%              square root covariance matrix, a measurement and the
%              lower-triangular square root covariance matrix of the
%              measurement. If only one type of measurement exists, then
%              this can be a function handle. Different measurement types
%              might be, for example, position with and without Doppler.
%        zMeas The zDimXnumMeas set of measurements.
%       SRMeas The zDimXzDimXnumMeas set of lower-triangular square root
%              covariance matrices associated with zMeas. If in a
%              particular measurement type only the first zpDim components
%              of zMeas are used, then one will typically only fill/ use
%              the zpDimXzpDim submatric of a particualar entry in SRMeas.
%              If all of the covariance matrices are the same, then SRMeas
%              can be a single zDimXzDim lower triangular matrix.
%           PD The target detection probabilities. This can either be a
%              numTarX1 vector if the targets have different probabilities,
%              or this can be a scalar value if all of the target detection
%              probabilities are the same. This value does not include the
%              effects of gating (set by gammaVal) eliminating measurements
%              from consideration.
%       lambda The false alarm density. If the false alarm density is the
%              same for all measurement types, then this is a scalar value.
%              If it differs by measurement type (e.g. some measurements
%              lack a range rate component, others have one), then this is
%              a numTypesX1 or 1XnumTypes vector containing the different
%              values. When given in Cartesian coordinates, this has units
%              of # false alarm/volume, for example, false alarms/m^3. When
%              given in the coordinate system of the measurement
%              (anticipating a conversion at a single point using a
%              Jacobian), this has units of inverse of the "volume" in the
%              measurement coordinate system. For example, for polar
%              measurements, this might be inverse meters*radians. If this
%              is not in Cartesian coordinates, then the measJacobDet must
%              be provided.
%        rPred When using a single-scan tracking algorithm with integrated
%              target existence probabilities, this is the numTarX1 set of
%              target existence probabilities. If omitted or an empty
%              matrix is passed, then rPred=1 is used (everything exists).
%     gammaVal This is the threshold for gating. If the gating value is the
%              same for all measurement types, then this is a scalar value.
%              If it differs, then this is a numTypesX1 or 1XnumTypes
%              vector. Gating is performed by a brute-force evaluation of
%              Mahalanobis distances between predicted target locations and
%              the measurements. gammaVal is related to the probability
%              region about the target as described in Chatper 2.3.2 of
%              [2]. For a zDim-dimensional measurement, one can obtain the
%              value of gammaVal for a 99.97% probability region using
%              ChiSquareD.invCDF(0.9997,zDim). The default if this
%              parameter is omitted or an empty matrix is passed is Inf,
%              meaning that everything gates.
%    measTypes A numMeasX1 or 1XnumMeas vector that specifies which
%              measurement update function in measUpdateFuns and which
%              value in lambda should be used for each measurement. The
%              values here range from 1 to numMeas. If there is only one
%              type of measurement, then an empty matrix can be passed.
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
%           sqrtKalmanUpdate function), one for each association hypothesis
%           in A. The last hypothesis (missed detection) is the same as the
%           prediction value.
%      PHyp A numDimXnumDImXnumTarX(numMeas+1) set of conditionally updated
%           target covariance matrices corresponding to the state estimate
%           in xHyp. Hypotheses are not provided in square root form so as
%           to simplify their usage in JPDAF-style filters.
%   GateMat A numTarXnumMeas boolean matrix indicating whether each target
%           gates with each measurement. 
%
%This function builds the dimensionless score function from [1]. Note,
%however that the missed detection probability with gating is 1-PD*PG and
%not just 1-PD When multiple thresholds [gammaVals] are given, the largest
%value of PG is seelected to use). This comes from normalizing the PDF of
%the measurement to the gate. See Section 3.5.3 of [2], where Equation
%3.5.3-6 has this normalizing term. This normalization term means that one
%can't just replace PD with PD*PG. This formulation in terms of likelihood
%ratios eliminates the need to compute the volumes of detection gates. Of
%course, if these likelihood ratios are computed in Cartesian coordinates,
%but lambda is in a measurement coordinate system, this is just an
%approximation. Square root updates are assumed because they tend to be the
%most numerically stable.
%
%Though the gate probability is used in the likelihood computation, it is
%not used in the computation of the covariance matrix for the missed
%detection hypothesis. The fact that a measurement is not seen in the gate
%can potentially mean that the target was detected and its measurement is
%outside of the gate. The issue is discussed in [3] (See Equation 28).
%After using this function, one might wish to use calcMissedGateCov to
%adjust the gate size of the missed detection hypothesis. This doesn't
%really matter if the detection probability is not near 1 but the gate
%probability is near 1.
%
%REFERENCES:
%[1] Y. Bar-Shalom, S. S. Blackman, and R. J. Fitzgerald, "Dimensionless
%    score function for multiple hypothesis tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 43, no. 1, pp. 392-400, Jan.
%    2007.
%[2] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%[3] X. R. Li, "Tracking in clutter with strongest neighbor measurements -
%    Part I: Theoretical analysis," IEEE Transactions on Automatic Control,
%    vol. 43, no. 11, pp. 1560-1578, Nov. 1998.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);
numTar=size(xPred,2);
zDim=size(zMeas,1);
numMeas=size(zMeas,2);
numHyp=numMeas+1;%The extra one is the missed detection hypothesis.

if(isa(measUpdateFuns,'function_handle'))
    temp=measUpdateFuns;
    measUpdateFuns=cell(1,1);
    measUpdateFuns{1}=temp;
end

numMeasUpdateFuns=length(measUpdateFuns);

if(nargin<10||isempty(measTypes))
   measTypes=ones(numMeas,1);
end

if(isscalar(lambda))
    lambda=lambda*ones(numMeasUpdateFuns,1);
end

if(nargin<9||isempty(gammaVal))
    gammaVal=Inf;
end

if(isscalar(gammaVal))
    gammaVal=gammaVal*ones(numMeasUpdateFuns,1);
end

if(nargin<8||isempty(rPred))
    %Assume all targets exist with probability 1.
    rPred=ones(numTar,1);
end

if(isscalar(PD))
    PD=PD*ones(numTar,1);
end

%The gate probability (from an inverse chi-squared PDF). We use the maximum
%gate probability when multiple ones are present.
PG=ChiSquareD.CDF(max(gammaVal),zDim);

%Whether or not the target is detected is affected by whether it exists
%(the rPred).
PD=PD(:).*rPred(:);

if(nargin<11)
    measJacobDet=[];
else
    %They must all be positive.
    measJacobDet=abs(measJacobDet);
end

%If only one is given, assume they are all the same.
if(size(SRMeas,3)==1)
    SRMeas=repmat(SRMeas,[1,1,numMeas]);
end

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

    for curMeas=1:numMeas
        curMeasType=measTypes(curMeas);
        
        curMeasUpdateFun=measUpdateFuns{curMeasType};
        curLambda=lambda(curMeasType);
        
        [xUpdate,SUpdate,innov,Szz]=curMeasUpdateFun(xPredCur,SPredCur,zMeas(:,curMeas),SRMeas(:,:,curMeas));
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

        %To just evaluate the likelihood ratio, we could do
        % A(curTar,curMeas)=PD(curTar)/(JDet*curLambda)*GaussianD.PDFS(innov,zeros(zDim,1),Szz);
        %However, we want to also perform gating of each target to each
        %measurement and not evaluate everything if the measurement does
        %not gate. Thus, we shall do this in two parts. First, we evaluate
        %the Mahalanobis distance necessary for gating, then we evaluate
        %the likelihood ratio above if the target gates.

        %Note that (Szz*Szz')^(-1)=(Szz')^(-1)*Szz^(-1)
        diff=Szz\innov;
        mahabDist=sum(diff.*diff,1);
        
        if(mahabDist<=gammaVal(curMeasType))%If the measurement gates
            GateMat(curTar,curMeas)=true;
            zDim=size(Szz,1);
            A(curTar,curMeas)=PD(curTar)/(JDet*curLambda)*(1/((2*pi)^(zDim/2)*abs(det(Szz))))*exp(-0.5*mahabDist);
        else%It does not gate; assign a zero likelihood ratio.
            A(curTar,curMeas)=0;
        end
    end
    
    %The missed detection likelihood (includes the gating probability).
    A(curTar,numMeas+curTar)=(1-PD(curTar)*PG);
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
