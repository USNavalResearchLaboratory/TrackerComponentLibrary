function [xUpdate,PUpdate,rUpdate,probNonTargetMeas]=singleScanUpdateWithExistence(xHyp,PHyp,PD,rPred,A,algSel1,algSel2,param3)
%%SINGLESCANUPDATEWITHEXISTENCE Perform the measurement update step in a
%                  single-scan tracking algorithm that uses Gaussian  
%                  approximations to represent the target state before and
%                  after the measurement update given a matrix of
%                  likelihoods and a set of track update hypotheses as well
%                  as probabilities that the targets actually exist. This
%                  is essentially the update step in a type of Joint
%                  Integrated Probabilistic Data Association Filter
%                  (JIPDAF), except different options for the algorithms
%                  allow it to also apply to a JIPDAF* as well as other
%                  variant filters. In all instances, it is assumed that
%                  each measurement can be assigned to at most one target,
%                  and each target can be assigned to at most one
%                  measurement.
%
%INPUTS: xHyp An xDimXnumTarXnumHyp set of track states updated with each
%              of the (numHyp-1) measurements for each of the numTar
%              targets, with the last hypothesis being the update for the
%              missed detection hypothesis. The order of the measurement
%              hypotheses is the same as the ordering of the likelihoods in
%              the columns of the A matrix.
%         PHyp An xDimXxDimXnumTarXnumHyp set of covariance matrices for
%              each of the track states for each of the targets updated
%              conditioned on each of the numHyp measurements with the last
%              one being for the missed detection hypothesis.
%           PD The a priori detection probabilities of the targets given
%              that they exist. This can either be a numTarX1 vector if the
%              targets have different detection probabilities. Otherwise,
%              this can be a scalar value if all targets have the same
%              detection probability.
%        rPred A numTarX1 vector of the predicted probabilities of target
%              existence.
%           A  A matrix of positive likelihoods or likelihood ratios (NOT
%              log-likelihood ratios). A is a numTar X (numMeas+numTar)
%              matrix of all-positive likelihoods or likelihood ratios for
%              assigning the target specified by the row to the measurement
%              specified by the column. Columns > numMeas hold
%              missed-detection likelihoods. Thus, off-diagonal terms for
%              columns > numMeas should be set to 0 and the diagonal terms
%              set to the costs of a missed detection for each given
%              target. These are computed as one would for the input to the
%              function singleScanUpdate,except to deal with the likelihood
%              of a target's non-existence, all occurences of the detection
%              probability must be multiplied by the probability that the
%              target exists (the appropriate entry in rHyp). For example,
%              a typical missed detection likelihood ratio for the ith
%              target might be 1-Pd*rHyp(i). where PD is the detection
%              probability. The function makeStandardCartOnlyLRMatHyps can
%              be used to make such likelihood ratio matrices in certain
%              Gaussian/Cartesian problems.
%      algSel1 An optional parameter that selects the assignment algorithm
%              to use. Approximate variants of the algorithms also use
%              algSel2. If this parameter and the next parameter are both
%              omitted, then the algorithm is chosen as described below.
%              Possible values are
%              0) GNN-JIPDA
%              1) JIPDA
%              2) JIPDA*
%              3) Approximate GNN-JIPDA
%              4) Approximate JIPDA
%              5) Naive nearest neighbor JIPDA.
%              6) Approximate naïve nearest neighbor JIPDA.
%      algSel2 An optional parameter that further specifies the algorithm
%              used when algSel1=3,4. If omitted but algSel1 is specified,
%              a default value of 0 is used. If both algSel1 and algSel2
%              are omitted, the algorithm is chosen as described below.
%              algSel2 chooses the approximation for the beta terms in the
%              function calc2DAssignmentProbsApprox (the approxType input)
%              when using algSel1=3,4; the comments to that function
%              describe the values it can take.
%       param3 For the case where algSel1=3,4 and algSel2=0, this is the
%              optional delta input to the function
%              calc2DAssignmentProbsApprox. If omitted or an empty matrix
%              is passed, the default value for that function is used.
%
%OUTPUTS: xUpdate An xDimXnumTar matrix of updated target states.
%         PUpdate An xDimXxDimXnumtar matrix of updated covariance matrices
%                 for the targets.
%         rUpdate A numTarX1 set of updated probabilities of target
%                 existence.
% probNonTargetMeas A numMeasX1 set of posterior probabilities that the
%                 measurements do not come from any of the known targets.
%                 These probabilities play a role in track initiation.
%
%This function basically implements the measurement update step for the
%JIPDAF (and related filters), updating the known targets, their existence
%probabilities, and providing probabilities of the measurements being from
%new tracks. The JIPDAF is described in [1].
%
%The JIPDAF is essentially a JPDAF which also updates probabilities of
%target existence. rPred is a vector of prior target existence
%probabiltiies. The target existence probability update given measurements
%is given in Equation 9.14 of Chapter 9 of [2]. If betas(:,1:(end-1)) holds
%posterior target-measurement association probabilities for known targets
%and betas(:,end) is the set of probabilities that the targets either were
%not detected or do not exist, then the update for the probability of
%target existence is derived as
% P(exists|Z)=P(exists and observed|Z)+P(exists and not observed |Z)
% Well, P(exists and observed|Z) is just sum(betas(:,1:(end-1)),2)
% On the other hand
% P(exists and not observed |Z)=P(exists|Z,not observed)*P(not observed|Z)
% Well, P(not observed|Z)=betas(:,end)
% That leaves 
% P(exists|Z,not observed)=P(exists|not observed)=P(exists and not observed)/P(not observed)
% Well,
% P(exists and not observed)=P(not observed|exists)*P(exists)=(1-PD)*rPred
% and similarly,
% P(not observed)=1-P(observed|exists)*P(exists)=1-PD*rPred
%
%Next, the other difference is in the computation of the target-measurement
%association probabilities. Here, it is simpler to look at [3], where the 
%Track-Oriented Marginal Multiple Target Multi-Bernoulli-Poisson (TOMB/P)
%filter is presented. As discussed in Section IV-A of [3], the TOMB/P
%filter differs from the JIPDA only in that it includes an extra term to
%account for "undetected targets" that one is modeling based on there being
%an expected number of targets present. Removing that term from Equation 68
%of [3], the result is identical to the association probability computation
%in the JPDAF, except one multiplied PD (the a priori detection
%probability) by rPred, the a priori target existence probability. Thus,
%one can use all of the methods available for computing and approximating
%target-measurement association probabilities in the JPDAF. Thus, the
%implementation of the JIPDAF*, which is derived in [4], requires nothing
%special compared to the standard JPDAF*.
%
%After using this function to perform a measurement update, one will
%typically use one-point differencing on all of the measurements to create
%potential new track hypotheses. The initial probabilities of these
%hypotheses are usually taken to be the posterior probability that the
%measurement is not from a known target times another value. Thus, this
%function returns those probabilities for later use. In Equation
%9.14 in Ch. 9.4.1 of [2], the value by which the posterior non-target
%measurement probabilities are multiplied is just taken to be an ad-hoc
%design parameter. This is probably the best way to go. Equation 9.13 of
%Ch. 9.4.1 of [2] and Eq. 69 of [3] take this other value to be something
%based on an assumes spatial density of new objects being discovered.
%
%The prediction step in filtering using such target existence probabilities
%is exactly the same as in any other single-Gaussian filter except the
%existence probabilities rUpdate are typically multiplied by a transition
%probability (to go from existence to non-existence) before the next
%update. This value is a design parameter.
%
%Note that the type of integrated tracking filter described here
%corresponds to one using Markov Chain 1 in [1]. The other option contains
%both probabilities that the target exists as well as probabilities that
%the target is just not visible at all (but still exists).
%
%REFERENCES:
%[1] D. Musicki and R. Evans, "Joint integrated probabilistic data
%    association: JIPDA," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 3, pp. 1093-1099, Jul. 2004.
%[2] S. Challa, M. R. Morelande, D. Musicki, and R. J. Evans, Fundamentals
%    of Object Tracking. Cambridge: Cambridge University press, 2011.
%[3] J. L. Williams, "Marginal multi-Bernoulli filters: RFS derivation of
%    MHT, JIPDA, and association-based MeMBer," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, Jul.
%    2015.
%[4] H. A. P. Blom, E. A. Bloem, and D. Musicki, "JIPDA*: Automatic target
%    tracking avoiding track coalescence," IEEE Transactions on Aerospace
%    and Electronic Systems, vol. 51, no. 2, pp. 962-974, Apr. 2015.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8)
    param3=[];
end

xDim=size(xHyp,1);
numTar=size(A,1);
numMeas=size(A,2)-numTar;
numHyp=numMeas+1;

if(isscalar(PD))
    PD=PD*ones(numTar,1); 
end

if(nargin<6||isempty(algSel1))
    if(numTar==1)
        algSel1=1;
    else
        algSel1=3;
        algSel2=4;
    end
end

if(nargin==5||isempty(algSel2))%If only algSel2 was omitted.
    algSel2=0;
end

%Compute the target-measurement association probabilities under different
%assumptions/ approximations. Note that the missed detection probabilities
%should also include the possibility that the target was not detected
%because it does not exist.
switch(algSel1)
    case 0%GNN-JIPDAF
        betas=calc2DAssignmentProbs(A,true);
    case 1%JPDAF
        betas=calc2DAssignmentProbs(A,true);
    case 2%JIPDAF*
        betas=calcStarBetasBF(A);
    case 3%Approximate GNN-JIPDAF
        betas=calc2DAssignmentProbsApprox(A,algSel2,true,param3);
    case 4%Approximate JIPDAF
        betas=calc2DAssignmentProbsApprox(A,algSel2,true,param3);
    case 5%Naïve nearest neighbor JIPDA
        betas=calc2DAssignmentProbs(A,true);
    case 6%Approximate naïve nearest neighbor JIPDA
        betas=calc2DAssignmentProbsApprox(A,algSel2,true,param3);
    otherwise
        error('Unknown algorithm selected')
end

rPred=rPred(:);

%Updated probabilities of target existence.
rUpdate=((1-betas(:,end))+((1-PD).*rPred)./(1-PD.*rPred).*betas(:,end)).';

%The probability that each measurement did not originate from a known
%target.
probNonTargetMeas=1-sum(betas(:,1:(end-1)),1)';

%Allocate space for the return variables.
xUpdate=zeros(xDim,numTar);
PUpdate=zeros(xDim,xDim,numTar);
if(algSel1==0||algSel1==3)%If a GNN estimate is used in place of the mean
    %Perform 2D assignment
    ALog=log(A);
    tar2Meas=assign2D(ALog,true);
    logLikes=zeros(numTar,1);
    for curTar=1:numTar
        logLikes(curTar)=ALog(curTar,tar2Meas(curTar));
    end

    %Adjust for the index of the missed detection hypothesis.
    sel=tar2Meas>numMeas;
    tar2Meas(sel)=numMeas+1;

    %Set the target estimate to the ML estimate
    for curTar=1:numTar
        maxIdx=tar2Meas(curTar);
        xUpdate(:,curTar)=xHyp(:,curTar,maxIdx);
    end
    
    %Compute the covariance matrix in the manner of the JPDAF, but
    %about the ML estimate, so the result is a MSE matrix estimate.
    for curTar=1:numTar
        x=reshape(xHyp(:,curTar,:),[xDim,numHyp,1]);
        P=reshape(PHyp(:,:,curTar,:),[xDim,xDim,numHyp]);
        [~,PUpdate(:,:,curTar)]=calcMixtureMoments(x,betas(curTar,:),P,xUpdate(:,curTar));
    end
elseif(algSel1==5||algSel1==6)%If a naïve nearest neighbor estimate is used
                   %in place of the mean.
    %The assignment for each target
    logLikes=zeros(numTar,1);
    for curTar=1:numTar
        [maxVal,maxIdx]=max(A(curTar,:));

        %If the missed detection hypothesis is the most likely.
        if(maxIdx>numMeas)
            maxIdx=numMeas+1;
        end

        xUpdate(:,curTar)=xHyp(:,curTar,maxIdx);
        logLikes(curTar)=log(maxVal);
    end

    %Compute the covariance matrix in the manner of the JPDAF, but
    %about the ML estimate, so the result is a MSE matrix estimate.
    for curTar=1:numTar
        x=reshape(xHyp(:,curTar,:),[xDim,numHyp,1]);
        P=reshape(PHyp(:,:,curTar,:),[xDim,xDim,numHyp]);
        [~,PUpdate(:,:,curTar)]=calcMixtureMoments(x,betas(curTar,:),P,xUpdate(:,curTar));
    end
else%JIPDAF/JIPDAF*/approximate JIPDAF
    %This is just the mixture mean and covariance matrix (conditioned on
    %the target existing).
    for curTar=1:numTar
        x=reshape(xHyp(:,curTar,:),[xDim,numHyp,1]);
        P=reshape(PHyp(:,:,curTar,:),[xDim,xDim,numHyp]);
        [xUpdate(:,curTar),PUpdate(:,:,curTar)]=calcMixtureMoments(x,betas(curTar,:),P);
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
