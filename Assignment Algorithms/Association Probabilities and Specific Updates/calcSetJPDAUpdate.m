function [xCum,PCum,measOrigProb]=calcSetJPDAUpdate(xEst,PEst,A,posDim)
%%CALCSETJPDAUPDATE Compute the measurement update step in the approximate
%              set joint probabilistic data association filter (Set 
%              JPDAF), which is described in [1]. The implementation here
%              is taken from [2], which used sequential 2D assignment to
%              solve a problem that appears to otherwise be NP-hard.
%              Unlike other variants of the JPDAF, one cannot extract
%              simple target-measurement association probabilities.
%
%INPUTS: xEst An xDimXnumTarX(numMeas+1) set of posterior state estimates
%             for each of numtar targets conditioned on which measurement
%             originated from the target. The entry
%             xEst(:,curTar,numMeas+1) is the missed detection hypothesis.
%        PEst An xDimXxDimXnumTarX(numMeas+1) set of covariance matrices
%             associated with the target state estimates in xEst.
%          A  A matrix of positive likelihoods or likelihood ratios (NOT
%             log-likelihood ratios). A is a numTar X (numMeas+numTar)
%             matrix of all-positive likelihoods or likelihood ratios for
%             assigning the target specified by the row to the measurement
%             specified by the column. Columns > numMeas hold
%             missed-detection likelihoods. Thus, off-diagonal terms for
%             columns > numMeas should be set to 0 and the diagonal terms
%             set to the costs of a missed detection for each given target.
%             NOTE: The missed detection hypotheses cannot have zero
%             likelihood.
%      posDim The number of dimensions of xEst that consist of position.
%             The optimization in the Set JPDAF only really makes sense
%             when performed only using the position components of the
%             state, which are assumed to be the first posDim components of
%             the state. If this parameter is omitted or an empty matrix is
%             passed, then posDim=xDim is used.
%
%OUTPUTS: xCum The posterior state estimate of the targets.
%         PCum The posterior covariance matrices associated with the states
%              in xCum.
% measOrigProb A numMeasX2 matrix such that measOrigProb(i,1) is the
%              probability that the ith measurement came from a target and
%              measOrigProb(i,2) is the probability that the ith
%              measurement is a false alarm.
%
%The Set JPDAF is introduced in [1] and a method of efficiently
%approximating it using 2D assignment algorithms is given in [2] and is
%implemented here. Unlike other variants of the JPDAF, the target states in
%the output are considered "unordered", making the act of drawing lines
%between tracks over time difficult. Here, the all targets missed detection
%hypothesis has the original target ordering in the computations. Thus, one
%can often assume that the ordering of the targets in xCum corresponds to
%the ordering in xEst.
%
%The algorithm must generate all target-measurement association hypotheses.
%The total number of hypotheses is given by
%totalNumHyps=num2DTarMeasHyps(numTar,numMeas);
%When that number is large, this algorithm is slow.
%
%REFERENCES:
%[1] L. Svensson, D. Svensson, M. Guerriero, and P. Willett, "Set JPDA
%    filter for multitarget tracking," IEEE Transactions on Signal
%    Processing, vol. 59, no. 10, pp. 4677-4691, Oct. 2011.
%[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target recognition XXII, vol. 874504, Baltimore, MD, 29 Apr. 2013.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xEst,1);
if(nargin<4||isempty(posDim))
   posDim=xDim; 
end

numTar=size(A,1);

numMeas=size(A,2)-numTar;
%Divide each row by the likelihoods of the missed detections and then get
%rid of those columns.
ANew=A(:,1:numMeas);
for curTar=1:numTar
   ANew(curTar,:)=ANew(curTar,:)/A(curTar,numMeas+curTar);
end
A=ANew;

%The cumulative sum of joint association event likelihoods. These are
%needed for normalization.
thetaJAESum=0;
%The first hypothesis is the missed detection hypothesis. If everything in
%A is properly normalized, then the likelihood of a missed detection for
%each target is 1, thus the product (the joint association event) of
%nothing being observed is also 1.
thetaJAECur=1;%Joint association event likelihood.

thetaJAESum=thetaJAESum+thetaJAECur;
curJAE=2;

%measOrigProbs(i,1) will be the probabilities that measurement i comes from
%a target. measOrigProbs(i,2) is the probability that it is a false alarm.
%Here, we initialize it with the likelihood ratio contribution of the first
%joint association event.
measOrigProb=[zeros(numMeas,1),thetaJAECur*ones(numMeas,1)];

%In this approximate implementation of the MMOSPA optimization for the Set
%JPDAF, we SEQUENTIALLY perform the optimization as in [2]. 
xCum=zeros(xDim,numTar);
PCum=zeros(xDim,xDim,numTar);

%Add in the all-missed detection hypothesis
for curTar=1:numTar
    xCum(:,curTar)=xEst(:,curTar,numMeas+1)*thetaJAECur;
    PCum(:,:,curTar)=(PEst(:,:,curTar,1)+xEst(:,curTar,numMeas+1)*xEst(:,curTar,numMeas+1)')*thetaJAECur;
end

%Obtain the likelihood values for all targets.
for numTargetsObserved=1:min(numTar,numMeas)
    %Now, go through all combinations of which targets are observed. The
    %first combination is just the first numTargetsObserved targets.
    obsTar=(0:(numTargetsObserved-1))';
    
    %missedTarSel selects which targets are involved in missed detection
    %hypotheses.
    missedTarSel=true(numTar,1);
    missedTarSel(obsTar+1)=false;
    missedTar=find(missedTarSel);
    while(~isempty(obsTar))
        %For each set of observed targets, determine which measurements are
        %target-originated.
        tarMeas=(0:(numTargetsObserved-1))';
        
        %falseAlarmSel selects which measurements are false alarms.
        falseAlarmSel=true(numMeas,1);
        falseAlarmSel(tarMeas+1)=false;
        while(~isempty(tarMeas))
            %Next, permute the assignment of targets to measurements. The
            %cost is the product of the relevant entries in the association
            %matrix. tarMeas(permIdx(curTar))+1 is assigned to
            %obsTar(curTar)+1
            numPerm=factorial(numTargetsObserved);
            for curPerm=0:(numPerm-1)
                permIdx=unrankPermutation(curPerm,numTargetsObserved);

                %Get indices for the selected elements in the assignment
                %matrix.
                idx=tarMeas(permIdx)*numTar+obsTar+1;
                %The unnormalized likelihood ratio. For this particular
                %hypothesis.
                thetaJAECur=prod(A(idx));
                thetaJAESum=thetaJAESum+thetaJAECur;

                JAEIdx=zeros(numTar,1);
                %Targets not observed get the index of the missed detection
                %hypothesis in this joint association event.
                JAEIdx(missedTar)=numMeas+1;
                
                %Record which measurements are assigned to which targets.
                JAEIdx(obsTar+1)=tarMeas(permIdx)+1;
                
                %Fill the cost matrix using the current partial solution.
                %Rows are from xCum.
                %Columns are from xHyp for the current set of assigned
                %targets.
                c=zeros(numTar,numTar);
                for curRow=1:numTar
                   for curCol=1:numTar
                       curTarEstInJAE=JAEIdx(curCol);

                       c(curRow,curCol)=sum(xCum(1:posDim,curRow).*xEst(1:posDim,curCol,curTarEstInJAE));
                   end
                end
                
                %Use the shortest path algorithm to find the optimal 2D assignment.
                %This is the assignment of each target in xCum to the targets selected
                %in xEst(:,:,JAEIdx(:,curJAE));
                orderList=assign2D(c,true);
                
                for curTar=1:numTar
                    curCol=orderList(curTar);
                    curTarEstInJAE=JAEIdx(curCol);
                    xCum(:,curTar)=xCum(:,curTar)+xEst(:,curTar,curTarEstInJAE)*thetaJAECur;
                    PCum(:,:,curTar)=PCum(:,:,curTar)+(PEst(:,:,curTar,curTarEstInJAE)+xEst(:,curTar,curTarEstInJAE)*xEst(:,curTar,curTarEstInJAE)')*thetaJAECur;
                end
                
                %Add in the likelihood to the false alarm probability
                %values.
                measOrigProb(falseAlarmSel,2)=measOrigProb(falseAlarmSel,2)+thetaJAECur;
                measOrigProb(tarMeas+1,1)=measOrigProb(tarMeas+1,1)+thetaJAECur;
                curJAE=curJAE+1;
            end

            tarMeas=getNextCombo(tarMeas,numMeas);
        end
        obsTar=getNextCombo(obsTar,numTar);
        missedTarSel=true(numTar,1);
        missedTarSel(obsTar+1)=false;
        missedTar=find(missedTarSel);
    end
end

%Normalize the estimates with respect to the association probabilties.
xCum=xCum/thetaJAESum;
PCum=PCum/thetaJAESum;

%Normalize the false alarm probabilities
measOrigProb=measOrigProb./sum(measOrigProb,2);

%Update the covariances (the final term in the sum) for each target and
%deal with finite precision errors.
for curTar=1:numTar
    PCum(:,:,curTar)=PCum(:,:,curTar)-xCum(:,curTar)*xCum(:,curTar)'; 
    
    %Now we shall deal with the possibility that precision errors gave us a
    %negative eigenvalue somewhere.
    [V,D] = eig(PCum(:,:,curTar));
    D(D<0)=0;
    PCum(:,:,curTar)=V*D/V;    
    PCum(:,:,curTar)=(PCum(:,:,curTar)+PCum(:,:,curTar)')/2;%To deal with symmetry errors.
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
