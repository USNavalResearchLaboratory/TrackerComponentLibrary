function beta=calcStarBetasBF(A)
%%CALCSTARBETASBF Calculate target-measurement association probabilities,
%                 as in the JPDAF* given an association matrix using a
%                 brute-force method.                
%
%INPUTS: A A matrix of positive likelihoods or likelihood ratios (NOT log-
%          likelihood ratios). A is a numTar X (numMeas+numTar) matrix of
%          all-positive likelihoods or likelihood ratios for assigning the
%          target specified by the row to the measurement specified by the
%          column. Columns > numMeas hold missed-detection likelihoods.
%          Thus, off-diagonal terms for columns > numMeas should be set to
%          0 and the diagonal terms set to the costs of a missed detection
%          for each given target. NOTE: The missed detection hypotheses
%          cannot have zero likelihood.
%
%OUTPUTS: beta A numTar X (numMeas+1) matrix of probabilities of assigning
%              the target given by the row to the measurement given by the
%              column. The final column is a set of missed detection
%              probabilities.
%
%The JPDAF* is the same as the JPDAF except certain joint association
%events are thrown out. The JPDAF* is originally from [1].
%
%REFERENCES:
%[1] H. A. Blom and E. A. Bloem, "Probabilistic data association avoiding
%    track coalescence," IEEE Transactions on Automatic Control, vol. 45,
%    no. 2, pp. 247-259, Feb. 2000.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTar=size(A,1);
numMeas=size(A,2)-numTar;
%Divide each row by the likelihoods of the missed detections and then get
%rid of those columns.
ANew=A(:,1:numMeas);
for curTar=1:numTar
   ANew(curTar,:)=ANew(curTar,:)/A(curTar,numMeas+curTar);
end
A=ANew;

%Allocate space for the target-measurement association probabilities.
%The value in beta(:,numMeas+1) gives the probability of a missed
%detection. The non-normalized betas shall be computed, then the betas
%shall be appropriately normalzed.
beta=zeros(numTar,numMeas+1);

%Assuming that all of the likelihood ratios are properly normalized, the
%contribution to the betas for everything being missed detections is just
%one.
beta(:,end)=1;

for numTargetsObserved=1:min(numTar,numMeas)
    %Now, go through all combinations of which targets are observed. The
    %first combination is just the first numTargetsObserved targets.
    obsTar=(0:(numTargetsObserved-1))';
    while(~isempty(obsTar))
        %For each set of observed targets, determine which measurements are
        %target-originated.
        tarMeas=(0:(numTargetsObserved-1))';
        while(~isempty(tarMeas))
            %Next, permute the assignment of targets to measurements. The
            %cost is the product of the relevant entries in the association
            %matrix. tarMeas(permIdx(curTar))+1 is assigned to
            %obsTar(curTar)+1
            numPerm=factorial(numTargetsObserved);
            maxLike=0;
            maxLikeIdx=0;
            for curPerm=0:(numPerm-1)
                permIdx=unrankPermutation(curPerm,numTargetsObserved);

                %Get indices for the selected elements in the assignment
                %matrix.
                idx=tarMeas(permIdx)*numTar+obsTar+1;
                %The unnormalized likelihood ratio.
                likeRat=prod(A(idx));
                if(likeRat>=maxLike)
                    maxLike=likeRat;
                    maxLikeIdx=permIdx;
                end                
            end
            
            %Now, add the likelihood ratio contribution of the most likely
            %permutation of measurements to targets to the appropriate beta
            %term for each target.
            obsTarIdx=1;
            for curTar=1:numTar
                %If the target was observed.
                if(obsTarIdx<=numTargetsObserved&&obsTar(obsTarIdx)+1==curTar)
                    %A target observed hypothesis for this target.
                    measSel=tarMeas(maxLikeIdx(obsTarIdx))+1;
                    beta(curTar,measSel)=beta(curTar,measSel)+maxLike;
                    obsTarIdx=obsTarIdx+1;
                else
                    %A missed detection hypothesis for this target.
                    beta(curTar,end)=beta(curTar,end)+maxLike;
                end
            end
            
            tarMeas=getNextCombo(tarMeas,numMeas);
        end
        obsTar=getNextCombo(obsTar,numTar);
    end
end

%Normalize the betas so that sum of the probabilities of each target being
%assigned to a measurement or a missed detection sums to one.
beta=bsxfun(@rdivide,beta,sum(beta,2));

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
