function beta=calcStarBetaskBest(A,K)
%%CALCSTARBETASKBEST Calculate target-measurement association probabilities,
%                    as in the JPDAFStar given an association matrix using
%                    only the k-best hypotheses for the computation of the
%                    target-measurement association probabilities rather
%                    than goign through all joint association events.
%
%INPUTS:   A  A numTar X numMeas matrix of all-positive likelihood
%              ratios for assigning the target specified by the row to the
%              measurement specified by the column.
%          K  The number of hypotheses to generate.
%
%OUTPUTS:  beta A numTar X (numMeas+1) matrix of probabilities of assigning
%               the target given by the row to the measurement given by the
%               column. The final column is a set of missed detection
%               probabilities.
%
%The general idea behind this function is the same as that of the function
%calcStarBetasBF. However, rather than generating all possible
%joint association events, only the k-best events are generated.
%
%To reformulate the problem such that a 2D assignment algorithm can be
%used, the logarithm of the likelihood ratios must be taken, then the cost
%values are exponentiated to get traditional probabilities.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%In order to be able to use the k-Best 2D assignment algorithm, the
%association matrix must be augmented with extra columns that represent
%missed detection events. When the likelihood ratios in A have been
%appropriately computed, the likelihood ratio of a missed detection event
%is equal to one.

numTar=size(A,1);
numMeas=size(A,2);

%Adjust the cost matrix so that 2D assignment can be use.
A=log(A);

%Alocate space for the results.
beta=zeros(numTar,numMeas+1);

%Augment the A matrix to handle missed detection hypotheses. These are
%supposed to have a likelihood ratio of 1, so  the logarithm should be
%zero.
AClut=1-eye(numTar);
AClut(AClut==1)=-inf;
[col4rowBest,row4ColBest,gainBest]=kBest2DAssignment([A,AClut],K,true);

%Undo the logarithm from before.
gainBest=exp(gainBest);

%Adjust K to reflect the fact that fewer than K hypotheses might have been
%present.
K=size(col4rowBest,2);

%For each hypothesis, determine the set of measurements that are
%target-originated and save them order in tarMeas. Also find the set of
%targets that are observed and save them in obsTar.
tarMeas=zeros(numTar,K);
obsTar=zeros(numMeas,K);
for curHyp=1:K
    tarMeas(:,curHyp)=col4rowBest(:,curHyp);

    %Remove elements that are missed detections.
    misSel=col4rowBest(:,curHyp)>numMeas;
    tarMeas(misSel,curHyp)=0;

    %Sort the result.
    tarMeas(:,curHyp)=sort(tarMeas(:,curHyp));
    obsTar(:,curHyp)=sort(row4ColBest(1:numMeas,curHyp));
end

%For each set of observed targets and measurements that are
%target-originated, go through and only keep the most likely hypothesis.
%Since everything is already ordered in decreasing gain, the first
%occurrence of each set of observed targets and target-originated
%measurements is the most likely.
validHyps=true(K,1);
for curHyp=1:K
    if(validHyps(curHyp)==false)
        continue;
    end
    validHyps(curHyp)=false;
    
    for HypS=(curHyp+1):K
        if(isequal(tarMeas(:,curHyp),tarMeas(:,HypS))&&isequal(obsTar(:,curHyp),obsTar(:,HypS)))
            validHyps(HypS)=false;
        end
    end
    
    %Add the ML hypothesis with the given set of observed targets and
    %target-originated measurements to the beta terms.
    for curTar=1:numTar
        curMeas=col4rowBest(curTar,curHyp);
        if(curMeas>numMeas)
            %If it is a missed detection.
             beta(curTar,end)=beta(curTar,end)+gainBest(curHyp);
        else
            beta(curTar,curMeas)=beta(curTar,curMeas)+gainBest(curHyp);
        end
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
