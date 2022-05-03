function [xPostSet,PPostSet,muPost,xMerged,PMerged]=multipleModelUpdate(AlgSel,xPredSet,PPredSet,z,R,measUpdateFuns,muPrev,Lambda,numStateDims,numMixDims)
%%MULTIPLEMODELUPDATE Given a set of predicted target state estimates and
%                     covariance matrices under different models,
%                     approximate the first two moments of the posterior
%                     distributions using the selected algorithm for
%                     handling multiple models. For purposes of computing
%                     the necessary likelihoods, it is assumed that the
%                     target state can be approximated as being jointly
%                     Gaussian with the measurement. With the exception of
%                     the GPB1 algorithm, The models can have different
%                     dimensionalities with only the first few states
%                     playing a role in the mixing. The function
%                     multipleModelPred should be used for predicting the
%                     models.
%
%INPUTS: algSel A parameter indicating the algorithm to use for updating
%               the models. The same algorithm must be used with this
%               update as with the prediction using the multipleModelPred
%               function. Possible values are:
%               'IMM'  The interacting multiple model estimator described
%                      in Chapter 11.6.6 of [3] with the modifications in
%                      [1], [2]  and [4] to make it handle multiple models
%                      of different dimensionalities.
%               'AMM'  The autonomous multiple model estimator, which is
%                      the same as the static model estimator of Chapter
%                      11.6.2 of [3].
%               'GPB1' The generalized pseudo-Bayesian estimator 1,
%                      described in Chapter 11.6.4 of [3]. All models must
%                      have the same dimensionality.
%               'GPB2' The generalized pseudo-Bayesian estimator 2,
%                      described in Chapter 11.6.5 of [3], with changes
%                      similar to that of [1], [2], and [4] to make it work
%                      with models having different dimensionalites.
%      xPredSet The set of predicted state estimates for the different
%               models, stored in a matrix. The models can have different
%               dimensionalities but the number of rows of the matrix is 
%               the dimensionality of the largest model. Models with fewer
%               state elements are zero padded. For the IMM, AMM, and GPB1
%               estimators, this is a xDimMaxXnumModels matrix, where
%               xDimMax is the length of the longest state vector in a
%               model. For the GBP2 estimator, this is a
%               xDimMaxXnumModelsXnumModels matrix as it contains all
%               hypotheses of model switches. The function
%               multipleModelPred can be used to get the prediction from a
%               past update step.
%      PPredSet The set of covariance matrices associated with the
%               solutions in xPredSet. For the IMM, AMM, and GBP1
%               estimators, this is a xDimMaxXxDimMaxXnumModels matrix,
%               where extra elements for models with fewer than xMaxDim
%               dimensions are filled with zeros. For the GBP2 estimator,
%               this is xDimMaxXxDimMaxXnumModelsXnumModels for all
%               possible combinations of model switches.
%             z A zDimX1 measurement that will be used to update all of the
%               models. If a missed detection event occurred, this should
%               be an empty matrix.
%             R A zDimXzDim measurement covariance matrix that will be used
%               when updating all of the models. An empty matrix can be
%               passed if one wishes to perform the hypothesis mixing but
%               not an update (a missed detection event).
% measUpdateFuns A function handle used for updating all of the models or a
%               numModelsX1 cell array of numModels function handles for
%               updating each model when differences in the types of states
%               necessitate different function handles. Each function
%               provided has a calling form of
%               [xUpdated,PUpdated,innov,S]=h(x,P,z,R);
%               where the first two outputs are the updated state, and the
%               second two are the innovation and innovation covariance
%               matrix. The innovation and innovation covariance matrix are
%               necessary for the computation of likelihoods.
%        muPrev The numModelsX1 set of prior mixing probabilities. 
%        Lambda When using the GPB1, IMM, and GPB2 algorithms, this is the
%               numModelsXnumModels matrix of model transition
%               probabilities. For example, one could use the function call
%               Lambda=getMarkovPTransProbMat(A,T) to get such a transition
%               matrix for a particular model. When using the AMM, this
%               input is not used, and an empty matrix can be passed.
%  numStateDims An optional numModelsX1 vector that lists the number of
%               dimensions in the state in each model. This lets the
%               function be used with models where the states have
%               different dimensionalities. If omitted, it is assumed that
%               all models have a dimensionality of size(xPredSet,1). This
%               parameter should be omitted when using the GPB1 algorithm,
%               as all dynamic models should have the same dimensionality
%               and mixing should be over all elements of the state.
%    numMixDims An optional numModelsX1 or 1XnumModels vector that lists
%               the number of dimensions in each state that are involved in
%               mixing. For example, when designing a model set including a
%               model with a turn rate and a model without a turn rate,
%               usually the turn rate will not be included in the mixing.
%               Similarly, one might mixed models with position and
%               velocity with those including position velocity and
%               acceleration, but if there is also a model with position,
%               velocity and drag, the drag term would not be appropriate
%               to mix with anything else. If omitted,
%               numMixDims=numStateDims is used. This parameter is not used
%               and should be omitted when using the GPB1 algorithm. 
%  
%OUTPUTS: xPostSet The updated set of target states for the models. For all
%                 of the methods except the GPB1, this is xDimMaxXnumModels
%                 in size. For the GBP1, it is xDimMaxX1 as the GPB1 merges
%                 everything into a single hypothesis. Note that the
%                 numModels^2 propagated hypotheses in the GPB2 get merged
%                 into just numModels hypotheses on the output of the
%                 update step.
%        PPostSet The covariance matrices associated with the updated set
%                 of target states. This is xDimMaxXxDimMaxXnumModels in
%                 size, except when using the GPB1, in which case it is
%                 just xDimMaxXxDimMax in size.
%          muPost The numModelsX1 set of updated model probabilities.
%         xMerged The min(numMixDims)X1 state estimate obtained by merging
%                 all of the models according to the method used in the
%                 selected algorithm. This is generally what might be used
%                 for display to an operator.
%         PMerged The min(numMixDims)Xmin(numMixDims) covariance matrix
%                  associated with xMerged. 
%
%Multiple model routines for target tracking are discussed in general in
%Chapter 11.6 of [3]. The IMM tends to be the most widely used multiple
%model routine, while the GPB2 algorithm can be somewhat better, but at a
%higher computational cost. The AMM algorithm does not include model
%transition probabilities (Lambda). This means that the model probabilities
%(muPost) can get very small for models that are not currently in use, thus
%preventing it from ever switching. An ad-hoc fix would be to take the
%output of this function and set some minimum probabilities for models that
%are too unlikely, so that it can still switch in the future.
%
%When using multiple models of different dimensionalities, rather than
%storing the states/ covariances in a cell array, they are stored in
%matrices with the extra elements set to zero. This is much faster in
%Matlab than using a cell array to store states of the exact sizes of the
%models.
%
%Note that with each of the algorithms, this function can be called with an
%empty measurement to cause the model mixing to occur without actually
%performing a measurement update.
%
%REFERENCES:
%[1] J. D. Glass, W. D. Blair, and Y. Bar-Shalom, "IMM estimators with
%    unbiased mixing for tracking targets performing coordinated turns," in
%    Proceedings of the IEEE Aerospace Conference, Big Sky, MT, 2-9 Mar.
%    2013.
%[2] T. Yuan, Y. Bar-Shalom, P. Willett, E. Mozeson, S. Pollak, and D.
%    Hardiman, "A multiple IMM estimation approach with unbiased mixing for
%    thrusting projectiles," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 48, no. 4, pp. 3250-3267, Oct. 2012.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%[4] K. Granstroem, P. Willett, and Y. Bar-Shalom, "Systematic approach
%    to IMM mixing for unequal dimension states," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 51, no. 4, Oct. 2015, pp.
%    2975-2986.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numModels=size(xPredSet,2);

%If the number of dimensions in the states is omitted, assume all the
%states are the same size.
if(nargin<9||isempty(numStateDims))
	xDim=size(xPredSet,1);
	numStateDims=repmat(xDim,[numModels,1]);
end

%If the number of dimensions involved in mixing is omitted, assume that all
%of the dimensions in each vector mix.
if(nargin<10||isempty(numMixDims))
    numMixDims=numStateDims;
end

numMergeDims=min(numMixDims);

%The function is written to assume a cell array was passed, so if the same
%thing is used for all functions, just duplicate it in a cell array.
if(isa(measUpdateFuns,'function_handle'))
    measUpdateFuns=repmat({measUpdateFuns},[numModels,1]);
end

switch(AlgSel)
     case 'AMM'%Autonomous multiple model.
        if(isempty(z))
            %Missed detection event.
            xPostSet=xPredSet;
            PPostSet=PPredSet;
            muPost=muPrev;
        else
            %Perform mode-matched filtering
            [xPostSet,PPostSet,likelihoods]=updateModels(xPredSet,PPredSet,measUpdateFuns,z,R,numStateDims);
            muPost=muPrev.*likelihoods;
            muPost=muPost/sum(muPost);
            %Deal with all of the likelihoods being too small by assuming a
            %uniform distribution if nonfinite numbers arise.
            if(any(~isfinite(muPost)))
                muPost=ones(numModels,1)*(1/numModels);
            end
        end

        %The combined output is just the weighted mixture of the
        %components; no interaction at all between the components is
        %involved in the update of each filter.
        [xMerged,PMerged]=calcMixtureMoments(xPostSet,muPost,PPostSet,[],[],numMergeDims);
     case 'GPB1'%Generalized pseudo-Bayesian 1
        if(isempty(z))
            %If no measurement is provided, but we just want to perform the
            %mixing.
            xPostSet=xPredSet;
            PPostSet=PPredSet;
            %Uninformative likelihoods.
            likelihoods=ones(numModels,1);
        else
            %%Perform mode-matched filtering
            [xPostSet,PPostSet,likelihoods]=updateModels(xPredSet,PPredSet,measUpdateFuns,z,R,numStateDims);
        end
        %%Update the mode probabilities.
        muPost=zeros(numModels,1);
        for curModel=1:numModels
            muPost(curModel)=likelihoods(curModel)*Lambda(:,curModel)'*muPrev;
        end
        muPost=muPost/sum(muPost);

        %Deal with all of the likelihoods being too small by assuming a
        %uniform distribution if nonfinite numbers arise.
        if(any(~isfinite(muPost)))
            muPost=ones(numModels,1)*(1/numModels);
        end
        
        %Combine the estimates
        [xMerged,PMerged]=calcMixtureMoments(xPostSet,muPost,PPostSet,[],[],numMergeDims);

        %The merged output is the only output.
        xPostSet=xMerged;
        PPostSet=PMerged;
    case 'IMM'%Interacting multiple model estimator
        %%Compute the mixing probabilities
        muMerge=bsxfun(@times,Lambda,muPrev(:));
        cBar=zeros(numModels,1);
        for curModel=1:numModels
            %Normalize the row; the normalizing constants get reused
            %during the mode probability update.
            cBar(curModel)=sum(muMerge(:,curModel));
            muMerge(:,curModel)=muMerge(:,curModel)/cBar(curModel);
        end
        
        %%Mix the predicted estimates. The method of [1] is used to make it
        %%unbiased when the states have different dimensionalities.
        [xPredSet,PPredSet]=IMMMixing(muMerge,xPredSet,PPredSet,numStateDims,numMixDims);

        if(isempty(z))
            %If no measurement is provided, but one wanted to do the IMM
            %mixing at this time (a missed detection event).
            xPostSet=xPredSet;
            PPostSet=PPredSet;
            muPost=muPrev;
        else
            %%Perform mode-matched filtering on the mixed estimates.
            [xPostSet,PPostSet,likelihoods]=updateModels(xPredSet,PPredSet,measUpdateFuns,z,R,numStateDims);

            %%Mode probability update
            muPost=cBar.*likelihoods;
            muPost=muPost/sum(muPost);
            %Deal with all of the likelihoods being too small by assuming a
            %uniform distribution if nonfinite numbers arise.
            if(any(~isfinite(muPost)))
                muPost=ones(numModels,1)*(1/numModels);
            end
        end

        %%Combine the estimates
        [xMerged,PMerged]=calcMixtureMoments(xPostSet,muPost,PPostSet,[],[],numMergeDims);
    case 'GPB2'%Generalized Pseudo-Bayesian 2 estimator
        %Unlike the other filters, xPredSet, PPredSet contains predicted
        %states for every pair of prior and current models in the GBP2
        %filter. If model i has n mixing states and model j and m>n mixing
        %states, then the (m-n) extra states are set to those for the
        %transition of model j to model j. The extra covariance elements
        %are also set to those of the j to j transition. This is in an
        %attempt to apply to the GPB2 estimator the same method of making
        %the IMM less biased when handling states of different
        %dimensionalitites.

        if(isempty(z))
            %If no measurement is provided, but we just want to perform the
            %mixing.
            xPostHyps=xPredSet;
            PPostHyps=PPredSet;
            %Uniform likelihoods means uninformative.
            likelihoods=ones(numModels,numModels);
        else
        %%%Perform mode-matched filtering for all combinations of prior and
        %%%posterior models.
            [xPostHyps,PPostHyps,likelihoods]=updateModelsGPB2(xPredSet,PPredSet,measUpdateFuns,z,R,numStateDims);
        end

        %%%Calculation of the merging probabilities.
        muMerge=likelihoods.*Lambda.*repmat(muPrev(:),[1,numModels]);
        cBar=zeros(numModels,1);
        for curModel=1:numModels
            %Normalize the row; the normalizing constants get reused
            %during the mode probability update.
            cBar(curModel)=sum(muMerge(:,curModel));
            muMerge(:,curModel)=muMerge(:,curModel)/cBar(curModel);
            
            %Deal with all of the likelihoods in the row being too small by
            %assuming a uniform distribution if nonfinite numbers arise.
            if(any(~isfinite(muMerge(:,curModel))))
                cBar(curModel)=(1/numModels);
                muMerge(:,curModel)=ones(1,numModels)*cBar(curModel);
            end
        end
        
        %%Merge the hypotheses across all prior models for each posterior
        %model. Because the prediction step copied the appropriate elements
        %to make the transition less biased (in a manner analogous to that
        %used for the IMM), no special type of mixing has to be performed.
        [xPostSet,PPostSet]=GPB2Merging(xPostHyps,PPostHyps,muMerge,numStateDims);
                
        %%%Compute the updated mode probabilities.
        muPost=cBar;
        muPost=muPost/sum(muPost);
        %Deal with all of the likelihoods being too small by assuming a
        %uniform distribution if nonfinite numbers arise.
        if(any(~isfinite(muPost)))
            muPost=ones(numModels,1)*(1/numModels);
        end
        
        %Combine the results into a state estimate and covariance matrix.
        [xMerged,PMerged]=calcMixtureMoments(xPostSet,muPost,PPostSet,[],[],numMergeDims);
    otherwise
        error('Unknown algorithm selected')
end
end

function [xSet,PSet,likelihoods]=updateModels(xSet,PSet,measUpdateFuns,z,R,numStateDims)
%Update all of the numModels hypotheses with the measurements.
    numModels=size(xSet,2);

    likelihoods=zeros(numModels,1);
    for curModel=1:numModels
        idx=1:numStateDims(curModel);
        xPred=xSet(idx,curModel);
        PPred=PSet(idx,idx,curModel);
        h=measUpdateFuns{curModel};
        [xSet(idx,curModel),PSet(idx,idx,curModel),innov,S]=h(xPred,PPred,z,R);
        %Assume a Gaussian likelihood.
        likelihoods(curModel)=(det(2*pi*S))^(-1/2)*exp(-0.5*invSymQuadForm(innov,S));
    end
end

function [xSet,PSet,likelihoods]=updateModelsGPB2(xSet,PSet,measUpdateFuns,z,R,numStateDims)
%Update all of the numModels^2 hypotheses with the measurement.
    numModels=size(xSet,2);
    
    likelihoods=zeros(numModels,numModels);
    for curPrevModel=1:numModels
        for curPostModel=1:numModels
            idx=1:numStateDims(curPostModel);
            xPred=xSet(idx,curPrevModel,curPostModel);
            PPred=PSet(idx,idx,curPrevModel,curPostModel);
            h=measUpdateFuns{curPostModel};
            [xSet(idx,curPrevModel,curPostModel),PSet(idx,idx,curPrevModel,curPostModel),innov,S]=h(xPred,PPred,z,R);
            %Assume a Gaussian likelihood.
            likelihoods(curPrevModel,curPostModel)=(det(2*pi*S))^(-1/2)*exp(-0.5*invSymQuadForm(innov,S));
        end
    end
end

function [xPostSet,PPostSet]=GPB2Merging(xPostHyps,PPostHyps,muMerge,numStateDims)
    xDimMax=size(xPostHyps,1);
    numHyps=size(xPostHyps,2);
    
    %Allocate space
    xPostSet=zeros(xDimMax,numHyps);
    PPostSet=zeros(xDimMax,xDimMax,numHyps);

    %First, do the means.
    for curPostHyp=1:numHyps
        numDims=numStateDims(curPostHyp);
        mixIdx=1:numDims;
        xPostSet(mixIdx,curPostHyp)=sum(bsxfun(@times,muMerge(:,curPostHyp)',xPostHyps(mixIdx,:,curPostHyp)),2);
    end
    
    %Next, do the covariances, taking into account the spread of the means
    %terms.
    for curPostHyp=1:numHyps
       numDims=numStateDims(curPostHyp);
       mixIdx=1:numDims;
        
       for curPreHyp=1:numHyps
           PCur=PPostHyps(mixIdx,mixIdx,curPreHyp,curPostHyp);
           diff=xPostHyps(mixIdx,curPreHyp,curPostHyp)-xPostSet(mixIdx,curPostHyp);
           PPostSet(mixIdx,mixIdx,curPostHyp)=PPostSet(mixIdx,mixIdx,curPostHyp)+muMerge(curPreHyp,curPostHyp)*(PCur+diff*diff');
       end
    end
end

function [xOut,POut]=IMMMixing(muMerge,xSet,PSet,numStateDims,numMixDims)
    numModels=size(xSet,2);
    xDimMax=size(xSet,1);

    %Allocate space for the outputs
    xOut=zeros(xDimMax,numModels);
    POut=zeros(xDimMax,xDimMax,numModels);

    %First, compute the means.
    for curOut=1:numModels
        %The hypothesis for it not changing states. The components of this
        %are used for the non-biased mixing with states with fewer mixing
        %components.
        xOutCur=xSet(:,curOut);
        numOutDim=numStateDims(curOut);
        numOutMixDims=numMixDims(curOut);
        
        idxNoMix=(numOutMixDims+1):numOutDim;
        xOut(idxNoMix,curOut)=xOutCur(idxNoMix);
        
        %The dimensions that are set during the mixing.
        idxAllSet=1:numOutMixDims;
        for curIn=1:numModels
            curInMixDims=numMixDims(curIn);
            dims2Mix=min(numOutMixDims,curInMixDims);
            idxMixDims=1:dims2Mix; 
            
            %Extra dimensions that have to be copied from xOutCur to keep
            %the estimate unbiased when the state being mixed in lacks
            %those dimensions.
            idxExtraDims=(dims2Mix+1):numOutMixDims;
            
            xOut(idxAllSet,curOut)=xOut(idxAllSet,curOut)+muMerge(curIn,curOut)*[xSet(idxMixDims,curIn);xOutCur(idxExtraDims)];
        end
    end

    %Next, compute the covariances, including the unbiased spread of the
    %means terms.
    for curOut=1:numModels
        %The hypothesis is for it not changing states. 
        xOutCur=xSet(:,curOut);
        POutCur=PSet(:,:,curOut);
        numOutDim=numStateDims(curOut);
        numOutMixDims=numMixDims(curOut);
        
        %The dimensions that are set during the mixing.
        idxAllSet=1:numOutMixDims;
        
        %Dimensions that are not involved in the mxing and that just have
        %to be copied from POutCur. The cross terms with these dimensions
        %also have to be set.
        idxNoMix=(numOutMixDims+1):numOutDim;
        POut(idxNoMix,idxNoMix,curOut)=POutCur(idxNoMix,idxNoMix);
        %Now the cross terms, which must be scaled by the probability of
        %the larger hypothesis, or the output matrix runs the risk of not
        %being positive definite.
        POut(idxAllSet,idxNoMix,curOut)=POutCur(idxAllSet,idxNoMix)*muMerge(curOut,curOut);
        POut(idxNoMix,idxAllSet,curOut)=POutCur(idxNoMix,idxAllSet)*muMerge(curOut,curOut);
        
        %The mean term that is used in the spread of the means terms when
        %mixing.
        xMixOut=xOut(idxAllSet,curOut);
        
        %Perform the mixing.
        for curIn=1:numModels
            curInMixDims=numMixDims(curIn);
            dims2Mix=min(numOutMixDims,curInMixDims);
            idxMixDims=1:dims2Mix;
            
            %Extra dimensions that have to be copied from xOutCur  and
            %POutCur to keep the estimate unbiased when the state being
            %mixed in lacks those dimensions.
            idxExtraDims=(dims2Mix+1):numOutMixDims;
            
            PMixCur=zeros(numOutMixDims,numOutMixDims);
            PMixCur(idxMixDims,idxMixDims)=PSet(idxMixDims,idxMixDims,curIn);
            %The non-mixing terms are set to the values in POutCur, but the
            %cross terms are kept all zero.
            PMixCur(idxExtraDims,idxExtraDims)=POutCur(idxExtraDims,idxExtraDims);
            
            diff=[xSet(idxMixDims,curIn);xOutCur(idxExtraDims)]-xMixOut;
            
            POut(idxAllSet,idxAllSet,curOut)=POut(idxAllSet,idxAllSet,curOut)+muMerge(curIn,curOut)*(PMixCur+diff*diff');
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
