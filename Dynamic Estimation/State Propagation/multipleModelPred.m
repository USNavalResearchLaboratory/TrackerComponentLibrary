function [xPredSet,PPredSet]=multipleModelPred(AlgSel,xSet,PSet,transFuns,numStateDims,numMixDims)
%%MULTIPLEMODELPRED Given a set of multiple model states (the outputs of
%                   the multipleModelUpdate function can be used),
%                   propagate the estimates forward in time using the
%                   different dynamic models given in transFuns.
%
%INPUTS: algSel A parameter indicating the algorithm to use for updating
%               the models. The same algorithm must be used with this
%               prediction as with the update using the multipleModelUpdate
%               function. Possible values are:
%               'IMM'  The interacting multiple model estimator described
%                      in Chapter 11.6.6 of [3] with the modifications in
%                      [1] and [2] to make it handle multiple models of
%                      different dimensionalities.
%               'AMM'  The autonomous multiple model estimator, which is
%                      the same as the static model estimator of Chapter
%                      11.6.2 of [3].
%               'GPB1' The generalized pseudo-Bayesian estimator 1,
%                      described in Chapter 11.6.4 of [3]. All models have
%                      to have the same dimensionality.
%               'GPB2' The generalized pseudo-Bayesian estimator 2,
%                      described in Chapter 11.6.5 of [3], with changes
%                      similar to that of [1] and [2] to make it work with
%                      models have different dimensionalites.
%          xSet The set of states that are to be predicted for the
%               different models. The models can have different
%               dimensionalities but the number of rows of the matrix is 
%               the dimensionality of the largest model. Models with fewer
%               state elements are zero padded. For all algorithms except
%               the GPB1, xSet is a xDimMaxXnumModels matrix. For the GPB1
%               algorithm, it is just a xDimMaxX1 vector, because the GPB1
%               algorithm merges everything into a single state after each
%               measurement update. Thus, the GPB1 algorithm (unlike the
%               others) requires that all models have the same
%               dimensionality.
%          PSet The set of covariance matrices associated with the states
%               in xSet. This is a xDimMaxXxDimMaxXnumModels matrix for all
%               algorithms except the GPB1, where it is just an
%               xDimMaxXxDimMax matrix.
%     transFuns A numModelsX1 cell array of the state transition functions
%               of the different models. The calling convention for each
%               function in the array is [xPred,PPred]=f(x,P);
%  numStateDims An optional array specifying the number of dimensions in
%               each state. If omitted, it is assumed that all states have
%               the same size, and that size is size(xSet,1); This is not
%               used in the GPB1, because all states have the same size in
%               that algorithm.
%    numMixDims This input is only used in the GPB2 algorithm in the
%               prediction step, though it is used in the update steps for
%               all algorithms except the GPB1. If omitted,
%               numMixDims=numStateDims is used, indicating that all
%               elements in a state mix. This lists the number of
%               dimensions in each state that are involved in mixing.
%               For example, when designing a model set including a model
%               with a turn rate and a model without a turn rate, usually
%               the turn rate will not be included in the mixing.
%               Similarly, one might mixed models with position and
%               velocity with those including position velocity and
%               acceleration, but if there is also a model with position,
%               velocity and drag, the drag term would not be appropriate
%               to mix with anything else.
%
%OUTPUTS: xPredSet The predicted set of target states. For all algorithms
%                 except the GPB2, this is xDimMaxXnumModels in size. When
%                 using the GPB2 algorithm, this is
%                 xDimMaxXnumModelsXnumModels in size as that method makes
%                 all possible model transition hypotheses. These can be
%                 given to multipleModelUpdate with other parameters to
%                 update the model set when measurements become available.
%        PPredSet The predicted set of covariance matrices associated with
%                 the predicted target states. For all algorithms
%                 except the GPB2, this is xDimMaxXxDimMaxXnumModels in
%                 size. When using the GPB2 algorithm, this is
%                 xDimMaxXxDimMaxXnumModelsXnumModels in size.
%
%Multiple model routines for target tracking are discussed in general in
%Chapter 11.6 of [3]. The multiple model algorithms assume that model swaps
%only occur at discrete times during the measurement update. Thus, while
%continuous-time propagation routines can be used here, the model swapping
%method is inherently discrete. One way to handle that would be to set a
%maximum propagation time, and if no measurement arrives in that timespan,
%run multipleModelUpdate with empty matrices for the measurement and its
%covariance and use a measurement update function that just returns the
%prediction. Thus, the uncertainty for model swapping could be better taken
%into account.
%
%When using multiple models of different dimensionalities, rather than
%storing the states/covariances in a cell array, they are stored in
%matrices with the extra elements set to zero. This is much faster in
%Matlab than using a cell array to store states of the exact sizes of the
%models.
%
%The IMM prediction does not perform the mixing of the models. The mixing
%is done at the beginning of the measurement update step. If one wishes to
%perform mixing without a measurement update (which could be done in a
%missed detection case), that can be done after prediction using
%multipleModelUpdate with an empty matrix passed for the measurement.
%
%REFERENCES:
%[1] J. D. Glass, W. D. Blair, and Y. Bar-Shalom, "IMM estimators with
%    unbiased mixing for tracking targets performing coordinated turns," in
%    Proceedings of the IEEE Aerospace Conference, Big Sky, MT, 2-9 Mar.
%    2013.
%[2] T. Yuan, Y. Bar-Shalom, P. Willett, E. Mozeson, S. Pollak, and D.
%    Hardiman, "A multiple IMM estimation approach with unbiased mixing for
%    thrusting projectiles," IEEE Transactions on Aerospace and Electronic
%    Systems, no. 4, pp. 3250-3267, Oct. 2012.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numModels=length(transFuns);

if(nargin<5)
   xDim=size(xSet,1);
   numStateDims=repmat(xDim,[numModels,1]); 
end

if(nargin<6)
    numMixDims=numStateDims;
end

switch(AlgSel)
    case 'AMM'
        %We are given numModels states. Each one is just predicted forward
        %using the appropriate transition function.
        [xPredSet,PPredSet]=predModelsIndependently(xSet,PSet,transFuns,numStateDims);
    case 'GPB1'
        %We are given one state. It is predicted forward in each of the
        %dynamic models. Thus, it is assumed that the size/ meaning of the
        %elements in the state vector is the same for all models.
        xDim=size(xSet,1);
        xPredSet=zeros(xDim,numModels);
        PPredSet=zeros(xDim,xDim,numModels);
        
        for curModel=1:numModels
            f=transFuns{curModel};
            [xPredSet(:,curModel),PPredSet(:,:,curModel)]=f(xSet,PSet);
        end
    case 'IMM'
        %We are given numModels states. Each one is just predicted forward
        %using the appropriate transition function.
        [xPredSet,PPredSet]=predModelsIndependently(xSet,PSet,transFuns,numStateDims);
    case 'GPB2'
        %We are given numModels states. Each state must be predicted
        %forward for all of the dynamic models. This can be difficult if
        %the types of state elements differ between models. A method
        %similar to that used for making the IMM unbiased when handling
        %mixed states is used here. If model i has n mixing states and
        %model j and m>n mixing states, then the (m-n) extra states are set
        %to those for the transition of model j to model j. The extra
        %covariance elements are also set to those of the j to j
        %transition. Basically, if a model lacks components and has to be
        %predicted through a different model, then the components from the
        %actual correct prediction for that model are used in place of
        %zeros or something arbitrary.
        
        xDimTotal=size(xSet,1);
        
        %Allocate space
        xPredSet=zeros(xDimTotal,numModels,numModels);
        PPredSet=zeros(xDimTotal,xDimTotal,numModels,numModels);
        
        for curOut=1:numModels
            numOutDim=numStateDims(curOut);
            numOutMixDims=numMixDims(curOut);
            
            xIn4Out=xSet(:,curOut);
            PIn4Out=PSet(:,:,curOut);
            
            %The state transition function
            f=transFuns{curOut};
            
            %The dimensions that are set during the prediction.
            idxAllSet=1:numOutDim;
            for curIn=1:numModels 
                %We have to construct the state that will be passed to the 
                %state prediction function.

                if(curIn==curOut)
                    %If the model does not switch, then just pass the state
                    %directly to the prediction function.
                    
                    x2Pred=xIn4Out(idxAllSet);
                    P2Pred=PIn4Out(idxAllSet,idxAllSet);
                else
                    %Otherwise, construct a state similar to the unbiased
                    %mixing method in the IMM in [1]. Take the state of the
                    %input and augment it with the parts from the output
                    %model that are not included in the mixing.
                    curInMixDims=numMixDims(curIn);
                    dims2Mix=min(curInMixDims,numOutMixDims);
                    idxMix=1:dims2Mix;
                    idxNoMix=(dims2Mix+1):numOutDim;
                    
                    x2Pred=[xSet(idxMix,curIn);xIn4Out(idxNoMix)];
                    P2Pred=zeros(numOutDim,numOutDim);
                    P2Pred(idxMix,idxMix)=PSet(idxMix,idxMix,curIn);
                    P2Pred(idxNoMix,idxNoMix)=PIn4Out(idxNoMix,idxNoMix);
                    %The cross terms are kept zero.
                end
                
                [xPredSet(idxAllSet,curIn,curOut),PPredSet(idxAllSet,idxAllSet,curIn,curOut)]=f(x2Pred,P2Pred);
            end
        end
    otherwise
        error('Unknown algorithm selected')
end

end

function [xPredSet,PPredSet]=predModelsIndependently(xSet,PSet,transFuns,numStateDims)
    numModels=length(transFuns);
    for curModel=1:numModels
        xDim=numStateDims(curModel);
        idx=1:xDim;
        f=transFuns{curModel};
        [xSet(idx,curModel),PSet(idx,idx,curModel)]=f(xSet(idx,curModel),PSet(idx,idx,curModel));
    end
    xPredSet=xSet;
    PPredSet=PSet;
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
