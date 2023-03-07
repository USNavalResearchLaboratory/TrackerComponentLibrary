classdef InformationCriteria
%%INFORMATIONCRITERIA A collection of static function for computing
%                     various information criteria as given in [1]. 
%
%Currently implemented criteria are: AIC, AICc, BIC, CLC, AWE, NEC,
%                                    KIC, KICc, AKICc, TIC
%
%REFERENCES:
%[1] S. Akogul and M. Erisoglu, "A Comparison of Information Criteria in
%    Clustering Based on Mixture of Multivariate Normal Distributions,"
%    Mathematical and Computational Applications, vol. 21, no. 3, p. 34,
%    Aug. 2016.
%[2] M. Dixon and T. Ward, "Information-corrected estimation: A 
%    generalization error reducing parameter estimation method,"
%    Entropy, vol. 23, no. 11, p. 1419, 2021.
%
%August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    methods(Static)
        function val = AIC(LL,paramDim)
            %AIC Computes the Akaike information criterion (AIC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+2*paramDim;
        end
        function val = AICc(LL,paramDim,numSamp)
            %%AICc Computes the corrected Akaike information criterion
            %      (AICc). This is recommended over the Akaike
            %      information criterion when paramDim>>numSamp.
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            % numSamp: The number of samples used to compute the
            %          log-likelihood.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+2*paramDim*numSamp/(numSamp-paramDim-1);
        end
        function val = BIC(LL,paramDim,numSamp)
            %%BIC Computes the Bayesian information criterion (BIC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            % numSamp: The number of samples used to compute the
            %          log-likelihood.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+paramDim*log(numSamp);
        end
        function val = CLC(LL,ENtau)
            %%CLC Computes the classification likelihood criterion (CLC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % ENtau: The entropy of the fuzzy classification matrix.
            %        This is defined as:
            %           ENtau = sum(sum(tau_{i,j} * log(tau_{i,j})))
            %        where the first sum is for i = 1:numClusters and the
            %        second sum is for j = 1:numSamples.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %Note that a negative sign was changed from [1] so that the
            %entropy ENtau would be positive.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL-2*ENtau;
        end
        function val = AWE(CLL,paramDim,numSamp)
            %%AWE Computes the approximate weight of evidence criterion
            %     (AWE).
            %
            %INPUTS:
            % CLL: The completed log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            % numSamp: The number of samples used to compute the
            %          log-likelihood.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*CLL+2*paramDim*(3/2+log(numSamp));
        end
        function val = NEC(LL,LL1,ENtau)
            %%NEC Computes the normalized entropy criterion (NEC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % LL1: The log-likelihood value when the number of clusters is
            %      1.
            % ENtau: The entropy of the fuzzy classification matrix.
            %        This is defined as:
            %           ENtau = sum(sum(tau_{i,j} * log(tau_{i,j})))
            %        where the first sum is for i = 1:numClusters and the
            %        second sum is for j = 1:numSamples.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %Note that a negative sign was changed from [1] so that the
            %entropy ENtau would be positive.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -ENtau/(LL-LL1);
        end
        function val = KIC(LL,paramDim)
            %%KIC Computes the Kullback information criterion (KIC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL-3*(paramDim+1);
        end
        function val = KICc(LL,paramDim,numSamp)
            %%KICC Computes the bias-corrected Kullback information
            %      criterion (KICc).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            % numSamp: The number of samples used to compute the
            %          log-likelihood.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+(2*numSamp*(paramDim+1)/(numSamp-paramDim-2))...
                       -numSamp*psi((numSamp-paramDim)/2)+...
                       +numSamp*log(numSamp/2);
        end
        function val = AKICc(LL,paramDim,numSamp)
            %%AKICC Computes the approximate bias-corrected Kullback
            %       information criterion (AICc).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % paramDim: The cardinality of the model's parameter space.
            % numSamp: The number of samples used to compute the
            %          log-likelihood.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+((paramDim+1)*(3*numSamp-paramDim-2)/(numSamp-paramDim-2))...
                       +(paramDim/(numSamp-paramDim));
        end
        function val = TIC(LL,J,K)
            %%TIC Computes the Takeuchi information criterion (TIC).
            %
            %INPUTS:
            % LL: The log-likelihood value.
            % J: The second vector derivative of log(f) with respect to the
            %    parameter vector evaluated at the parameters used to
            %    compute LL. The Fisher information matrix.
            % K: The first vector derivative of log(f) with respect to the
            %    parameter vector evaluated at the parameters used to
            %    compute LL. The expected Hessian.
            %
            %OUTPUTS:
            % val: The criterion value.
            %
            %Note: This equals the AIC if J=K.
            %
            %January 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
            val = -2*LL+2*trace(J\K);
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
