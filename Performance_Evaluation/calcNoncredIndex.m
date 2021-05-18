function NCI=calcNoncredIndex(xTrue,xEst,PEst,type,Sigma)
%%CALCNONCREDINDEX Given a batch of estimates and associated true target
%         states, compute one of the variants of the noncredibility index
%         defined in [1].
%
%INPUTS: xTrue The truth data. If the truth is the same for all samples of
%              the estimate, then this is an xDimX1 vector. Otherwise, this
%              is an xDimXN matrix where N is the number of samples in the
%              batch.
%         xEst An xDim X N set of estimates.
%         PEst An xDim X xDim XN set of covariance matrices associated
%              with the estimates. If the covariance matrix is supposed to
%              be the same for all estimates, then PEst is just xDimXxDim.
%              These should not be singular.
%         type An optional parameter selecting the type of noncredibility
%              index to use. Possible values are
%              1 Use the type-1 index defined in Equation 3 in [1]. This
%                requires either that Sigma be provided or that the sample
%                estimate of Sigma be nonsingular.
%              2 (The default if this parameter is omitted or an empty
%                matrix is passed) Use the type-2 index defined in
%                Equations 4 and 6 of [1].
%              3 Use the type-3 index defined in Equations 5 and 7 of [1].
%        Sigma This parameter is only used if type=='I' and is optional.
%              This parameter is an xDimXxDimXN set of matrices of the
%              covariance matrix for the "optimal" estimator. If all of the
%              matrices are the same, then a single xDimXxDim matrix can be
%              passed. If this paraeter is omitted or an empty matrix is
%              passed, then a mean-squared error based on the samples is
%              computed. Matrices in Sigma cannot be singular.  
%
%OUTPUS: NCI The scalar noncredibility index value. Credible estimators
%            have NCI values near 0.
%
%The noncredibility index is discussed in [1] and [2].
%
%EXAMPLE 1:
%Here we compute the noncredibility index of a single-point Cartesian
%initialization algorithm for a given target state estimate. We will see
%that the results are not credible, because the velocity standard deviation
%is just set to cover all uncertainty. We cannot use the type 1 estimator,
%because the sampe Sigma will be singular (resulting in a NaN output.
% xTrue=[100;1000;20;40];
% zTrue=xTrue(1:2);
% sigmaV=300;%To cover the unknown velocity.
% R=[14, 6;
%     6, 14];
% SR=chol(R,'lower');
% numRuns=1000;
% %Generate all Monte Carlo measurements
% zMeas=bsxfun(@plus,zTrue,SR*randn(2,numRuns));
% [xEst,P]=onePointCartInit(zMeas,SR,sigmaV);
% %Turn the covariance matrices into lower-triangular square root covariance
% %matrices.
% S=applyFunToEachMatrix(@(X)chol(X,'lower'),P);
% 
% %Turn lower-triangular square root matrices into regular matrices.
% PEst=applyFunToEachMatrix(@(X)(X*X'),S);
% NCI=calcNoncredIndex(xTrue,xEst,PEst,2)
% NCI=calcNoncredIndex(xTrue,xEst,PEst,3)
% %The non-credibility index will be negative, indicating the
% %inconsistency.
%
%EXAMPLE 2:
%Here we compute the noncredibility index for two-point initialization.
%Here, the results are credible.
% xTrue=[100;1000;20;40];
% zTrue=xTrue(1:2);
% R=[14, 6;
%     6, 14];
% SR=chol(R,'lower');
% T=1;
% F=FPolyKal(T,4,1);
% numRuns=10000;
% %Generate all the Monte Carlo measurements at the first time.
% zMeas(:,1,:)=reshape(bsxfun(@plus,zTrue,SR*randn(2,numRuns)),[2,1,numRuns]);
% %Got to the next time-step
% xTrue=F*xTrue;
% zTrue=xTrue(1:2);
% %Generate all the Monte Carlo measurements at the second time.
% zMeas(:,2,:)=reshape(bsxfun(@plus,zTrue,SR*randn(2,numRuns)),[2,1,numRuns]);
% 
% [xEst,PEst]=twoPointDiffInit(T,zMeas,R);
% 
% NCI1=calcNoncredIndex(xTrue,xEst,PEst,1)
% NCI2=calcNoncredIndex(xTrue,xEst,PEst,2)
% NCI3=calcNoncredIndex(xTrue,xEst,PEst,3)
%The non-credibility indices should all be close to zero indicating
%consistency. However, the third type has a larger variance and is thus
%more likely to deviate from zero.
%
%REFERENCES:
%[1] X. R. Li, Z. Zhao, and V. P. Jilkov, "Estimator's credibility and its
%    measures," in Proceedings of the 15th Triennial Word Congress of the
%    International Federation of Automatic Control, Barcelona, Spain, 21-
%    26 Jul. 2002, pp. 1-6.
%[2] X. R. Li and Z. Zhao, "Measuring estimator's credibility:
%    Noncredibility index," in Proceedings of the 9th International
%    Conference on Information Fusion, Florence, Italy, 10-13 Jul. 2006.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xEst,1);
M=size(xEst,2);

if(size(xTrue,2)==1)
    xTrue=repmat(xTrue,[1,M]);
end

if(size(PEst,3)==1)
    PEst=repmat(PEst,[1,1,M]);
end

if(nargin<4||isempty(type))
   type=2; 
end

switch(type)
    case 1%Equation 3 in [1].
        if(nargin<5||isempty(Sigma))
            
            %Compute Sigma as a covariance matrix
            Sigma=zeros(xDim,xDim);
            for curEst=1:M
                diff=xTrue(:,curEst)-xEst(:,curEst);
                Sigma=Sigma+diff*diff';
            end
            Sigma=Sigma/M;
        end

        if(size(Sigma,3)==1)
            Sigma=repmat(Sigma,[1,1,M]); 
        end

        NCI=0;
        for curEst=1:M
            diff=xTrue(:,curEst)-xEst(:,curEst);
            NCI=NCI+log10(invSymQuadForm(diff,PEst(:,:,curEst)))-log10(invSymQuadForm(diff,Sigma(:,:,curEst)));
        end
        NCI=(10/M)*NCI;
    case 2%Equation 4 with Equation 6
        NCI=0;
        for curEst=1:M
            diff=xTrue(:,curEst)-xEst(:,curEst);
            NCI=NCI+log10(invSymQuadForm(diff,PEst(:,:,curEst)));
        end
        NCI=(10/M)*NCI;
        
        NCI=NCI-(10/log(10))*(log(2)+polygamma(xDim/2));
    case 3%Equation 5 with Equation 7
        NCI=0;
        for curEst=1:M
            diff=xTrue(:,curEst)-xEst(:,curEst);
            NCI=NCI+log10(sum(chol(PEst(:,:,curEst),'lower')\diff)^2);
        end
        NCI=(10/M)*NCI;
        
        gammaVal=-polygamma(1);
        NCI=NCI-(10/log(10))*(log(xDim/2)-gammaVal);
    otherwise
        error('Unknown type of noncredibility index selected')
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
