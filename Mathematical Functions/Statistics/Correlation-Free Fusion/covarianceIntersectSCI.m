function [PSCI,xSCI]=covarianceIntersectSCI(CovHyp,covType,xHyp,u,numSamp)
%%COVARIANCEINTERSECTSCI Perform sampling covariance intersection. This is
%                   a method of fusing the first two moments of estimates
%                   when the correlation between the estimates is unknown.
%                   Given the z(inverse) covariance matrices of the
%                   estimates, this function returns a covariance matrix.
%                   Alternatively, if the estimates themselves are given,
%                   this function can also return the merged estimate.
%                   Unlike covarianceIntersect, this function has a random
%                   component and does not specifically optimize any
%                   particular function. Compare this function to
%                   ellipsoidIntersect.
%
%INPUTS: CovHyp The xDimXxDimXN covariance matrices or inverse covariance
%           matrices of the values to be merged.
%   CovType An optional input specifying whether covariance matrices or
%           inverse covariance matrices are in CovHyp. Possible values are:
%           0 (The default if omitted or an empty matrix is passed) CovHyp
%             contains covariance matrices.
%           1 CovHyp contains inverse covariance matrices.
%      xHyp The optional xDimXN set of vectors to merge. These are only
%           needed if xSCI is requested on the output.
%         u A parameter between 0 and 1 that affects the performance of the
%           algorithm. The default if omitted or an empty matrix is passed
%           is 0.5.
%   numSamp The algorithm is stochastic. This optional input is the number
%           of samples to use. The default if this parameter is omitted or
%           an empty matrix is passed is max(100*(N-1),1);
%
%OUTPUTS: PSCI The xDimXxDim fused covariance matrix.
%         xSCI The merged estimate. This requires xHyp to be given on
%              the input. The covariance matrix associated with the
%              merged estimate is PSCI.
%
%This function implements the sampling covariance intersection algorithm of
%[1]. If only a single estimate is passed, then it is just returned. This
%algorithm was developed to fuse tracks from multiple sensors in a
%networked tracking environment.
%
%EXAMPLE:
% This is similar to the example given in the paper for when standard
% covariance intersection might be bad. , except the means are
% not identical.
% x=zeros(2,2);
% P=zeros(2,2,2);
% x(:,1)=[1;-2];
% P(:,:,1)=[1,0;0,100];
% x(:,2)=[-2;-1];
% P(:,:,2)=[100,0;0,1];
% 
% [PM,xM]=covarianceIntersectSCI(P,[],x);
% figure()
% clf
% hold on
% drawEllipse(x(:,1),inv(P(:,:,1)),[],'--r')
% drawEllipse(x(:,2),inv(P(:,:,2)),[],'--g')
% drawEllipse(xM,inv(PM),[],'-b')
%
%REFERENCES:
%[1] X. Tian, Y. Bar-Shalom, and G. Chen, "A no-loss covariance
%    intersection algorithm for track-to-track fusion," in Proceedings of
%    SPIE: Signal and Data Processing of Small targets, vol. 7698,
%    Orlando, FL, 5 Apr. 2010.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(CovHyp,1);
numHyp=size(CovHyp,3);

if(nargin<5||isempty(numSamp))
    numSamp=max(100*(numHyp-1),1);
end

if(nargin<4||isempty(u))
    u=0.5;
end

if(nargin<2||isempty(covType))
   covType=0; 
end

if(covType==0)%If covariance matrices are passed.
    PInvHyp=zeros(numDim,numDim,numHyp);

    for curHyp=1:numHyp
        PInvHyp(:,:,curHyp)=inv(CovHyp(:,:,curHyp));
    end
else%If inverse covariance matrices are passed.
    PInvHyp=CovHyp;
end

if(numHyp>0)
    P0=zeros(numDim,numDim);
    for curHyp=1:numHyp
        P0=P0+PInvHyp(:,:,curHyp);
    end
    P0Inv=P0;
    P0=inv(P0);

    S0=chol(P0,'lower');

    x=S0*randn(numDim,numSamp);

    maxVals=zeros(numSamp,1);
    for curSamp=1:numSamp
        xCur=x(:,curSamp);

        maxVal=-1;
        for curHyp=1:numHyp
            val=xCur'*PInvHyp(:,:,curHyp)*xCur;
            maxVal=max(val,maxVal);
        end
        maxVals(curSamp)=maxVal;
    end

    rMax=-1;
    rMin=Inf;
    for curSamp=1:numSamp
        xCur=x(:,curSamp);
        curRat=xCur'*P0Inv*xCur/maxVals(curSamp);
        rMax=max(curRat,rMax);
        rMin=min(curRat,rMin);
    end

    PSCI=P0./(u*rMin+(1-u)*rMax);
else%Given 1 or 0, don't change it.
    if(numHyp==0)
        xSCI=[];
        PSCI=[];
        return;
    end
    return
end

if(nargin>2&&~isempty(xHyp))
    xSCI=zeros(numDim,1);
    
    for curHyp=1:numHyp
        xSCI=xSCI+PInvHyp(:,:,curHyp)*xHyp(:,curHyp);
    end
    xSCI=P0*xSCI;
else
    xSCI=[];
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
