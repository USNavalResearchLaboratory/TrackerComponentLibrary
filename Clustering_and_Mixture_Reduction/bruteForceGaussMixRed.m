function [wRed,muRed,PRed,setParMin]=bruteForceGaussMixRed(w,mu,P,k,sequential)
%%BRUTEFORCEGAUSSMIXRED Perform brute-force Gaussian mixture reduction to
%               minize the integrated squared error (ISE). This function
%               reduces the number of components in the original mixture
%               from n to k by computing the ISE of all reduced
%               distributions obtained by merging components (not by
%               explicitly optimizing over the final values.
%
%INPUTS: w An nX1 or 1Xn vector of weights of the components of the
%          original Gaussian mixture.
%       mu An xDimXn matrix of the means of the vector components of the
%          original Gaussian mixture.
%        P An xDimXxDimXn hypermatrix of the covariance matrices for the
%          components of the original Gaussian mixture.
%        k The number of components desired in the mixture after reduction.
% sequential When reducing the number of components to near half the
%          original value, the brute-force algorithm can be slow. If this
%          parameter is true, instead of reducing the entire distribution
%          all at once, do sequential reductions by 1 component until the
%          number of components has gotten down to k. The default if this
%          parameter is omitted or an empty matrix is passed is false.
%
%OUTPUTS: wRed The kX1 weights of the mixture after reduction.
%        muRed The xDimXk means of the mixture after reduction.
%         PRed The xDimXxDimXk covariance matrices of the mixture after
%              reduction.
%    setParMin The nX1 partition of the original components that produced
%              the reduced distribution. Sets are numbered starting from 1.
%              Components within a common set were merged. If
%              sequential=true, this output is just set to the empty
%              matrix. If n<=k, then this is also an empty matrix.
%
%In the general case, all possible reductions are obtained by going through
%all length-k set partitions using the function getNextSetMPartition. This
%is more efficient than the brute-force algorithm described in [1], which
%generates the length-m set partitions in a very inefficient manner.
%
%An optimization is used to avoid needless recomputation of component costs
%when reducing the number of components by 1. However, in the general case,
%the implementation does not currently precompute merged components nor
%does it precompute PDF values in the ISE cost. However, this could offer a
%substational speed improvement at the expense of a very large increase in
%memory usinge. Most execution time is currently spent in the function
%computeGaussMixISE computing the costs.
%
%EXAMPLE:
%We reduce a scalar PDF from 9 components to 5 components all at once and
%then sequentially and then we plot the results. We also record the
%execution times. Note that the full reduction on many systems might take
%over 10 seconds, but the sequential reduction will typically be a fraction
%of a second.
% w1=[0.03,0.18,0.12,0.19,0.16,0.06,0.1,0.08,0.06];
% w1=w1/sum(w1);
% n1=length(w1);
% mu1=[1.45,2.20,0.67,0.48,0.91,1.01,1.42,2.77,0.89];
% P1=[0.0487,0.0305,0.1171,0.0174,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P1=reshape(P1,[1,1,n1]);
% 
% k=5;
% %Reduction all at once:
% tic
% [wRed,muRed,PRed]=bruteForceGaussMixRed(w1,mu1,P1,k,false);
% toc
% %Sequential reduction.
% tic
% [wRed1,muRed1,PRed1]=bruteForceGaussMixRed(w1,mu1,P1,k,true);
% toc
% 
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w1,mu1,P1);
% PDFVals2=GaussianMixtureD.PDF(xVals,wRed,muRed,PRed);
% PDFVals3=GaussianMixtureD.PDF(xVals,wRed1,muRed1,PRed1);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'--r','linewidth',2)
% plot(xVals,PDFVals3,'-g','linewidth',1)
% legend('Full Mixture','Full Reduction','Sequential Reduction')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the reduced PDFs are both essentially the same However,
%if the number of components in the reduced Gaussian were increased to 6,
%then we would see a larger difference between the methods.
%
%REFERENCES:
%[1] D. F. Crouse, "A look at Gaussian mixture reduction algorithms," in
%    Proceedings of the 14th International Conference on Information
%    Fusion, Chicago, IL, 5-8 Jul. 2011.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(sequential))
    sequential=false;
end

n=length(w);

if(n<=k)
    %If no reduction is necessary.
    wRed=w;
    muRed=mu;
    PRed=P;
    setParMin=[];
   return 
end


if(sequential==false)
    if(k==n-1)%If it is just reducing it by one component.
        [wRed,muRed,PRed,setParMin]=bruteForceGaussMixRedBy1(w,mu,P);
    else
        [wRed,muRed,PRed,setParMin]=bruteForceGaussMixRedFull(w,mu,P,k);
    end
else
    wRed=w;
    muRed=mu;
    PRed=P;
    %The special case.
    if(k==n)
        return;
    end
    
    for curComp=(n-1):-1:k
        [wRed,muRed,PRed]=bruteForceGaussMixRedBy1(wRed,muRed,PRed);
    end
    setParMin=[];
end
end

function [wRed,muRed,PRed,setParMin]=bruteForceGaussMixRedFull(w,mu,P,k)
%%BRUTEFORCEGAUSSMIXREDFULL The full brute-force Gaussian mixture reduction
%           algorithm without any particular optimization relating to how k
%           is related to n, the length of w.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(w);
xDim=size(mu,1);

%Special cases
if(n==k)
    wRed=w;
    muRed=mu;
    PRed=P;
    return
elseif(k==1)
    wRed=1;
    [muRed,PRed]=calcMixtureMoments(mu,w,P);
    return
end

[curSetPart,recurSetData]=getNextSetMPartition(n,k);

costMin=Inf;

wRed=zeros(k,1);
muRed=zeros(xDim,k);
PRed=zeros(xDim,xDim,k);
while(~isempty(curSetPart))
    for setIdx=1:k
        sel=(curSetPart==setIdx);
        wSelNorm=w(sel);
        wSum=sum(wSelNorm);
        wSelNorm=wSelNorm/wSum;
        
        [muMerged,PMerged]=calcMixtureMoments(mu(:,sel),wSelNorm,P(:,:,sel));
        
        wRed(setIdx)=wSum;
        muRed(:,setIdx)=muMerged;
        PRed(:,:,setIdx)=PMerged;
    end
    costVal=computeGaussMixISE(w,mu,P,wRed,muRed,PRed,1);
    
    if(costVal<costMin)
        costMin=costVal;
        setParMin=curSetPart;
    end

    [curSetPart,recurSetData]=getNextSetMPartition(recurSetData);
end

%Return the best reduced distribution.
for setIdx=1:k
    sel=(setParMin==setIdx);
    wSelNorm=w(sel);
    wSum=sum(wSelNorm);
    wSelNorm=wSelNorm/wSum;
    
    [muMerged,PMerged]=calcMixtureMoments(mu(:,sel),wSelNorm,P(:,:,sel));

    wRed(setIdx)=wSum;
    muRed(:,setIdx)=muMerged;
    PRed(:,:,setIdx)=PMerged;
end
end

function [wRed,muRed,PRed,setParMin]=bruteForceGaussMixRedBy1(w,mu,P)
%%BRUTEFORCEGAUSSMIXREDBY1 Brute-force Gaussian mixture redution by 1
%       component using the ISE cost function avoiding needless
%       recomputation of the PDF values.
%
%To understand the optimization, look at the computation of various terms
%incomputeGaussMixISE.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(w);
xDim=size(mu,1);

%Get all cross PDF values of the unmerged components and themselves. We
%will only use the upper-triangular part of the matrix.
PDFValsUnmerged=zeros(n,n);
for i1=1:n
    %The following line is equivalent to 
    %w(i1)^2*GaussianD.PDF(mu2(:,n1),mu2(:,n1),2*P2(:,:,n1));
    %because the argument of the exponent is zero.
    PDFValsUnmerged(i1,i1)=(det(4*pi*P(:,:,i1)))^(-1/2);
    for i2=(i1+1):n
        PDFValsUnmerged(i1,i2)=GaussianD.PDF(mu(:,i1),mu(:,i2),P(:,:,i1)+P(:,:,i2));
        PDFValsUnmerged(i2,i1)=PDFValsUnmerged(i1,i2);
    end
end

ISEMin=Inf;

%To hold temporary values and then for the return weights.
wRed=zeros(n-1,1);
PDFValsCross=zeros(n,1);

%Now, go through all possibilities of which two elements are merged. 
for i1=1:n
    for i2=(i1+1):n
        %Get the current merged components.
        wMerged=w(i1)+w(i2);
        muMerged=(w(i1)*mu(:,i1)+w(i2)*mu(:,i2))/wMerged;
        diff=mu(:,i1)-muMerged;
        PMerged=w(i1)*(P(:,:,i1)+diff*diff');
        diff=mu(:,i2)-muMerged;
        PMerged=PMerged+w(i2)*(P(:,:,i2)+diff*diff');
        PMerged=PMerged/wMerged;

        %Get the full set of weights of the reduced distribution.
        curIdx=1;
        for i=1:n
            if(i~=i1&&i~=i2)
                wRed(curIdx)=w(i);
                curIdx=curIdx+1;
            end
        end
        wRed(n-1)=wMerged;

        for k1=1:n
            PDFValsCross(k1)=GaussianD.PDF(mu(:,k1),muMerged,P(:,:,k1)+PMerged);
        end
        %The simplified gaussian PDF with itself.
        PDFValMerged=(det(4*pi*PMerged))^(-1/2);
  
        Jhr=0;
        for k1=1:n
            %The index in the reduced distribution.
            curIdx=1;
            for k2=1:n
                if(k2==i1||k2==i2)
                    continue;
                end 
                Jhr=Jhr+w(k1)*wRed(curIdx)*PDFValsUnmerged(k1,k2);
                curIdx=curIdx+1;
            end
            %Add in the cross term for the merged component
            Jhr=Jhr+w(k1)*wRed(curIdx)*PDFValsCross(k1);
        end
        
        Jrr=0;
        curIdx1=0;
        for k1=1:n
            if(k1==i1||k1==i2)
                continue;
            end
            curIdx1=curIdx1+1;
            
            Jrr=Jrr+wRed(curIdx1)^2*PDFValsUnmerged(k1,k1);
            
            curIdx2=curIdx1;
            for k2=(k1+1):n
                if(k2==i1||k2==i2)
                    continue;
                end
                curIdx2=curIdx2+1;
                
                Jrr=Jrr+2*wRed(curIdx1)*wRed(curIdx2)*PDFValsUnmerged(k1,k2);
            end
            %Add in the final term.
            Jrr=Jrr+2*wRed(curIdx1)*wRed(n-1)*PDFValsCross(k1);
        end
        %Add in the final term with itself.
        Jrr=Jrr+wRed(n-1)^2*PDFValMerged;
        
        %The ISE omitting the constant term.
        ISEVal=Jrr-2*Jhr;
        
        if(ISEVal<ISEMin)
            ISEMin=ISEVal;
            minIdx1=i1;
            minIdx2=i2;
        end
    end
end

muRed=zeros(xDim,n-1);
PRed=zeros(xDim,xDim,n-1);

%Get the optimal merged components.
wMerged=w(minIdx1)+w(minIdx2);
muMerged=(w(minIdx1)*mu(:,minIdx1)+w(minIdx2)*mu(:,minIdx2))/wMerged;
diff=mu(:,minIdx1)-muMerged;
PMerged=w(minIdx1)*(P(:,:,minIdx1)+diff*diff');
diff=mu(:,minIdx2)-muMerged;
PMerged=PMerged+w(minIdx2)*(P(:,:,minIdx2)+diff*diff');
PMerged=PMerged/wMerged;

%Get the reduced Gaussian mixture with the merged components at the end.
curIdx=1;
for i=1:n
    if(i~=minIdx1&&i~=minIdx2)
        wRed(curIdx)=w(i);
        muRed(:,curIdx)=mu(:,i);
        PRed(:,:,curIdx)=P(:,:,i);
        curIdx=curIdx+1;
    end
end
wRed(n-1)=wMerged;
muRed(:,n-1)=muMerged;
PRed(:,:,n-1)=PMerged;

if(nargout>3)
%If desired, return a set partition indicating which compoinents were
%merged. The merged components ar eput in set 1 and the other components
%are assigned their own sets.
    setParMin=zeros(n,1);
    
    setParMin(minIdx1)=1;
    setParMin(minIdx2)=1;
    
    curIdx=2;
    for k=1:n
        if(k==minIdx1||k==minIdx2)
            continue;
        end
        setParMin(k)=curIdx;
        curIdx=curIdx+1;
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
