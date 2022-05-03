function [w,mu,P]=RunnalsGaussMixRed(w,mu,P,K,gammaBound,KMax)
%%RUNNALSGAUSSMIXRED Perform Gaussian mixture reduction using the greedy
%                    merging algorithm by Runnals in [1].
%
%INPUTS: w An NX1 or 1XN vector of weights of the components of the
%          original Gaussian mixture.
%       mu An xDimXN matrix of the means of the vector components of the
%          original Gaussian mixture.
%        P An xDim XxDim XN hypermatrix of the covariance matrices for the
%          components of the original Gaussian mixture.
%        K The minimum desired number of components desired in the mixture
%          after reduction. This will be always achieved unless gammaBound
%          is set. Often, if gammaBound is set, one will just make K=1.
% gammaBound If the number of components is less than or equal to KMax,
%          then reduction will continue until the lowest cost between two
%          components is greater than gammaBound or there are K components
%          left. The default if this parameter is omitted or an empty
%          matrix is passed is Inf, which means that reduction will
%          continue until K components is reached.
%     KMax Regardless of the gammaBound setting, it won't be allowed that
%          more than KMax components are present. The default if omitted or
%          an empty matrix is passed is Inf.
%
%OUTPUTS: w The length K vector of weights of the mixture after reduction.
%           (the length can be less than K if the input mixture had fewer
%           than K components).
%        mu The xDimXK means of the mixture after reduction.
%         P The xDimXxDimXK covariance matrices of the mixture after
%           reduction.
%
%This implements the suboptimal Gaussian mixture reduction algorithm of
%[1]. Care is taken in the implementation to avoid copies and reallocations
%associated with actually resizing the cost array every time something is
%eliminated.
%
%EXAMPLE 1:
%This scalar example plots a full 9-component PDF and then the PDF reduced
%to 5 components using sequential brute-force reduction and Runnals'
%algorithm. Increasing the number of components maintained to 6 makes
%the results of Runnals' algorithm much closer to that of the brute-force
%solution.
% w=[0.03, 0.18, 0.12, 0.19, 0.02, 0.16, 0.06, 0.1, 0.08, 0.06];
% mu=[1.45, 2.20, 0.67, 0.48, 1.49, 0.91, 1.01, 1.42, 2.77, 0.89];
% P=reshape([0.0487, 0.0305, 0.1171, 0.0174, 0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679],[1,1,10]);
% 
% K=6;
% %Runnals' cost function.
% [wRed,muRed,PRed]=RunnalsGaussMixRed(w,mu,P,K,[],[],0);
% %Sequential brute-force reduction.
% [wRedBF,muRedBF,PRedBF]=bruteForceGaussMixRed(w,mu,P,K,true);
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w,mu,P);
% PDFVals2=GaussianMixtureD.PDF(xVals,wRedBF,muRedBF,PRedBF);
% PDFVals3=GaussianMixtureD.PDF(xVals,wRed,muRed,PRed);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'--r','linewidth',2)
% plot(xVals,PDFVals3,'-m','linewidth',1)
% legend('Full Mixture','Brute-Force','Cost Function of Runnalls')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%This is similar to example 1, but we only consider Runnal's cost function
%and we set K=1 and use gammaVal to control how many elements it ultimately
%has in the end. In this example, it reduced to from 10 to KNew=6
%components.
% w=[0.03, 0.18, 0.12, 0.19, 0.02, 0.16, 0.06, 0.1, 0.08, 0.06];
% mu=[1.45, 2.20, 0.67, 0.48, 1.49, 0.91, 1.01, 1.42, 2.77, 0.89];
% P=reshape([0.0487, 0.0305, 0.1171, 0.0174, 0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679],[1,1,10]);
% 
% K=1;
% gammaVal=0.05;
% %Runnals' cost function.
% [wRed,muRed,PRed]=RunnalsGaussMixRed(w,mu,P,K,gammaVal);
% KNew=length(wRed)%The number of reduced components.
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w,mu,P);
% PDFVals2=GaussianMixtureD.PDF(xVals,wRed,muRed,PRed);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'-m','linewidth',1)
% legend('Full Mixture','Reduced Mixture')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Transactions on Aerospace and Electronic Systems,
%    vol. 43, no. 3, pp. 989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(gammaBound))
        gammaBound=Inf;
    end

    if(nargin<6||isempty(KMax))
        KMax=Inf;
    end
    
    N=length(w);
    
    %If no reduction is necessary.
    if(N<=K)
        return;
    end

    %We will only be using one triangular of this matrix.
    M=Inf*ones(N,N);%This is the cost matrix.
    wLogDetPVals=zeros(N,1);
    for k=1:N
        wLogDetPVals(k)=w(k)*log(det(P(:,:,k)));
    end
    
    %We shall fill the cost matrix with the cost of all pairs.
    for cur1=1:(N-1)
        for cur2=(cur1+1):N
            M(cur1,cur2)=BDist(w(cur1),w(cur2),mu(:,cur1),mu(:,cur2),P(:,:,cur1),P(:,:,cur2),wLogDetPVals(cur1),wLogDetPVals(cur2));
        end
    end
    
    KCur=N;
    selIdxPresent=true(N,1);
    for mergeRound=1:(N-K)
        [Mm,minRows]=min(M);%Minimize over the rows
        [minVal,minCol]=min(Mm);%Minimize over the columns.
        minRow=minRows(minCol);
        
        if(minVal>gammaBound&&KCur<=KMax)
            %If the distribution has been sufficiently reduced in terms of
            %cost.
            break;
        end

        %Now we know which two hypotheses to merge: The ones with indices
        %minRow and minCol. We will merge those hypotheses and put the
        %results in minRow.
        wSum=w(minRow)+w(minCol);
        w1=w(minRow)/wSum;
        w2=w(minCol)/wSum;
        muMerged=w1*mu(:,minRow)+w2*mu(:,minCol);
        diff1=mu(:,minRow)-muMerged;
        diff2=mu(:,minCol)-muMerged;
        PMerged=w1*(P(:,:,minRow)+diff1*diff1')+w2*(P(:,:,minCol)+diff2*diff2');
        w(minRow)=wSum;
        mu(:,minRow)=muMerged;
        P(:,:,minRow)=PMerged;
        wLogDetPVals(minRow)=wSum*log(det(PMerged));
        
        %The column is removed.
        selIdxPresent(minCol)=false;
        
        %Make all the costs Inf so that nothing will be assigned to the
        %removed column.
        M(minCol,:)=Inf;
        M(:,minCol)=Inf;
        
        %We must now fill in the costs for the merged estimate, which is in
        %minRow. We shall fill the cost matrix with the cost of all pairs.           
        for cur1=1:(minRow-1)
            if(selIdxPresent(cur1))
                M(cur1,minRow)=BDist(w(cur1),w(minRow),mu(:,cur1),mu(:,minRow),P(:,:,cur1),P(:,:,minRow),wLogDetPVals(cur1),wLogDetPVals(minRow));
            end
        end

        for cur2=(minRow+1):N
            if(selIdxPresent(cur2))
                M(minRow,cur2)=BDist(w(minRow),w(cur2),mu(:,minRow),mu(:,cur2),P(:,:,minRow),P(:,:,cur2),wLogDetPVals(minRow),wLogDetPVals(cur2));
            end
        end

        KCur=KCur-1;
    end
    
    %The merged items are marked with selIdxPresent and can now be all
    %grouped together to be returned.
    w=w(selIdxPresent);
    w=w/sum(w);%Guarantee continued normalization.
    mu=mu(:,selIdxPresent);
    P=P(:,:,selIdxPresent);
end

function val=BDist(w1,w2,mu1,mu2,P1,P2,wLogDetP1,wLogDetP2)
%%BDIST This is the distance measure given in Equation 21 in Section VI-B
%       in [1]. This is a bound related the Kullback-Leibler divergence.
%
%REFERENCES:
%[1] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Transactions on Aerospce and Electronic Systems,
%    vol. 43, no. 3, pp. 989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    wSum=w1+w2;
    w1m=w1/wSum;
    w2m=w2/wSum;
    
    %In [1], there is a constant factor of (1/2) out front, but since all
    %components have it, we omit it here. Also, the dterminant is a slow
    %step, so we put in explicit determinants for things less than 4D.
    switch(size(mu1,1))
        case 1
            diff=mu1-mu2;
            P=w1m*P1+w2m*P2+w1m*w2m*(diff*diff);
            val=wSum*log(P)-wLogDetP1-wLogDetP2;
        case 2
            diff1=mu1(1)-mu2(1);
            diff2=mu1(2)-mu2(2);
            wmProd=w1m*w2m;
            val=wSum*log((w1m*P1(1,1)+w2m*P2(1,1)+wmProd*(diff1*diff1))*(w1m*P1(2,2)+w2m*P2(2,2)+wmProd*(diff2*diff2))-(w1m*P1(1,2)+w2m*P2(1,2)+wmProd*(diff1*diff2)).^2)-wLogDetP1-wLogDetP2;
        case 3
            diff=mu1-mu2;
            P=w1m*P1+w2m*P2+w1m*w2m*(diff*diff');
            detVal=-P(1,3)*P(2,2)*P(3,1)+P(1,2)*P(2,3)*P(3,1)+P(1,3)*P(2,1)*P(3,2)-P(1,1)*P(2,3)*P(3,2)-P(1,2)*P(2,1)*P(3,3)+P(1,1)*P(2,2)*P(3,3);
            val=wSum*log(detVal)-wLogDetP1-wLogDetP2;
        otherwise
            diff=mu1-mu2;
            P12=w1m*P1+w2m*P2+w1m*w2m*(diff*diff');
            val=wSum*log(det(P12))-wLogDetP1-wLogDetP2;
    end
    
    %Deal with the case where w1 and w2 are both essentially zero.
    if(~isfinite(val))
        val=0;
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
