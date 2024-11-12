function [mu,P,w,cost,minCostClustPartition,minCostClustSizes,clusterInfo]=kMeanspp(z,K,maxIter,RInv,maxInits)
%%KMEANSPP Run the k-means++ algorithm for clustering a set of more than K
%          data points into k clusters. The k-means++ algorithm is the
%          same as the general k-means algorithm, except a different
%          initialization routine is used. The initialization has also been
%          modified to make use of covariance matrix information, if
%          provided, in which case the iteration run is the kError
%          clustering algorithm of [2].
%
%INPUTS: z A zDimXnumPoints set of vectors that are to be clustered.
%        K The number of clusters to form. K>=numPoints.
%  maxIter The maximum number of iterations to perform for clustering. If
%          omitted, maxIter=1000.
%     RInv If available, this is a zXimXzDimXnumPoints matrix of inverse
%          covariance matrices, one for each point in z. If this matrix is
%          provided (and is not empty), then the kError algorithm of [2] is
%          run instead of the k means algorithm of [1].
% maxInits This is a value >=1 that gives the maximum number of times a new
%          initialization is tried if the previous one iterated until there
%          was an empty cluster. If maxInits initializations have been
%          tried and an empty cluster arises, then the parameters for the
%          previous iteration of k-means are used (which should have only
%          non-empty clusters). The default if omitted or an empty matrix
%          is passed is 3.
%
%OUTPUTS: mu A zDim X K set of the cluster means. 
%          P A zDim X zDim X K set of sample covariance matrices for each
%            of the k clusters, where mu(:,n) is the mean of the nth
%            cluster.
%          w A KX1 vector of weights such that w(n) is the fraction of the
%            total original points assigned to the nth cluster.
%       cost The cost of the k-means assignment. This is the sum of the
%            squared distances between the points assigned to a cluster and
%            the cluster mean or is RInv is provided, this is a sum of
%            Mahalanobis distances.
% minCostClustPartition A numPointsX1 vector holding the index of each
%            cluster (from 1 to numClust) to which each measurement is
%            assigned.
% minCostClustSizes A KX1 vector holding the number of measurements in each
%            cluster.
% clusterInfo This holds the indices of which measurements are in which
%            clusters. clusterInfo(1:minCostClustSizes(k),k) holds the
%            indices of the measurements in the kth cluster.
%
%The k-means algorithm is a suboptimal algorithm that tries to find a set
%of k-means such that the sum of the squared distances from the points to
%the closest means is minimized. The implementation here uses a better
%initialization the k-means++ initialization from [1], which improves the
%probability that the algorithm will converge to a good solution. When RInv
%is provided, a heuristic modification of the k-mean++ algorithm is used
%and the kError algorithm of [1] is then run.
%
%Note that the k-means algorithm used a randomized initialization, so one
%will not always get the same results when the function is run twice on the
%same data.
%
%REFERENCES:
%[1] D. Arthur and S. Vassilvitskii, "k-means++: The advantages of careful
%    seeding," in Proceedings of the Eighteenth Annual ACM-SIAM Symposium
%    on Discrete Algorithms, New Orleans, LA, Jan. 2007, pp. 1027-1035.
%[2] M. Kumar and N. R. Patel, "Clustering data with measurement errors,"
%    Computational Statistics and Data Analysis, vol. 51, pp. 6084-6101,
%    2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(maxIter))
    maxIter=1000; 
end

if(nargin<4)
    RInv=[];
end

if(nargin<5||isempty(maxInits))
    maxInits=3;
end

%This uses the k-mean++ initalization algorithm on the given data.
%K is the number of clusters to form.

zDim=size(z,1);
numPoints=size(z,2);

%The k-mean++ initialization.
mu=zeros(zDim,K);
for curReset=1:maxInits
    shouldReset=false;
    
    %The first center is chosen randomly.
    randIdx=randi([1;numPoints]);
    mu(:,1)=z(:,randIdx);
    diff=bsxfun(@minus,mu(:,1),z);
    if(isempty(RInv))
        nearestMuDist=sum(diff.*diff,1);%The squared distances from this point.
    else
        PInv=zeros(zDim,zDim,K);
        nearestMuDist=zeros(1,numPoints);
        for curPoint=1:numPoints
            nearestMuDist(curPoint)=diff(:,curPoint)'*RInv(:,:,curPoint)*diff(:,curPoint);
        end
    end
    
    for k=2:K
        %The points are chosen with a probability equal to the ratio of the
        %squared distance to the sum of all distances.
        %The PMF from which we will draw.
        PMF=nearestMuDist/sum(nearestMuDist);
        CMF=cumsum(PMF);
        
        %We have to find the first element in the CMF that is greater than
        %a random draw; that will be the next center.
        idx=sum(CMF<rand(1))+1;
        %idx=randi([1;numPoints]);%For the regular k-means.
        %Set the center and update the nearest center distances.
        mu(:,k)=z(:,idx);
        diff=mu(:,k)*ones(1,numPoints)-z;
       if(isempty(RInv))
            nearestMuDist=min(nearestMuDist,sum(diff.*diff,1));
       else
           for curPoint=1:numPoints
               nearestMuDist(curPoint)=min(nearestMuDist(curPoint),diff(:,curPoint)'*RInv(:,:,curPoint)*diff(:,curPoint));
           end
       end
    end
    
    %The k-means refinement.
    muPrev=mu;
    numIter=0;
    while(numIter==0||sum(sum(muPrev~=mu))~=0&&numIter<maxIter)
        %Calculate the distances between the points and the current set of
        %means.
        minCosts=ones(1,numPoints)*Inf;
        selClust=zeros(1,numPoints);
        for k=1:K
            diff=bsxfun(@minus,z,mu(:,k));
            if(isempty(RInv))
                dist=sum(diff.*diff,1);%The squared distances from this center.
            else
                dist=zeros(1,numPoints);
                for curPoint=1:numPoints
                    dist(curPoint)=diff(:,curPoint)'*RInv(:,:,curPoint)*diff(:,curPoint);
                end
            end
            
            sel=dist<minCosts;
    
            if(all(sel==0))
                shouldReset=true;
                break;
            end
    
            minCosts(sel)=dist(sel);
            %Assign those points to this cluster if it has the lowest cost.
            selClust(sel)=k;
        end

        if(shouldReset)
            break;
        end

        %Check for empty clusters. If any are empty, then the algorithm
        %should reset.
        for k=1:K
            if(~any(selClust==k))
                shouldReset=true;
                break;
            end
        end

        if(shouldReset)
            break;
        end

        %Now, all of the points are assigned to clusters, and we can calculate
        %a new set of means.
        muPrev=mu;
        selClustPrev=selClust;
        if(isempty(RInv))
            for k=1:K
                mu(:,k)=mean(z(:,selClust==k),2);
            end
        else
            for k=1:K
                sel=(selClust==k);
                zSel=z(:,sel);
                RInvSel=RInv(:,:,sel);
    
                numMeasSel=size(zSel,2);
                RzMerged=zeros(zDim,1);
                RInvSum=zeros(zDim,zDim,1);
                for curMeas=1:numMeasSel
                    RzMerged=RzMerged+RInvSel(:,:,curMeas)*zSel(:,curMeas);
                    RInvSum=RInvSum+RInvSel(:,:,curMeas);
                end
                mu(:,k)=RInvSum\RzMerged;
                PInv(:,:,k)=RInvSum;
            end
        end
        numIter=numIter+1;
    end
    
    if(shouldReset==false)
        break;
    end
end

if(shouldReset)
    %If it should have reset but ran out of resets, so it couldn't, then
    %just use the previously found solution.
    mu=muPrev;
    selClust=selClustPrev;
end

%The clustering is done. We will now also return the weights and covariance
%estimates. The weights will be the fraction of the total points assigned
%to a given cluster. The covariance will be the sample covariance, unless
%RInv is provided, in which case the covariance matrices of independent
%merged measurements is returned (and will be inconsistent due to
%misassociations).
if(nargout>1)
    w=zeros(K,1);
    P=zeros(zDim,zDim,K);
    minCostClustSizes=zeros(K,1);
    clusterInfo=zeros(numPoints,K);
    if(isempty(RInv))
        for k=1:K
            sel=selClust==k;
            w(k)=sum(sel);
            minCostClustSizes(k)=sum(sel);
            clusterInfo(1:minCostClustSizes(k),k)=find(sel);

            zSel=z(:,sel);
            diff=bsxfun(@minus,zSel,mu(:,k));    
            P(:,:,k)=diff*diff'/w(k);
        end
    else
        for k=1:K
            sel=selClust==k;
            w(k)=sum(sel);
            minCostClustSizes(k)=sum(sel);
            clusterInfo(1:minCostClustSizes(k),k)=find(sel);
            P(:,:,k)=inv(PInv(:,:,k));
        end
    end
    w=w/sum(w);
    cost=sum(minCosts);

    %Size to fit.
    maxNumMeas=max(minCostClustSizes);
    clusterInfo=clusterInfo(1:maxNumMeas,:);

    minCostClustPartition=selClust(:);
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
