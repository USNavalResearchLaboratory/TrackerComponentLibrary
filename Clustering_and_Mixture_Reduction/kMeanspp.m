function [mu,P,w,cost]=kMeanspp(z,K,maxIter)
%%KMEANSPP Run the k-means++ algorithm for clustering a set of more than K
%          data points into k clusters. The k-means++ algorithm is the
%          same as the general k-means algorithm, except a different
%          initialization routine is used.
%
%INPUTS: z A zDim X numPoints set of vectors that are to be clustered.
%        K The number of clusters to form. K>=numPoints.
%  maxIter The maximum number of iterations to perform for clustering. If
%          omitted, maxIter=1000;
%
%OUTPUTS: mu A zDim X K set of the cluster means. 
%          P A zDim X zDim X K set of sample covariance matrices for each
%            of the k clusters, where mu(:,n) is the mean of the nth
%            cluster.
%          w A KX1 vector of weights such that w(n) is the fraction of the
%            total original points assigned to the nth cluster.
%       cost The cost of the k-means assignment. This is the sum of the
%            squared distances between the points assigned to a cluster and
%            the cluster mean.
%
%The k-means algorithm is a suboptimal algorithm that tries to find a set
%of k-means such that the sum of the squared distances from the points to
%the closest means is minimized. The implementation here uses a better
%initialization the k-means++ initialization from [1], which improves the
%probability that the algorithm will converge to a good solution.
%
%Note that the k-means algorithm used a randomized initialization, so one
%will not always get the same results when the function is run twice on the
%same data.
%
%REFERENCES:
%[1] D. Arthur and S. Vassilvitskii, "k-means++: The advantages of careful
%    seeding," in Proceedings of the Eighteenth Annual ACM-SIAM Symposium
%    on Discrete Algorithms, New Orleans, LA, Jan. 2007, pp. 1027-1035.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   maxIter=1000; 
end

%This uses the k-mean++ initalization algorithm on the given data.
%K is the number of clusters to form.

zDim=size(z,1);
numPoints=size(z,2);

mu=zeros(zDim,K);

%Use k-mean++ initialization.
%The first center is chosen randomly.
mu(:,1)=z(:,randi([1;numPoints]));
diff=mu(:,1)*ones(1,numPoints)-z;
nearestMuDist=sum(diff.*diff,1);%The squared distances from this point.

for k=2:K
    %The points are chosen with a probability equal to the ratio of the
    %squared distance to the sum of all distances.
    PMF=nearestMuDist/sum(nearestMuDist);%The PMF from which we will draw.
    CMF=cumsum(PMF);
    
    %We have to find the first element in the CMF that is greater than a
    %random draw; that will be the next center.
    idx=sum(CMF<rand(1))+1;
    %idx=randi([1;numPoints]);%For the regular k-means.
    %Set the center and update the nesrest center distances.
    mu(:,k)=z(:,idx);
    diff=mu(:,k)*ones(1,numPoints)-z;
    nearestMuDist=min(nearestMuDist,sum(diff.*diff,1));
end

%The k-means refinement.
muPrev=ones(zDim,K)*Inf;
numIter=0;
while(sum(sum(muPrev~=mu))~=0&&numIter<maxIter)
    %Calculate the distances between the points and the current set of
    %means.
    minCosts=ones(1,numPoints)*Inf;
    selClust=zeros(1,numPoints);
    for k=1:K
        diff=bsxfun(@minus,z,mu(:,k));%mu(:,k)*ones(1,numPoints)-z;
        dist=sum(diff.*diff,1);%The squared distances from this center.
        
        sel=dist<minCosts;
        minCosts(sel)=dist(sel);
        %Assign those points to this cluster if it has the lowest cost.
        selClust(sel)=k;
    end
    
    %Now, all of the points are assigned to clusters, and we can calculate
    %a new set of means.
    muPrev=mu;
    for k=1:K
        mu(:,k)=mean(z(:,selClust==k),2);
    end
    numIter=numIter+1;
end

%The clustering is done. We will now also return the weights and covariance
%estimates. The weights will be the fraction of the total points assigned
%to a given cluster. The covariance will be the sample covariance.
if(nargout>1)
    w=zeros(K,1);
    P=zeros(zDim,zDim,K);
    for k=1:K
        sel=selClust==k;
        w(k)=sum(sel);

        zSel=z(:,sel);

        diff=bsxfun(@minus,zSel,mu(:,k));    
        P(:,:,k)=diff*diff'/w(k);
    end
    w=w/sum(w);

    cost=sum(minCosts)/numPoints;
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
