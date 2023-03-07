function [xClust,xCov,clusterList]=distBasedClustering(x,threshold,distFunc,mergeType,weights,covMats)
%%DISTBASEDCLUSTERING Cluster vectors and merge the cluster (also return
%           cluster information). Two vectors are said to belong to the
%           same cluster if they are sufficiently close together. "Close"
%           is determined by a distance metric specified by distFun. The
%           vectors in a cluster can be merged in a number of different
%           ways.
%
%INPUTS:  x An xDimXN set of N vectors to cluster and possibly merge.
% threshold The threshold for declaring two targets in the same cluster. If
%           their distance is less than this value, then they are declared
%           in the same cluster. The default if this parameter is omitted
%           is (sqrt(2)+sqrt(3))/2, which if the points are given on an
%           integer grid is sufficent to take neighboring points including
%           diagonals (the threshold is larger than needed to avoid finite
%           precision errors making any diagonals not equal).
%  distFunc This specified the distance metric to use to determine how
%           close the points are to each other. If a scalar integer p>=0 is
%           passed, then norm(x(:,i)-x(:,j),p) is used as the distance
%           metric. If the string 'Mahalanobis' is passed, then the input
%           covMats is used and the distance is 
%           (x(:,i)-x(:,j))'*inv(covMats(:,:,i)+covMats(:,:,j))*(x(:,i)-x(:,j))
%           If a matrix or function handle is massed, it is assumed that
%           distFunc(i,j) gives the distance between x(:,i) and x(:,j). It
%           is assumed that distFunc(i,j)=distFunc(j,i). If this parameter
%           is omitted or an empty matrix is passed, then the l2 norm is
%           used.
% mergeType Once the vectors have been clustered, this specifies how the
%           values within each cluster will be merged together. Possible
%           values are:
%           0 Do not perform any merging. xClust and xCov outputs will be
%             empty and only clusterList is returned.
%           1 (The default if this parameter is omitted or an empty matrix
%             is passed) xClust holds the weighted means of the points in
%             each cluster and xCov is the covariances of the points in
%             each cluster.
%           2 The points are assumed to be in a Gaussian mixture with
%             each component having the covariance in covMats. Here, xClust
%             is the mean and xCov the covariance of the submixture that
%             goes into each cluster.
%           3 xClust holds the maximum weighted point in each cluster. xCov
%             hold mean squared error matrices when that point is used as
%             an estimate. If multiple points in x have the same weight,
%             this chooses the first one.
%           4 xClust holds the maximum weighted point in each cluster and
%             xCov is computed like in 2 but instead of the speak being
%             computed with respect to the mean of the mixture, it is
%             computed with respect to the maximum weighted point in each
%             cluster. This effectively provides a mean squared error
%             matrix.
%           5 The points in each cluster are treated as independent
%             Gaussian sampled and are fused. The weights are not used.
%     weights An NX1 or 1XN set of positive weights for each vector. if
%             this parameter is omitted or an empty matrix is passed, then
%             all ones will be used. This is not used if mergeType 5 is
%             chosen.
%     covMats This input is only used if distFunc='Mahalanobis' or
%             mergeType=2 or 4. This is an xDimXxDimXN set of positive
%             definite covariance matrices, one for each vector.
%
%OUTPUTS: xClust An xDimXnumClust set of merged cluster values or an empty
%                matrix if mergeType=0.
%           xCov An xDimXxDimXnumClust set of covariance (or mean squared
%                error) matrices for each of the merged clusters or an
%                empty matrix if mergeType=0.
%    clusterList A ClusterSet class such that clusterList(i,j) gives the
%                index of the vector in x that is the jth vector in the ith
%                cluster.
%
%The determination of which targets are in common clusters is done by
%looking at all pairs of measurements to determine which gate together (
%have distances less than or equal to the threshold). For those that 
%gate together, the gating information is passed to a DisjointSet data
%structure, which is used to produce clusterList.
%
%Given clusterList, the merging depends on mergeType. For mergeType=0, no
%merging is done. For mergeType=1, and  mergeType=2, the merging is done
%using the calcMixtureMoments function, which can handle the covariance
%weighting. For mergeType=3 and 4, the calcMixtureMoments function is used
%to find the mean square error matrix (instead of a covariance matrix),
%where the muHyp input of the function is used with the maximum weight
%point in each cluster for the covariance computation.
%
%EXAMPLE 1:
%Here we make a range-Doppler plot, detect three targets and then centroid
%the detections using this function.
% fc=1e9;%1GHz carrier frequency.
% B=2e6;%2Mhz bandwidth.
% %Baseband start and end frequencies.
% fStart=-B/2;
% fEnd=B/2;
% %Sampling rate is two times the Nyquist rate.
% T0=1/(2*2*fEnd);%Sampling period in seconds.
% T=2e-5;%Chirp duration in seconds.
% 
% PRF=2000;%Pulse repetition frequency (Hertz)
% TB=1/PRF;%The pulse repetition period.
% %The number of samples per PRI. The above parameters were chosen so that
% %this is an integer. Fix just deals with finite precision errors.
% Ns=fix(TB/T0);
% 
% %Generate the reference signal. This is an up-chirp.
% x=LFMChirp(T,fStart,fEnd,T0);
% x=x(:);
% 
% %We will use 64 pulse repetition intervals.
% NB=64;
% 
% %True target parameters.
% c=Constants.speedOfLight;
% rTrue=[40e3;50e3;60e3];
% tau=rTrue/c;%The true delay (s).
% 
% %The true range rate (m/s).
% rrTrue=[100;60;-100];
% a=rrTrue/c;
% 
% %Complex amplitudes
% A=[48;512*exp(1j*2*pi*rand(1));32*exp(1j*2*pi*rand(1))];
% 
% numTargets=length(rTrue);
% 
% %Allocate space for the received signal. The first dimensions is "fast
% %time"; the second dimension if "slow time".
% y=zeros(Ns,NB);
% 
% %Create the received signal, properly delayed and Doppler shifted for
% %each PRI. The same waveform is used in each PRI.
% t=0:T0:((Ns-1)*T0);%Sample times
% for i=0:(NB-1)
%     for curTar=1:numTargets
%         tCur=t-a(curTar)*t-tau(curTar)-i*TB;
%         %The signal simulated with range migration.
%         y(:,i+1)=y(:,i+1)+(A(curTar)*exp(-1j*2*pi*fc*(tau(curTar)+a(curTar)*t)).*LFMChirp(T,fStart,fEnd,tCur)).';
%     end
% 
%     y(:,i+1)=y(:,i+1)+ComplexGaussianD.rand(Ns).';
%     
%     t=t+TB;%Increment to the next time step.
% end
% 
% %Windowing in range and Doppler to lower sidelobes.
% Nx=length(x);
% wRange=windowFunSym(Nx,'Blackman',0);
% x=wRange.*x;
% wDoppler=windowFunSym(NB,'Nuttall',1);
% 
% [delayDopPlot,Doppler,delay]=delayDopplerPlotNBPulseDop(x,y,1,1,wDoppler,T0);
% 
% numGuardCells=[2;4];
% numAvgCells=[5;3];
% PFA=1e-7;
% 
% range=c*delay;
% rangeRate=Doppler*(c/fc);
% 
% DetectionList=CACFAR(delayDopPlot,numGuardCells,numAvgCells,PFA,[],0);
% idx=DetectionList(1).Index;
% 
% rVals=range(idx(1,:));
% RRVals=rangeRate(idx(2,:));
% 
% weights=abs(DetectionList(1).Value);
% 
% centIdxVals=distBasedClustering(idx,[],[],[],weights);
% 
% %Turn the centroided index values into range and range rate values.
% deltaR=range(2)-range(1);
% rValsCent=(centIdxVals(1,:)-1)*deltaR;
% 
% deltaRR=rangeRate(2)-rangeRate(1);
% minVal=rangeRate(1);
% RRValsCent=(centIdxVals(2,:)-1)*deltaRR+minVal;
% 
% %Display the CFAR detections and the centroided detections.
% figure(1)
% clf
% hold on
% scatter(RRVals,rVals,'ob')
% scatter(RRValsCent,rValsCent,'xr','linewidth',2)
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% axis([rangeRate(1), rangeRate(end), range(1), range(end)])
%
%EXAMPLE 2:
%Here, we have two Gaussians. We compare what one obtains using mergeType=1
%versus mergeType=5. With 1, the uncertainty region tends to encompass each
%cluster of points. With 5, it tends to represent an uncertianty region for
%the location of the mean of each Gaussian. The function is not good if the
%Gaussians are close enough that  measurements from Each will overlap (they
%will all just be clustered together due to the primitive clustering
%algorithm).
% x0=[-20;-20];
% R0=[4,2;
%     2,3];
% x1=[20;20];
% R1=[1,0;
%     0,1];
% S0=chol(R0,'lower');
% S1=chol(R1,'lower');
% 
% num0=100;
% num1=200;
% numTotal=num0+num1;
% x=zeros(2,numTotal);
% covMats=zeros(2,2,numTotal);
% 
% x(:,1:num0)=bsxfun(@plus,x0,S0*randn(2,num0));
% covMats(:,:,1:num0)=repmat(R0,[1,1,num0]);
% x(:,(num0+1):numTotal)=bsxfun(@plus,x1,S1*randn(2,num1));
% covMats(:,:,(num0+1):numTotal)=repmat(R1,[1,1,num1]);
% 
% %This is a gating threshold based on a 99.97% confidence region for a 2D
% %Gaussian PDF.
% threshold=ChiSquareD.invCDF(0.9997,2);
% 
% mergeType=1;
% [xClust1,xCov1]=distBasedClustering(x,threshold,[],mergeType,[],covMats);
% mergeType=5;
% [xClust5,xCov5]=distBasedClustering(x,threshold,[],mergeType,[],covMats);
% 
% figure(1)
% clf
% subplot(2,1,1)
% hold on
% scatter(x(1,:),x(2,:),200,'.k')
% scatter(xClust1(1,:),xClust1(2,:),400,'xr')
% numClust=size(xClust1,2);
% for curClust=1:numClust
%     drawEllipse(xClust1(:,curClust),inv(xCov1(:,:,curClust)),[],'-b','linewidth',2)
% end
% title('Merge Type=1')
% 
% subplot(2,1,2)
% hold on
% scatter(x(1,:),x(2,:),200,'.k')
% scatter(xClust5(1,:),xClust5(2,:),400,'xr')
% numClust=size(xClust5,2);
% for curClust=1:numClust
%     drawEllipse(xClust5(:,curClust),inv(xCov5(:,:,curClust)),[],'-b','linewidth',2)
% end
% title('Merge Type=5')
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(x,1);
numTargets=size(x,2);

if(nargin<2||isempty(threshold))
    threshold=(sqrt(2)+sqrt(3))/2;
end

if(nargin<3||isempty(distFunc))
    distFunc=2;%Use the l2 norm.
end

if(nargin<4||isempty(mergeType))
    mergeType=1;
end

if(nargin<5||isempty(weights))
    weights=ones(1,numTargets);
end

if(ischar(distFunc)&&strcmp(distFunc,'Mahalanobis'))
    distFunc=@(i,j)MahalanobisDist(x(:,i),x(:,j),covMats(:,:,i)+covMats(:,:,j));
elseif(isscalar(distFunc))
    %The distance function is a norm.
    p=distFunc;
    distFunc=@(i,j)norm(x(:,i)-x(:,j),p);
end
%If the distance function is anything else, assume that distFunc(i,j) gives
%the distance.

%Create the disjoint set for the clustering.
theSet=DisjointSet(numTargets);

numPoints=size(x,2);
for idx1=1:(numPoints-1)
    for idx2=(idx1+1):numPoints
        if(distFunc(idx1,idx2)<=threshold)
            theSet.unionFromList([idx1;idx2]);
        end
    end
end

%Pull out the clusters
clusterList=theSet.createClusterSet();
theSet.delete();

if(mergeType==0)%No centroiding. Only return the cluster information.
    xClust=[];
    xCov=[];
else
    numClust=clusterList.numClusters();
    xClust=zeros(xDim,numClust);
    xCov=zeros(xDim,xDim,numClust);

    weights=weights(:).';
    switch(mergeType)
        case 1%Weighted mean without using covMats
            for curClust=1:numClust
                idxInClust=clusterList(curClust,:);
                [xClust(:,curClust),xCov(:,:,curClust)]=calcMixtureMoments(x(:,idxInClust),weights(idxInClust)/sum(weights(idxInClust)));
            end
        case 2%Weighted mean using covMats
            weights=weights(:).';
            for curClust=1:numClust
                idxInClust=clusterList(curClust,:);
                [xClust(:,curClust),xCov(:,:,curClust)]=calcMixtureMoments(x(:,idxInClust),weights(idxInClust)/sum(weights(idxInClust)),covMats(:,:,idxInClust));
            end
        case 3%Take the point in the cluster with the maximum weight.
            for curClust=1:numClust
                idxInClust=clusterList(curClust,:);
                [~,maxIdx]=max(weights(idxInClust));
                xClust(:,curClust)=x(:,idxInClust(maxIdx));
                [~,xCov(:,:,curClust)]=calcMixtureMoments(x(:,idxInClust),weights(idxInClust)/sum(weights(idxInClust)),[],xClust(:,curClust));
            end
        case 4%Take the point in the cluster with the maximum weight and
              %use covMats
            for curClust=1:numClust
                idxInClust=clusterList(curClust,:);
                [~,maxIdx]=max(weights(idxInClust));
                xClust(:,curClust)=x(:,idxInClust(maxIdx));
                [~,xCov(:,:,curClust)]=calcMixtureMoments(x(:,idxInClust),weights(idxInClust)/sum(weights(idxInClust)),covMats(:,:,idxInClust),xClust(:,curClust));
            end
        case 5%Treat the points as independent samples and merge them.
               %The weighting is not used.
            for curClust=1:numClust
                idxInClust=clusterList(curClust,:);
                numStates=length(idxInClust);
                xClustCur=x(:,idxInClust(1));
                xCovCur=covMats(:,:,idxInClust(1));
                
                for curMeas=2:numStates
                    [xClustCur,xCovCur]=KalmanUpdate(xClustCur,xCovCur,x(:,idxInClust(curMeas)),covMats(:,:,idxInClust(curMeas)),eye(xDim));
                end

                xClust(:,curClust)=xClustCur;
                xCov(:,:,curClust)=xCovCur;
            end
        otherwise
            error('Unknown centroid type specified')
    end
end
end

function val=MahalanobisDist(x,y,Sigma)
    diff=x-y;
    val=diff'*inv(Sigma)*diff;
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
