function [rDisambig,sourceInfo]=disambiguateClust1D(rMeas,ampMeas,deltaR,rMax,clustDist,twoWayDisambig,threshold)
%%DISAMBIGUATECLUST1D This function uses a clustering algorithm to
%                   disambiguate measurements in one dimension. It can
%                   handle measurements from multiple simultaneous targets.
%                   Information on which detections went into the
%                   disambiguated measurements is also returned. The
%                   snapshots all have a different number of bins; the
%                   lengths in deltaR should be relatively prime (coprime)
%                   integers (possibly all times a scaling constant). When
%                   the "true" value of a detection is outside of the range
%                   0<=rMeas<deltaR(i) for the ith interval, it is aliased
%                   back in. This function tries to dealias the
%                   measurements and determine how many detections are
%                   present. This function only considers aliasing in a
%                   positive (+deltaR) direction and is thus best suited
%                   for dealiasing ranges.
%
%INPUTS: rMeas A numIntX1 or 1XnumInt cell matrix with the detections in
%              each interval. Each cell contains a numMeasX1 or 1XnumMeas
%              vector of the detections in that interval. Detections should
%      ampMeas A numIntX1 or 1XnumInt cell matrix with the positive real
%              amplitudes associated with the detections in each interval.
%              These are used to perform a weighted average to determine
%              the true target location.
%       deltaR The numIntX1 or 1XnumInt vector of the unambiguous range in
%              each interval. These values must be relatively prime
%              integers (times a common possibly non-integer scale factor)
%              for the algorithm to work.
%         rMax The maximum unambiguous range to consider for the
%              completely unaliased measurements (all dealiased values
%              are less than this. If this value is omitted or an empty
%              matrix is passed, then the default of prod(deltaR) is used,
%              which is the maximum from the Chinese remainder theorem.
%              However, as demonstrated below, the presence of multiple
%              targets can cause extra distant "ghost" detections to appear
%              due to the interference of the targets. Thus, it is often
%              good to design the maximum unambiguous value to be very big
%              and to limit the maximum number of bins actually considered.
%              Put another way, use more coprime intervals.
%    clustDist This is the distance to determine that two aliased-out
%              values are in the same cluster. This depends on the scaling
%              of the intervals. The default if this parameter is omitted
%              or an empty matrix is passed is 0.5, which is good if deltaR
%              are all integers.
% twoWayDisambig If this is false, then the integers over which aliasing
%              are performed are only positive. Otherwise, they can be
%              positive and negative, which can be more convenient for
%              Doppler disambiguation.
%    threshold An optional threshold for declaring a detection. At least
%              this many common intervals must cluster together. The
%              default if this parameter is omitted or an empty matrix is
%              passed is numInt.
%
%OUTPUTS: rDisambig The values of the detections, a numDetectX1 vector.
%              These are weighted averages of the values in rMeas that went
%              into the detections. If there are no detections, then this
%              will be an empty matrix.
%        sourceInfo Information allowing one to reconstruct which values in
%              rMeas went into forming each detection. This is a
%              numDetectX1 cell array. Each cell ion the array has a 
%              numMeasX2 matrix providing information on the measurements
%              in that cluster. The first column gives the number of the
%              interval that produced the measurement (index of rMeas) and
%              the second column provides information on which measurement
%              in that interval led to the detection.
%
%The general idea of such clustering is described in [1] and [2]. However,
%it is done differently here. The authors in [2] described a
%disambiguiation algorithm that is unusable if one of the intervals does
%not produce detections and the basic algorithm in [1] is only focussed on
%a single target. Here, we alias out the measurement across all ambiguity
%intervals and then determine which aliased out measurements gate. Two
%values gate if they are from different intervals and are within clustDist
%of each other. A DisjointSet data structure is used to take the gating
%ifnromation and cluster the aliased out measurements. Clusters of a
%sufficient size are counted as detections. Note that this clustering is
%not globally optimal as it is still possible for two detections from a
%single PRI to end up in one cluster, which is only particularly likely
%with closely spaced targets. An optimal clustering algorithm would be a
%multiframe assignment problem.
%
%Note that because linear averaging is performed, points very near the
%edges will not average to the correct values. However, since measurement
%noises are typically small and unambiguous intervals large, this is seldom
%a problem.
%
%EXAMPLE 1:
%This is similar to the toy example given in [2], except we add some noise
%to the measurements In this instance, the algorithm will produce results
%that are close to the true unaliased bins.
% deltaR=[7;8;11];
% maxR=24;
% 
% trueBins=[6;13];
% 
% rMeas=cell(3,1);
% for curInt=1:3
%     rMeas{curInt}=unique(wrapRange(trueBins,0,deltaR(curInt)));
%     %Add noise.
%     rMeas{curInt}=rMeas{curInt}+0.1*randn(size(rMeas{curInt}));
% end
%
%[rDisambig,sourceInfo]=disambiguateClust1D(rMeas,[],deltaR,maxR)
%
%EXAMPLE 2:
%This example is better suited for Doppler disambiguation Here, positive
%and negaitve values are present, which could correspond to positive and
%negative range rates. The disambiguation is thus two-sided and the maximum
%number of bins used is limited so that irrelevant extra-fast targets are
%not detected.
% deltaR=[7;8;11;13];
% maxR=100;
% %-50 to 49
% trueBins=[6;-13];
% 
% rMeas=cell(3,1);
% for curInt=1:4
%     rMeas{curInt}=unique(wrapRange(trueBins,0,deltaR(curInt)));
%     %Add noise.
%     rMeas{curInt}=rMeas{curInt}+0.1*randn(size(rMeas{curInt}));
% end
% disambiguateClust1D(rMeas,[],deltaR,maxR,[],true)
%
%REFERENCES:
%[1] G. Trunk and S. Brockett, "Range and velocity ambiguity resolution,"
%    in Record of the IEEE National Radar Conference, Lynnfield, MA, 20-22
%    Apr. 1993, pp. 146-149.
%[2] P. Stinco, M. Greco, F. Gini, A. Farina, and L. Timmoneri, "Analysis
%    and comparison of two disambiguity algorithms: The modified CA and
%    CRT," in Proceedings of the International Radar Conference -
%    Surveillance for a Safer World, Bordeaux, France, 12-16 Oct. 2009.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(rMax))
   rMax=prod(deltaR);
end

numPRI=length(rMeas);
if(nargin<5||isempty(clustDist))
   clustDist=0.5;
end

if(nargin<6||isempty(twoWayDisambig))
   twoWayDisambig=false; 
end

if(nargin<7||isempty(threshold))
   threshold=numPRI;
end

%Determine the maximum number of aliased out measurements in both
%directions.
numAliased=0;
for curPRI=1:numPRI
    aliasedMax=ceil(rMax/deltaR(curPRI));
    numAliased=numAliased+length(rMeas{curPRI})*aliasedMax;
end

%We will now alias out all of the measurements
measList=zeros(numAliased,1);
%souceData(i,1) holds the PRI that produced the aliased-out measurement.
%sourceData(i,2) holds the index of the measurement in the PRI.
sourceData=zeros(numAliased,2);
numAdded=0;
for curPRI=1:numPRI
    curMeas=rMeas{curPRI};
    numMeas=length(curMeas);

    aliasedMax=ceil(rMax/deltaR(curPRI));
    if(twoWayDisambig)
        valRange=(-ceil((aliasedMax)/2)):fix((aliasedMax)/2);
    else
        valRange=(0:(aliasedMax-1));
    end
      
    aliasedVals=bsxfun(@plus,curMeas(:),deltaR(curPRI)*valRange);
    origIdx=repmat(1:numMeas,1,aliasedMax+1);
    aliasedVals=aliasedVals(:);
    origIdx=origIdx(:);
    if(twoWayDisambig)
        %Get rid of values that aliased to values that are too low.
        %Get rid of values that aliased too low or too high
        sel=(aliasedVals>=-rMax/2)&(aliasedVals<=rMax/2);
        aliasedVals=aliasedVals(sel);
    else
        %Get rid of values that aliased to values that are too high.
        sel=aliasedVals<=rMax;
        aliasedVals=aliasedVals(sel);
    end
    
    origIdx=origIdx(sel);

    num2Add=length(aliasedVals);
    idx2Add=(numAdded+1):(numAdded+num2Add);
    
    measList(idx2Add)=aliasedVals;
    
    sourceData(idx2Add,1)=curPRI;
    sourceData(idx2Add,2)=origIdx;
    
    numAdded=numAdded+num2Add;
end

%Shrink to fit.
measList=measList(1:numAdded);
sourceData=sourceData(1:numAdded,:);

%Sort the values.
[measList,idx]=sort(measList,'ascend');
sourceData=sourceData(idx,:);

%Create a DisjointSet and determine which measurements gate with each other.
%The DisjointSet will cluster them. To gate, they must be sufficiently
%close and also not within the same interval. The loops go through all
%pairs of aliased-out measurements. We take advantage of having sorted the
%measurements so that once a distance comparison fails, no measurements are
%farther ranges are considered.
gateSet=DisjointSet(numAdded);
for idx1=1:(numAdded-1)
    for idx2=(idx1+1):numAdded
        %If they originate in different PRIs.
        if(sourceData(idx1,1)~=sourceData(idx2,1))
            if(abs(measList(idx1)-measList(idx2))<clustDist)
                %The two values gate together.
                gateSet.unionFromList([idx1;idx2]);
            else
                break;
            end
        end
    end
end

%Form clusters.
theClust=gateSet.createClusterSet();

%Count how many clusteres are large enough to keep.
numClusters=theClust.numClusters();

numPastThresh=sum(theClust.clustSize>=threshold);
rDisambig=zeros(numPastThresh,1);
sourceInfo=cell(numPastThresh,1);

curDetectClust=0;
for curClust=1:numClusters
    numMeasCur=theClust.clustSize(curClust);
    if(numMeasCur>=threshold)
        sourceDataCur=sourceData(theClust(curClust,:),:);
        
        rAmpCur=zeros(numMeasCur,1);
        if(~isempty(ampMeas))
            for curMeas=1:numMeasCur
                PRI=sourceDataCur(curMeas,1);
                idx=sourceDataCur(curMeas,2);
                idxR=theClust(curClust,curMeas);
                
                rAmpCur(curMeas,1)=measList(idxR);
                rAmpCur(curMeas,2)=ampMeas{PRI}(idx);
            end
        else
             for curMeas=1:numMeasCur
                idxR=theClust(curClust,curMeas);
                
                rAmpCur(curMeas,1)=measList(idxR);
                rAmpCur(curMeas,2)=1;
            end
        end

        curDetectClust=curDetectClust+1;
        
        %A simple weighted average.
        rDisambig(curDetectClust)=sum(rAmpCur(:,1).*rAmpCur(:,2))/sum(rAmpCur(:,2));
        sourceInfo{curDetectClust}=sourceDataCur;
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
