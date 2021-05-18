function [zDisambig,numDetect]=disambigCRT1D(binMeas,numBins,maxNumBins,twoWayDisambig,threshold,useSparse)
%%DISAMBIGCRT1D This function uses the Chinese remainder theorem to
%               disambiguate measurements in one dimensions (e.g.
%               range or Doppler). It can handle measurements from multiple
%               simultaneous targets. binMeas is a cell array containing
%               measurements (integers from 0 to numBins(i)-1 for the ith
%               snapshot). The snapshots all have a different number of
%               bins; the lengths in numBins should be relatively prime
%               (coprime) integers. When the "true" value of a detection
%               should have been outside of the range of 0->(numBins(i)-1),
%               it is aliased back in (for example, using wrapRange). This
%               function tries to dealias the measurements and determine
%               how many detections are present.
%
%INPUTS: binMeas A numIntX1 or 1XnumInt cell matrix with the detections in
%                each interval. Each cell contains a numMeasX1 or 1XnumMeas
%                vector of the detections in that interval. For the ith
%                interval, all detections must be integers from 0 to
%                numBins(i).
%        numBins The numIntX1 or 1XnumInt vector of the number of bins in
%                each interval. These values must be relatively prime for
%                the algorithm to work properly. A warning will be issued
%                if the values are not relatively prime.
%     maxNumBins This is the total number of bins considered for the
%                completely unaliased measurements. If this value is
%                omitted or an empty matrix is passed, then the default of
%                prod(numBins) is used, which is the maximum from the
%                Chinese remainder theorem. However, as demonstrated below,
%                the presence of multiple targets can cause extra distant
%                "ghost" detections to appear due to the interference of
%                the targets. Thus, it is often good to design the maximum
%                unambiguous value to be very big and to limit the maximum
%                number of bins actually considered. Put another way, use
%                more coprime intervals.
% twoWayDisambig If this is false, then the integers over which aliasing
%                are performed are only positive. Otherwise, they can be
%                positive and negative, which can be more convenient for
%                Doppler disambiguation.
%      threshold An optional threshold for declaring a detection. At least
%                this many common intervals must have values that overlap.
%                The default if this parameter is omitted or an empty
%                matrix is passed is numInt. A detection occurs when a
%                sufficient number of unwrapped measurements align.
%      useSparse Normally the algorithm allocates a vector of length
%                maxNumBins and accumulates values in it. If this parameter
%                is true, then a sparse vector is used. The default if this
%                parameter is omitted or an empty matrix is used is false.
%               
%OUTPUTS: zDisambig Integer values of the detections, a numDetectX1 vector.
%                   If there are no detections, then this will be an empty
%                   matrix. If twoWayDisambig is true, then the values can
%                   be positive and negative.
%         numDetect The number of overlapping values that produced the
%                   detection.
%
%The algorithm is described in [1], among other sources. However, it has a
%number of shortcomings. As given in [2], a clustering algorithm using
%non-binned measurements can sometimes perform better, though the authors
%of [1] would dispute that.
%
%EXAMPLE 1:
%Here, we use the toy example given in [1], but we get the (noise-free)
%measurements by directly aliasing the true unaliased bins to show how that
%works.
%There are three intervals
% numBins=[7;8;11];
% 
% trueBins=[6;13];
% 
% binMeas=cell(3,1);
% for curInt=1:3
%     binMeas{curInt}=unique(wrapRange(trueBins,0,numBins(curInt)));
% end
% 
% disambigCRT1D(binMeas,numBins,[])
%One will get back the true bins of 6 and 13. However, since we did not
%limit the maximum unambiguous range, we also get 237 and 398 due to the
%interference between the two targets. Thus, one either needs to limit the
%maximum range to reduce interference between targets or one needs to use
%more pulse repetition intervals.
%
%EXAMPLE 2:
%This example is better suited for Doppler disambiguation Here, positive
%and negaitve values are present, which could correspond to positive and
%negative range rates. The disambiguation is thus two-sided and the maximum
%number of bins used is limited so that irrelevant extra-fast targets are
%not detected.
% numBins=[7;8;11;13];
% maxNumBins=100;
% %With maxNumBins=100, the range that can be disambiguated is from -50 to
% %+49 when doing it two-sided.
% trueBins=[6;-13];
% 
% binMeas=cell(3,1);
% for curInt=1:4
%     binMeas{curInt}=unique(wrapRange(trueBins,0,numBins(curInt)));
% end
% disambigCRT1D(binMeas,numBins,maxNumBins,true)
%
%REFERENCES:
%[1] P. Stinco, M. Greco, F. Gini, A. Farina, and L. Timmoneri, "Analysis
%    and comparison of two disambiguity algorithms: The modified CA and
%    CRT," in Proceedings of the International Radar Conference -
%    Surveillance for a Safer World, Bordeaux, France, 12-16 Oct. 2009.
%[2] G. Trunk and S. Brockett, "Range and velocity ambiguity resolution,"
%    in Record of the IEEE National Radar Conference, Lynnfield, MA, 20-22
%    Apr. 1993, pp. 146-149.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(areRelativelyPrime(numBins)==false)
   warning('This algorithm may not work as the numbers of bins given are not relatively prime.') 
end

if(nargin<3||isempty(maxNumBins))
     maxNumBins=prod(numBins);
end

if(nargin<4||isempty(twoWayDisambig))
    twoWayDisambig=false;
end

numPRI=length(binMeas);
if(nargin<5||isempty(threshold))
   threshold=numPRI;
end

if(nargin<6||isempty(useSparse))
   useSparse=false;
end

%Bins for all of the aliased-out measurements.
if(useSparse)
    theBins=spalloc(maxNumBins,1,100*numPRI);
else
    theBins=zeros(maxNumBins,1);
end

for curPRI=1:numPRI
    curBins=binMeas{curPRI};
    
    aliasedMax=ceil((maxNumBins-1)/numBins(curPRI));
    if(twoWayDisambig)
        valRange=(-ceil((aliasedMax)/2)):fix((aliasedMax)/2);
    else
        valRange=(0:(aliasedMax-1));
    end
    
    aliasedVals=bsxfun(@plus,curBins(:),numBins(curPRI)*valRange);
    aliasedVals=aliasedVals(:);

    if(twoWayDisambig)
        %Get rid of values that aliased too low.
        aliasedVals=aliasedVals(aliasedVals>-ceil((maxNumBins)/2)-1);
        aliasedVals=wrapRange(aliasedVals,0,maxNumBins);
    else
        %Get rid of values that aliased to indices that are too high. This
        %will only occur in the last ambiguity interval.
        aliasedVals=aliasedVals(aliasedVals<maxNumBins);
    end
    
    %Increment bin values.
    theBins(aliasedVals+1)=theBins(aliasedVals+1)+1;
end

zDisambig=find(theBins>=threshold)-1;
numDetect=theBins(zDisambig+1);
if(twoWayDisambig)    
    zDisambig=wrapRange(zDisambig,-ceil((maxNumBins)/2),fix((maxNumBins)/2));
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
