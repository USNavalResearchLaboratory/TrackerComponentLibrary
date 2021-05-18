function [numMeasPerScan,measMap,invMeasMap]=getCompressedIdxMappings(theHypIdxHist,numSec,numPartsPerSec,secLength)
%%GETCOMPRESSEDIDXMAPPINGS Given a matrix holding a set of hypotheses
%       of assigned measurements at each dimension of an S-dimensional
%       assignment, obtain a mapping of the indices of the measurements
%       at each dimension to one that has no gaps. For example, if in
%       the first dimension, hypotheses contain measurements 1, 3, 4,
%       and 5, these would be remappted to indices 1,2,3, and 4 (the
%       gap between 1 and 3 is eliminated). At the same time, create a
%       matrix to perform the inverse mapping from the compressed
%       indices to the spread out indices. The zero index is assumed to
%       be the missed detection hypothesis and thus is not considered
%       for remapping (it implicitly does not change). Also, this function
%       can handle the case of multiple groups of hypotheses with regions
%       between them that should be ignored.
%
%INPUTS: theHypIdxHist A numDimsXnumHyps collection of numDims-length
%             assignment indices across multiple hypotheses. All of the
%             values must be >=0. If hypotheses come in groups with gaps
%             between them, then the next two inputs must be provided.
% numSec, numPartsPerSec, secLength If these two inputs are provided, it
%             means that theHypIdxHist does not contain all hypotheses in a
%             contiguous manner, but rather that there are numSec sections
%             of length secLength where the ith section contains
%             numPartsPerSec(i) hypotheses. To see how this is used, look
%             at the second example below.
%
%OUTPUTS: numMeasPerScan A numDimsX1 vector where numMeasPerScan(i) is
%             the number of unique indices in theHypIdxHist(i,:). Zero
%             indices present in theHypIdxHist do not count.
%     measMap A max(numMeasPerScan)XnumDims matrix of how compressed
%             indices map to the global indices in theHypIdxHist.
%             measMap(i,j) is how the ith compressed index in the jth
%             dimension maps to a measurement index in jth dimension of
%             theHypIdxHist.
%  invMeasMap This is the inverse map of measMap. This is a
%             maxMeasPerScanXnumDims matrix such that if one takes an
%             index idx from the jth dimension of theHypIdxHist, it
%             maps to the compressed indices as invMeasMap(idx,j),
%             assuming that idx~=0.
%
%Index compression can be useful when using algorithms that like
%sequential indexation of things, but where pruning may have removed a
%large number of indices.
%
%EXAMPLE 1:
%We compress and then decompress some indices where no sections are
%skipped.
% theHypIdxHist=[1,1,2,0;
%                0,2,2,0;
%                1,20,2,0;
%                4,6,3,0];
% [numMeasPerScan,measMap,invMeasMap]=getCompressedIdxMappings(theHypIdxHist);
% %Remap the indices in theHypIdxHist (ignoring zero indices) based on the
% %previded information.
% numDim=size(theHypIdxHist,1);
% numHyp=size(theHypIdxHist,2);
% 
% for curHyp=1:numHyp
%     for curDim=1:numDim
%         origIdx=theHypIdxHist(curDim,curHyp);
%         if(origIdx>0)
%             theHypIdxHist(curDim,curHyp)=invMeasMap(origIdx,curDim);
%         end
%     end
% end
% %View the remapped indices. See that there are no gaps in each dimension.
% theHypIdxHist
% 
% %Map the indices back to their original state.
% for curHyp=1:numHyp
%     for curDim=1:numDim 
%         modIdx=theHypIdxHist(curDim,curHyp);
%         if(modIdx~=0)
%             theHypIdxHist(curDim,curHyp)=measMap(modIdx,curDim);
%         end
%     end
% end
% %See that the original indices have returned.
% theHypIdxHist
%
%%EXAMPLE 2:
%This is an example of a scenario with two targets where we want the
%compressed index mapping. We set the unused parts of each section to NaNs.
%The NaNs are not even looked at by the algorithm and thus do not
%contribute to the unique index determination.
% numSec=2;
% numPartsPerSec=[4;2];
% secLength=5;
% theHypIdxHist=[1,  1, 2, 0, NaN, 1, 0, NaN, NaN, NaN;
%                0,  2, 2, 0, NaN, 6, 1, NaN, NaN, NaN;
%                1, 20, 2, 0, NaN, 3, 2, NaN, NaN, NaN;
%                4,  6, 3, 0, NaN, 9, 3, NaN, NaN, NaN];
% [numMeasPerScan,measMap,invMeasMap]=getCompressedIdxMappings(theHypIdxHist,numSec,numPartsPerSec,secLength);
% 
% %Remap the indices in theHypIdxHist (ignoring zero indices) based on the
% %previded information.
% numDim=size(theHypIdxHist,1);
% numHyp=size(theHypIdxHist,2);
% 
% for curHyp=1:numHyp
%     for curDim=1:numDim
%         origIdx=theHypIdxHist(curDim,curHyp);
%         if(origIdx>0&&~isnan(origIdx))
%             theHypIdxHist(curDim,curHyp)=invMeasMap(origIdx,curDim);
%         end
%     end
% end
% %View the remapped indices. See that there are no gaps in each dimension.
% theHypIdxHist
% 
% %Map the indices back to their original state.
% for curHyp=1:numHyp
%     for curDim=1:numDim 
%         modIdx=theHypIdxHist(curDim,curHyp);
%         if(modIdx~=0&&~isnan(modIdx))
%             theHypIdxHist(curDim,curHyp)=measMap(modIdx,curDim);
%         end
%     end
% end
% %See that the original indices have returned.
% theHypIdxHist 
%
%September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=size(theHypIdxHist,1);
hasSkipSections=(nargin>2&&~isempty(numSec));

%First, determine the maximum measurement index over all scans. This
%will tell us how large an array to assign. Alternatively, we could
%just use a fixed amount.
if(hasSkipSections)
    maxMeasPerScan=0;
    for curSec=1:numSec
        startIdx=1+(curSec-1)*secLength;
        
        for hypIdx=startIdx:(startIdx+numPartsPerSec(curSec)-1)
            maxMeasPerScan=max(maxMeasPerScan,max(theHypIdxHist(:,hypIdx)));
        end
    end
else
    maxMeasPerScan=max(theHypIdxHist(:));
end

%Initially, we use invMeasMat as a boolean matrix, but after that, we
%replace the nonzero elements with appropriate indices.
invMeasMap=zeros(maxMeasPerScan,numDims);
if(hasSkipSections)
    invMeasMap=zeros(maxMeasPerScan,numDims);
    for curDim=1:numDims
        for curSec=1:numSec
            startIdx=1+(curSec-1)*secLength;
            
            for hypIdx=startIdx:(startIdx+numPartsPerSec(curSec)-1)
                idx=theHypIdxHist(curDim,hypIdx);

                if(idx~=0)
                    invMeasMap(idx,curDim)=1;
                end
            end
        end
    end
else
    numHyps=size(theHypIdxHist,2);
    for curDim=1:numDims
        for curHyp=1:numHyps
            idx=theHypIdxHist(curDim,curHyp);

            if(idx~=0)
                invMeasMap(idx,curDim)=1;
            end
        end
    end
end

%We now have a boolean matrix marking which measurement indices are present
%at each scan.
numMeasPerScan=sum(invMeasMap,1);

%Next, we create a matrix such that measMat(i,curDim) will give the index
%of the ith measurement from the specified scan. At the same time, we
%replace the nonzero indices in invMeasMat with the corresponding indices
%in measMat. This lets one map from the larger set of measurement values at
%each scan in theHypIdxHist to a compressed set of values (no gaps).
measMap=zeros(max(numMeasPerScan),numDims);
for curDim=1:numDims
    if(numMeasPerScan(curDim)>0)
        measMap(1:numMeasPerScan(curDim),curDim)=find(invMeasMap(:,curDim));

        for curMeas=1:numMeasPerScan(curDim)
            invMeasMap(measMap(curMeas,curDim),curDim)=curMeas;
        end
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
