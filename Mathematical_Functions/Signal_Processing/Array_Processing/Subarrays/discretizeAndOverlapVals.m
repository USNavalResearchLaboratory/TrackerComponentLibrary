function T=discretizeAndOverlapVals(g1,g2,numLevels,adjMat)
%%DISCRETZIEANDOVERLAPVALS Given two length-N vectors of real values g1
%              and g2, discretize the values in each array into numLevels
%              discrete levels each. Then, when overlapping the discretized
%              values, collect elements that have the same combinations of
%              discretized values. The return result is a boolean matrix
%              that selects which of the N elements go into different
%              overlapped discretized subsets. This function is useful for
%              breaking different array configuations into subarrays based
%              on discretizing difference pattern weights as in [1]. The
%              option of an adjacency matrix allows non-adjancent subsets
%              to be separated.
%
%INPUTS: g1, g2 Two NX1 or 1XN sets of values that are to be discretized
%               and overlapped. It is assumes that greater than or equal to
%               numLevels unique values are in g1 and g2. If only g1 is
%               provided, then only discretization is performed and adjMat
%               is not needed.
%     numLevels The number of levels to use for the discretization. If this
%               parameter is omitted or an empty matrix is passed, a
%               default of 4 is used. If g1 or g2 are complex and not just
%               all real or all imaginary, then numLevels can be a 2X1 or
%               1X2 vector where numLevels(1) is the number of levels for
%               the real discretization and numLevels(2) is the number of
%               levels for the complex discretization. The discretization
%               is uniformly spaced across the range of real and imaginary
%               points.
%        adjMat This is an NXN adjacency matrix such that adjMat(i,j) is
%               nonzero if element i is adjacent to element j. The value of
%               adjMat(i,i) does not matter. Only the lower half of this
%               matrix is used. If this matrix is omitted or an empty
%               matrix is passed, all elements are taken to be adjacent to
%               each other. The purpose of this matrix is to make sure that
%               no disjoint subsets are formed. That is, a subset consists
%               of a continuum of adjacent elements with no breaks.
%
%OUTPUTS: T A numSubarraysXN boolean matrix (a matrix such that T(i,:) is a
%           set of boolean values indicating which elements are in each
%           subset). The subsets do not overlap.
%
%The need for the function can be seen from [1]. The functions
%findBaylissSubarrays and findMLSubarrays are examples of how this function
%is used.
%
%REFERENCES:
%[1] U. R. O. Nickel, "Subarray configurations for digital beamforming with
%    low sidelobes and adaptive interference suppression," in Record of the
%    IEEE International Radar Conference, Alexandria, VA, 8-11 May 1995,
%    pp. 714-719.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEls=length(g1);

if(nargin<3||isempty(numLevels))
   numLevels=4; 
end

%numLevels ultimately has to specify the number of levels in the real and
%the imaginary domains. 
if(isscalar(numLevels))
   numLevels=[numLevels;numLevels]; 
end

if(nargin<4||isempty(adjMat))
    %If no adjacency matrix is given, then everything is adjacent.
    adjMat=ones(numEls,numEls);
end

if(isempty(g2))
    gDisc=discComplexValue(g1,numLevels);
    
    diffIdx=find(diff(gDisc));
    numSubarrays=length(diffIdx)+1;

    %The tapering matrix for each subarray. A weight of zero means
    %that the element is not in the subarray.
    T=false(numSubarrays,numEls);
    startIdx=1;
    for curSubarray=1:(numSubarrays-1)
        T(curSubarray,startIdx:(diffIdx(curSubarray)))=true; 
        startIdx=diffIdx(curSubarray)+1;
    end

    %The last subarray.
    T(numSubarrays,startIdx:numEls)=true;
else
    gDisc1=discComplexValue(g1,numLevels);
    gDisc2=discComplexValue(g2,numLevels);

    %The subdivision into subarrays depends on which elements have the same
    %discretization value AND are neighboring. Thus, we need the adjMat matrix
    %to see which elements are next to each other. To determine which elements
    %cluster into subarrays, we will use the DisjointSet data structure.
    theSet=DisjointSet(numEls);

    %If equal non-adjacent values should form different subarrays.
    for curEl1=1:(numEls-1)
        for curEl2=(curEl1+1):numEls
            %If the two elements are adjacent and have the same
            %discretized values (horizontal and vertical.
            if(adjMat(curEl1,curEl2)~=0&&gDisc1(curEl1)==gDisc1(curEl2)&&gDisc2(curEl1)==gDisc2(curEl2))
                theSet.unionFromList([curEl1;curEl2])
            end
        end
    end

    %Now, extract what the subarrays are.
    subarraySet=theSet.createClusterSet();
    numSubarrays=subarraySet.numClusters();

    %The tapering matrix for each subarray. A weight of zero means
    %that the element is not in the subarray.
    T=false(numSubarrays,numEls);
    for curSubarray=1:numSubarrays
        els=subarraySet(curSubarray,:);

        T(curSubarray,els)=true;
    end
end
end

function gDisc=discComplexValue(g1,numLevels)
%%DISCOMPLEXVALUE This function discretizes the real and imaginary parts of
%                 g1 onto a grid with the specified number of levels in
%                 each dimension. The returned gDisc is a set of grid
%                 box numbers.

    g1=g1-min(real(g1));%Make all weights non-negative.
    g1=g1-1j*min(imag(g1));%Make all weights non-negative.
    
    maxValR=max(real(g1(:)));
    maxValI=max(imag(g1(:)));
    cellSizesR=maxValR/numLevels(1);
    cellSizesI=maxValI/numLevels(2);
    
    %Real and complex parts must be separated.
    if(cellSizesR~=0&&cellSizesI~=0)
        %The eps is to make sure that the top value is not aliased into the
        %bottom bin.
        cellSizesR=cellSizesR+eps(cellSizesR);
        cellSizesI=cellSizesI+eps(cellSizesI);
        gDisc=hashPoint2Grid([real(g1(:).');imag(g1(:).')],[cellSizesR;cellSizesI],numLevels);
    elseif(cellSizesR~=0&&cellSizesI==0)
        cellSizesR=cellSizesR+eps(cellSizesR);
        gDisc=hashPoint2Grid(real(g1(:).'),cellSizesR,numLevels(1));
    elseif(cellSizesR==0&&cellSizesI~=0)
        cellSizesI=cellSizesI+eps(cellSizesI);
        gDisc=hashPoint2Grid(imag(g1(:).'),cellSizesI,numLevels(2));
    else
        gDisc=ones(numEls,1);
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
