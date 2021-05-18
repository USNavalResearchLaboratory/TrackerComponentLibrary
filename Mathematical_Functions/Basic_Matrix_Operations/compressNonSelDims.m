function [C,selDims]=compressNonSelDims(C,selDims,justDims)
%%COMPRESSNONSELDIMS Given a hypermatrix, reduce the number of dimensions
%       keeping the selected dimensions as distrinct dimensions. For
%       example, if selDim=3, then a matrix of the form
%       C(:,:,selDim,:,:,:) (where we marked the selected dimensions), will
%       be reduced to the form C(:,selDim,:).
%
%INPUTS: C An S-dimensional matrix or, if justDims=true, size(C).
%  selDims A vector holding the indices of which dimensions of S should be
%          preserved. Alternatively, this could be a length-S boolean
%          vector where selDims(i) specifies whether the ith dimension of C
%          is fixed. If this is a boolean vector, then islogical(selDims)
%          must be true. If this is not a boolean vector, then
%          islogical(selDims) must be false.
% justDims An optional parameter specifying what the first input is. If
%          justDims=true, then instead of passing C, one is passing size(C)
%          and the C output is empty. The default if omitted or an empty
%          matrix is passed is false.
%
%OUTPUTS: C The input matrix reshaped such that non-selected dimensions are
%           compressed as much as possible. Note that the ordering of the
%           elements in C does not change; just the dimensions of C. If
%           justDims=true, then this will be size(C).
%   selDims This indicates which of the compressed dimensions correspond to
%           the original uncompressed dimensions. If selDims on the input
%           is an array of indices, then this will be an array of indices.
%           If selDims on the input was sorted in ascending order, then
%           selDims(i) is the dimension of C to which the ith index
%           originally passed on the input corresponds. If selDims on the
%           input was a boolean vector, then this is a boolean vector
%           indicating which of the compressed indices are fixed.
%
%EXAMPLE:
% C=randn(12,26,36,13,18,9,6);
% selDims=[1;2;5;7];
% [C,newSelDims]=compressNonSelDims(C,selDims);
% size(C)
% newSelDims
%It will compress the third and fourth dimensions. The final size will thus
%be [12, 26, 468, 18, 9, 6] and newSelDims=[1;2;4;6];
%
%November 2020 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(justDims))
    justDims=false; 
end

if(nargin<2)
    selDims=[];
end

%We overwrite nDims with the new dimensions.
if(justDims)
    nDims=C;
    S=length(nDims);
else
    nDims=size(C);
    S=ndims(C);
end

numNewDims=0;
if(islogical(selDims))    
    curDim=1;
    while(curDim<=S)
        if(selDims(curDim))%If it is a fixed dimension.
            numNewDims=numNewDims+1;
            nDims(numNewDims)=nDims(curDim);

            %Overwrite the old selDims.
            selDims(numNewDims)=true;
            curDim=curDim+1;
            continue;
        else
            %If is not a fixed dimension, we have to scan forward until we
            %either find the next fixed dimension or the end.
            prodVal=nDims(curDim);
            curDim=curDim+1;
            while(curDim<=S&&selDims(curDim)==false)
                prodVal=prodVal*nDims(curDim);
                curDim=curDim+1;
            end
            numNewDims=numNewDims+1;
            nDims(numNewDims)=prodVal;
            %This is not a selected dimension.
            selDims(numNewDims)=false;
            continue;
        end
    end
    %Shrink to fit.
    nDims=nDims(1:numNewDims);
    selDims=selDims(1:numNewDims);
else
    numNewDims=0;
    numSel=length(selDims);

    selDims=sort(selDims,'ascend');
    curDim=1;
    for curSelIdx=1:numSel
        curSelDim=selDims(curSelIdx);
        if(curDim<curSelDim)
            %If there is a gap between the current dimension and the
            %selected dimension, everything from curDim to curSelDim-1 gets
            %compressed together.
            prodVal=prod(nDims(curDim:(curSelDim-1)));
            numNewDims=numNewDims+1;
            nDims(numNewDims)=prodVal;
        end
            
        %When here, curDim=curSelDim. Thus, we add the selected dimension
        %and we also record the index in the newSelDims vector.
        numNewDims=numNewDims+1;
        nDims(numNewDims)=nDims(curSelDim);
        %Update the index of the currently selected dimension.
        selDims(curSelIdx)=numNewDims;
        curDim=curSelDim+1;
    end
    
    %If there are dimensions after the last selected dimension, compress
    %them.
    if(curDim<=S)
        prodVal=prod(nDims(curDim:S));
        numNewDims=numNewDims+1;
        nDims(numNewDims)=prodVal;
    end
    
    %Shrink to fit.
    nDims=nDims(1:numNewDims);
end

if(numNewDims==1)
    %Because the nDims input to reshape requires thing to be at least
    %2D.
    nDims=[nDims,1];
end

if(justDims)
    C=nDims;
else
    C=reshape(C,nDims);
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
