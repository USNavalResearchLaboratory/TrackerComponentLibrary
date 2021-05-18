function minVal=minOverDimIdx(C,selDim,selIdx)
%%MINOVERDIMIDX Given an S-dimensional matrix C, find the minimum value in
%       the matrix when the index in dimension selDim is fixed to selIdx.
%       For example, if C is 5 dimensional and selDim=4this function
%       evaluates the equivalent of min(vec(C(:,:,:,selIdx,:))). This
%       function can be useful when the dimensionality of the matrix in
%       question can vary.
%
%INPUTS: C An S-dimensional hypermatrix.
%   selDim The dimension of the hypermatrix that is selected.
%   selIdx The index in the selected dimension that is fixed.
%
%OUTPUTS: minVal The minimum value in the matrix when the specified
%                dimension is fixed to the specified index.
%
%For an arbitrary-dimensional matrix, the dimensions before the selected
%dimension can be collapsed into a single dimension and those after the
%selected dimension can also be collapsed into a single dimension. Thus,
%the problem reduces to the minimum of a 3D matrix where the middle index
%is fixed. In the even that one selects the first or last dimension, then
%the problem is the minimum of a 2D matrix where the first or last
%dimension is fixed.
%
%EXAMPLE:
%One will see that both of the value obtained by minOverDimIdx in this
%example is the same as the direct Matlab expression.
% C=randn(12,26,36,12,18);
% selDim=3;
% selIdx=8;
% minVal=minOverDimIdx(C,selDim,selIdx)
% min(vec(C(:,:,selIdx,:,:)))
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

nDims=size(C);
S=length(nDims);

if(S==1)
    %The special case of an array.
    minVal=C(selIdx);
elseif(selDim==1)
    startDim=nDims(1);
    endDim=prod(nDims(2:S));
    
    C=reshape(C,[startDim,endDim]);
    minVal=min(C(selIdx,:));
elseif(selDim==S)
    startDim=prod(nDims(1:(S-1)));
    endDim=nDims(S);
    
    C=reshape(C,[startDim,endDim]);
    minVal=min(C(:,selIdx));
else
    startDim=prod(nDims(1:(selDim-1)));
    endDim=prod(nDims((selDim+1):S));
    C=reshape(C,[startDim,nDims(selDim),endDim]);
    minVal=min(vec(C(:,selIdx,:)));
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
