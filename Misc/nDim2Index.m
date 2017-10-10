function idx=nDim2Index(dims,indices)
%%NDIM2INDEX Given the number of dimensions of a multidimensional array/
%            matrix and a list of indices for each dimension, find the
%            linear index of the point addressed by the indices. This
%            function is like the sub2ind function that comes with Matlab,
%            except the sub2ind function cannot take a vector of indices.
%            Rather, to address n dimensions, it takes n inputs. That is
%            not useful if one wishes to write a function where the number
%            of dimensions needed in a hypermatrix is not known a priori, 
%            such as when handling multivariate polynomials.
%
%INPUTS: dims A 1XnumDim or numDimX1 vector holding the size of each of the
%             nDim dimensions. The values are >=1. Alternatively, a scalar
%             value can be passed if all of the dimensions have the same
%             size.
%     indices A numDim1XnumIdx matrix of indices where each of the
%             corresponding numIdx linear indices is desired. The
%             indexation is Matlab-style, starting from 1, not 0. Normally,
%             numDim1=numDim. However, if numDim1<numDim, then it is
%             assumed that the omitted indices are all ones. This helps
%             deal with singleton dimensions in some instances.
%
%OUTPUTS: idx A numIdxX1 vector of linear indices corresponding to each of
%             the multidimensional index sets given in indices.
%
%EXAMPLE:
% dims=[4;3;2];
% indices=[2;3;1];
% idx0=nDim2Index(dims,indices)
% %The result above is the same as one would get using
% idx1=sub2ind(dims,2,3,1)
%
%The function index2NDim is the inverse of this function.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numDim1=length(indices);
    
    %If all of the dimensions have the same value.
    if(isscalar(dims))
       dims=repmat(dims,[numDim1,1]);
    end
    numDim=numel(dims);
    
    dims=dims(:);%Make it a row vector
    multIdx=[1;cumprod(dims(1:(end-1)))];
    numIdx=size(indices,2);
    
    %Assume omitted trailing indices are 1.
    if(numDim1<numDim)
       indices=[indices;ones(numDim-numDim1,numIdx)];
    end

    idx=ones(1,numIdx);
    for i=1:numDim
        idx=idx+(indices(i,:)-1)*multIdx(i);
    end
    idx=idx';
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
