function iqRVal=interquartileRangeWeighted(x,w,dim)
%%INTERQUARTILERANGEWEIGHTED Determine the interquartile range of the
%             values in x. This is the difference between the value in the
%             75th percentile and the 25th, when considering them as part
%             of a scalar empirical distribution.
%
%INPUTS: x A vector or a matrix of values.
%        w An optional matrix of weights associated with the samples in x.
%          w should have the same dimensionality as x, be all positive, and
%          should sum to 1 over dimension dim. If omitted or an empty
%          matrix is passed, the samples are assumed to be uniformly
%          weighted.
%      dim An optional parameter specifying the dimensions across which the
%          interquartile range is found. If this value is omitted or an
%          empty matrix is passed, then the interquartile range is computed
%          across the first non-singleton dimension of x.
%
%OUTPUTS: r Values of the interquartiles ranges of the data in x. The
%           dimensionality of r is the same as x, except the dimension over
%           which the ranges are taken is reduced to 1.
%
%The percentiles for the interquartile range are determined using
%EmpiricalD.invCDF. That is, given the points, one build an empirical
%distribution. If none of the discrete points equals 1/4 or 3/4 exactly,
%then choose the next highest point. When considering uniform weights, the
%results of this function are the same as the interquartileRange function
%with definition=6, even though the results are computed in a different
%manner.
%
%EXAMPLES:
% x=[33;-5;42;27;-28;-35;-38;-27;-27;-27;-5];
% iqRVal=interquartileRangeWeighted(x)
%The value obtained is 55. This is=27-(-28). To see why this makes sense,
%consider, the sorted x vector 
%x=[-38;-35;-28;-27;-27;-27;-5;-5;27;33;42];
%and the interquartiles range is 27-(-28).
%On the other hand, if the samples are strongly weighted toward the
%extremes, the result can be widened.
% x=[-38;-35;-28;-27;-27;-27;-5;-5;27;33;42];
% w=[4/10;1/45;1/45;1/45;1/45;1/45;1/45;1/45;1/45;1/45;4/10];
% iqRVal=interquartileRangeWeighted(x,w)
%Now one finds that iqRVal=80 as the very strong weights at either end
%pushed the 25% and 45% values to the extremes.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dimSizes=size(x);
numEls=numel(x);

if(isempty(x))
    iqRVal=[];
    return;
end

%Select the first non-singleton dimension if dim is not provided.
if(nargin<3||isempty(dim))
    dim=find(dimSizes>1,1);
    if(isempty(dim))
        dim=1;
    end
end

numElsSelDim=dimSizes(dim);
dimSizes(dim)=1;

if(nargin<2||isempty(w))
    w=ones(numElsSelDim)/numElsSelDim;
end

%Allocate space
iqRVal=zeros(dimSizes);

%We have to collect numElsSelDim values from dimension dim for every
%possible value of the other dimensions (tuples).
numTuples=numEls/numElsSelDim;
%Elements in each batch of values for dimension dim are separated by
%stride.
stride=prod(dimSizes(1:(dim-1)));
%We thus want to loop through all numTuples of the batches of values with
%consective indices in dimension dim. We want the order of the loop to
%match the order of the values in iqRVal.
numSkips=0;
for curBatch=1:numTuples
    %curBatch corresponds to the ith value in iqRVal. We have to determine
    %the correct initial offset to get to the first item. Subsequent ones
    %are separated by stride.
    offset=curBatch;
    offset=offset-numSkips*stride;
    if(offset>stride)
       numSkips=numSkips+1;
       offset=offset-stride;
    end
    offset=offset+numSkips*stride*numElsSelDim;
    
    idxVals=(0:(numElsSelDim-1))*stride+offset;
    
    iqRVal(curBatch)=diff(EmpiricalD.invCDF([1/4;3/4],x(idxVals),w(idxVals)));
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
