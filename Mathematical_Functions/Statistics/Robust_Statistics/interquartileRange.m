function IQR=interquartileRange(x,dim,definition)
%%INTERQUARTILERANGE Determine the interquartile range of the values in x.
%             This is notionally the difference between the value in the
%             75th percentile (the third quartile) and the 25th percentile
%             (the first quartile) of the data. However, the strict
%             definition varies, which is why there are multiple options.
%
%INPUTS: x A vector, matrix or hypermatrix of real values. The
%          interquartile range is only taken over a particular dimension.
%      dim An optional parameter specifying the dimensions across which
%          the interquartile range is found. If this value is omitted or an
%          empty matrix is passed, then the interquartile range is computed
%          across the first non-singleton dimension of x.
% definition This selects the definition fo the interquartile range that is
%          to be used. Let n be the number of points in dimension dim of x.
%          The following definitions are considering each independent
%          vector in dimension dim being processed. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) If an
%            odd number of points is passed, then the first quartile is the
%            median of ordered (ascending) points up to floor(n/2) and the
%            third quartile is the median of ordered points from ceil(n/2)
%            to n. Thus, the median is assigned to the upper half. If an
%            even number of points is passed, then the first quartile is
%            the median of the ordered points up to n/2 and the third
%            quartile is the median of ordered points n/2+1 to n. Thus, the
%            halves are split.
%          1 This is the same as 0, except when given an odd number of
%            points, the middle point is omitted, not assigned to the upper
%            quartile. If only 1 point is given, 0 is returned.
%          2 This is the same as 0 and 1 except when given an odd number of
%            points, the median point is included in both the lower and the
%            upper quartile. This type of interquartile range is also known
%            as the midhinge.
%          3 The first and third quartiles are given by percentileVal with
%            its definition=0.
%          4 The first and third quartiles are given by percentileVal with
%            its definition=1.
%          5 This uses the kthOrderStat to get the 1+floor(0.25*(n-1)) and
%            1+ceil(0.75*(n-1)) order statistics.
%          6 In the set of ordered points, let the split point be
%            k50=ceil(n/2). The lower quartile value is the one at index
%            ceil(k50/2) and the upper quartile value is the one at index
%            k50+ceil((n-k50)/2). This is identical to considering the
%            points part of an empirical distribution and choosing the
%            percentiles using EmpiricalD.invCDF. If none aling directly
%            with 1/4 and 3/4, choose the next highest point.
%
%OUTPUTS: IQR The interquartile range values. This is a matrix with the
%             same dimensions as x except dimension dim has become unitary.
%             If dimension dim is a singleton, then IRQ will be zero for
%             all definitions.
%
%EXAMPLE:
%Here, one sees some of the differences:
% for def=0:6
%     IQR=interquartileRange([1,2,3,4,5],[],def)
% end
% one will get 2.5, 3, 2, 2.5, 2, 2, and 2.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(x))
    IQR=[];
    return;
end

%Select the first non-singleton dimension if dim is not provided.
if(nargin<2||isempty(dim))
    dim=find(size(x)>1,1);
    if(isempty(dim))
        dim=1;
    end
end

if(nargin<3||isempty(definition))
    definition=0;
end

n=size(x,dim);

if(n==1)%If the interquartile range is a singleton, then zero is returned
        %for all methods.
    IQR=zeros(size(x));
    return;
end

switch(definition)
    case 0%The definition used in Mathematica and in Matlab's IQR function.
          %The middle point for an odd number of points is assigned to the
          %third quartile.
        idxVec=repmat({':'},1,ndims(x));

        x=sort(x,dim,'ascend');
        if(mod(n,2)==0)%Split the halves
            idxVec{dim}=1:(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=(n/2+1):n;
            Q3=median(x(idxVec{:}),dim);
        else%Assign the median point to the third quartile.
            idxVec{dim}=1:floor(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=ceil(n/2):n;
            Q3=median(x(idxVec{:}),dim);
        end
    case 1%Alternative successive median definition 1.
        idxVec=repmat({':'},1,ndims(x));
        
        x=sort(x,dim,'ascend');
        if(mod(n,2)==0)%Split the halves
            idxVec{dim}=1:(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=(n/2+1):n;
            Q3=median(x(idxVec{:}),dim);
        else%Omit the median point.
            idxVec{dim}=1:floor(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=(ceil(n/2)+1):n;
            Q3=median(x(idxVec{:}),dim);
        end
    case 2%Alternative successive median definition 2.
        idxVec=repmat({':'},1,ndims(x));
        
        if(mod(n,2)==0)%Split the halves
            idxVec{dim}=1:(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=(n/2+1):n;
            Q3=median(x(idxVec{:}),dim);
        else%Include the median point in both the upper and the lower
            %quartiles.
            idxVec{dim}=1:ceil(n/2);
            Q1=median(x(idxVec{:}),dim);
            
            idxVec{dim}=(ceil(n/2)):n;
            Q3=median(x(idxVec{:}),dim);
        end
        
    case 3%The percentile definition, option 0.
        Q1=percentileVal(x,25,dim,0);
        Q3=percentileVal(x,75,dim,0);
    case 4%The percentile definition, option 1.
        Q1=percentileVal(x,25,dim,1);
        Q3=percentileVal(x,75,dim,1);
    case 5%The order statistic definition.
        [Q1,x]=kthOrderStat(x,1+floor(0.25*(n-1)),dim);
        Q3=kthOrderStat(x,1+ceil(0.75*(n-1)),dim);
    case 6
        idxVec=repmat({':'},1,ndims(x));
        
        x=sort(x,dim,'ascend');
        
        k50=ceil(n/2);
        k25=ceil(k50/2);
        k75=k50+ceil((n-k50)/2);

        idxVec{dim}=k25;
        Q1=x(idxVec{:});
        
        idxVec{dim}=k75;
        Q3=x(idxVec{:});
    otherwise
        error('Unknown definition of the interquartile range specified.')
end

IQR=Q3-Q1;

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
