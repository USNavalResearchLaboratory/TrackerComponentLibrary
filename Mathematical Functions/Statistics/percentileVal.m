function val=percentileVal(data,p,dim,definition)
%%PERCENTILEVAL Given a set of data, return the values for a particular
%               percentile [0-100]. Two definitions of the percentile are
%               available. 
%
%INPUTS: data A matrix or hypermatrix.
%           p A percentile value from 0 to 100.
%         dim The optional parameter specifying the dimension over which
%             the percentile value is evaluated. The defaults if omitted
%             is to choose the first non-unitary dimension.
%  definition This selects the definition of the percentile that should be
%             used. Possible values are:
%             0 (The default if omitted or an empty matrix is passed) After
%               sorting the data, the value in the percentile of the ith
%               point along the dimension being considered is taken to be
%               100*(i-0.5)/N, where N is the size of dimension dim. Linear
%               interpolation is performed for values that do not align
%               perfectly with a point. Values of p less than 100*(1-0.5)/N
%               are clipped to the first ordered data value and values
%               greater than 100*(N-0.5)/N are clipped to the last ordered
%               data value.
%             1 After sorting the data, the percentile of the ith point
%               along the dimension being considered is 100*(i-1)/(N-1).
%               Thus, the first point is the zeroth percentile and the last
%               is the 100th. Interpolation is performed for points that do
%               not align ith a given percentile.
%
%OUTPUTS: val A matrix having the same dimensions as data except dimension
%             dim is now unitary. The values are the interpolated values of
%             data at the selected percentile.
%
%EXAMPLE:
%Here, we get the 75th percentile of a particular normal distribution using
%both definitions.
% rng('default');%To get the same results each time.
% x=3+1.5*randn(1,10000);
% val0=percentileVal(x,75,[],0)
% val1=percentileVal(x,75,[],1)
%One will get val=3.985566052076216.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(data))
   val=[];
   return;
end

if(nargin<4||isempty(definition))
   definition=0; 
end

if(nargin<3||isempty(dim))
    dim=find(size(data)>1,1);
    if(isempty(dim))
        dim=1; 
    end
end

N=size(data,dim);
p=p/100;%Turn % into fraction.
data=sort(data,dim,'ascend');

switch(definition)
    case 0
        idx=(1/2)+N*p;
        idxMin=floor(idx);
        idxMax=ceil(idx);

        idxVec=repmat({':'},1,ndims(data));

        %Deal with the edge cases.
        if(idxMin<1)
            idxVec{dim}=1;
            val=data(idxVec{:});
            return
        elseif(idxMax>N)
            idxVec{dim}=N;
            val=data(idxVec{:});
            return
        elseif(idxMin==idxMax)
            idxVec{dim}=idxMin;
            val=data(idxVec{:});
            return
        end

        pMin=(idxMin-0.5)/N;
        pMax=(idxMax-0.5)/N;
    case 1
        idx=1+p*(N-1);
        idxMin=floor(idx);
        idxMax=ceil(idx);

        idxVec=repmat({':'},1,ndims(data));
        %Edge cases are not an issue, but we must check for falling exactly
        %on a point.
        if(idxMin==idxMax)
            idxVec{dim}=idxMin;
            val=data(idxVec{:});
            return
        end

        pMin=(idxMin-1)/(N-1);
        pMax=(idxMax-1)/(N-1);
    otherwise
        error('Unknown definition of the percentile specified.')
end

idxVec{dim}=idxMin;
dataMin=data(idxVec{:});
idxVec{dim}=idxMax;
dataMax=data(idxVec{:});

val=dataMin+((p-pMin)/(pMax-pMin))*(dataMax-dataMin);

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
