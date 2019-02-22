function val=geometricMean(x,dim)
%%GEOMETRICMEAN Compute the geometric mean of a set of samples.
%
%INPUTS: x The set of samples over which the geometric mean is to be found.
%          How the function handles vectors and matrices is determined by
%          the parameter dim.
%      dim An optional parameter specifying over which dimension of x the
%          geometric mean is to be taken. If dim is omitted or an empty
%          matrix is passsed, then if x is a vector, the geometric mean is
%          found over all elements of the vector. If x is a matrix, then
%          the mean is found over the columns of the matrix (resulting in a
%          row vector), and if x is an n-dimensional (n>2) matrix, then the
%          geometric mean is found over the first non-singleton dimension
%          of the matrix.
%
%OUTPUT: val The geometric mean of x taken over the appropriate dimension.
%
%The geometric mean of n items is the nth-root of the product of those
%items. The computation is performed using logarithms to avoid overflows.
%This can mean that one does not obtain the desired complex nth root when
%negative values are complex numbers are the argument as there are multiple
%solutions to the nth root.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(x))
    val=[];
    return
end

xDims=size(x);
if(nargin<2||isempty(dim))
    if(length(xDims)<3)%If x is a vector or matrix.
        if(xDims(1)==1||xDims(2)==1)%If x is a vector
            x=x(:);
            dim=1;
            n=length(x);
        else%x is a matrix, go over the columns.
            dim=2;
            n=xDims(2);
        end
    else%If x is an n-dimensional matrix.
        dim=find(xDims>1,1);
        if(isempty(dim))
            dim=1;
            n=1;
        else
            n=xDims(dim);
        end
    end
else
    n=xDims(dim);
end

%We could do prod(x,dim).^(1/n);, but that is likely to overflow if there
%are many values. Thus, we take the exponent of the sum of the logarithm.
val=exp((1/n)*sum(log(x),dim));

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
