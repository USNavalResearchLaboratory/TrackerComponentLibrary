function val=meanDev(x,dim)
%%MEANDEV Determine the mean (average) deviation of the elements in x along
%         dimension dim, or along the first non-singleton dimension if dim
%         is omitted. When considering a vector, this is the mean over i of
%         abs(x(i)-mean(x)), whereas the variance is the mean of the square
%         of that term. This is not a robust estimator.
%
%INPUTS: x A vector, matrix or hypermatrix of values.
%      dim The dimension over which the mean deviation should be found. If
%          omitted or an empty matrix is passed, then the first
%          nonsingleton dimension of x will be used.
%
%OUTPUTS: val The mean deviation. This will have the same dimensions as x
%             except dimension dim will be unitary.
%
%EXAMPLE:
% x=magic(4);
% val=meanDev(x)
%One will get val=[4,4,4,4].
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(x))
    val=NaN;
    return;
end

%Select the first non-singleton dimension if dim is not provided.
if(nargin<2||isempty(dim))
    dim=find(size(x)>1,1);
    if(isempty(dim))
        dim=1;
    end
end

val=mean(abs(bsxfun(@minus,x,mean(x,dim))),dim);

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
