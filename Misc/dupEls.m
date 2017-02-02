function xDup=dupEls(x,numTimes)
%%DUPELS Return a vector with the elements of the vector x duplicated 
%        numTimes.
%
%INPUTS: x        A vectors whose elements are to be duplicated.
%        numTimes The number of times the elements in x are to be
%                 duplicated.
%
%OUTPUTS: xDup    A vector the same orientation (row or column vector) as x
%                 where the elements of x have been duplicated numTimes.
%                 For example, if x=[1;2;3], then dupEls(x,3) returns
%                 [1;1;1;2;2;2;3;3;3].
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=length(x);
xDup=zeros(xDim*numTimes,1);

%If the input was a row vector.
if(size(x,1)==1)
    xDup=xDup';
end

sel=(0:(xDim-1))*numTimes;
for curDup=1:numTimes
    xDup(sel+curDup)=x;
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
