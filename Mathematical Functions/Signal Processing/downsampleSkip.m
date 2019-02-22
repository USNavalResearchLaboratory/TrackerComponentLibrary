function y=downsampleSkip(x,P)
%%DOWNSAMPLESKIP Reduce the number of samples in x by skipping values.
%                This only takes the first and every Pth value.
%
%INPUTS: x An NX1 or 1XN vector or an NXnumVec matrix to downsample. If it
%          is a matrix, then the values are taken per column to downsample.
%          If it is a vector, then the output will have the same
%          orientation as the input (column or row).
%        P The positive integer factor to downsample. P=1 returns the
%          original signal. 
%
%OUTPUTS: y The downsamples signal. If x is a vector, this will be a vector
%           of the same orientation.
%
%EXAMPLE:
% x=1:10;
% P=3;
% y=downsampleSkip(x,P)
%One will get y=[1,4,7,10];
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

isTranspose=false;
if(isvector(x))
    if(size(x,2)>size(x,1))
        isTranspose=true;
        x=x(:);
    end
end

N=size(x,1);
sel=1:P:N;
y=x(sel,:);

if(isTranspose)
    y=y.';
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
