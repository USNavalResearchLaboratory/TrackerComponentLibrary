function y=upsampleZero(x,P)
%%UPSAMPLEZERO Increase the sampling rate of a signal by inserting zeros
%               in the new samples. This is often used as a first step in
%               interpolatory filters. Adding zeros between samples causes
%               the frequency response of the fft to be aliased. One can
%               then apply a window in the frequency domain and take the
%               ifft to get an interpolated signal. Just adding zeros
%               results in an interpolation with a lot of Gibbs phenomenon.
%
%INPUTS: x An NX1 or 1XN vector or an NXnumVec matrix to upsample. If it is
%          a matrix, then the values are taken per column to upsample. If
%          it is a vector, then the output will have the same orientation
%          as the input (column or row).
%        P The positive integer factor to upsample. P=1 returns the original
%          signal. 
%
%OUTPUTS: y The upsampled signal. If x is a vector, it will be (N*P)X1 or
%           1X(N*P), matcghing the input orientation. If a matrix, it will
%           be (N*P)XnumVec in size.
%
%EXAMPLE:
% x=1:5;
% P=3;%So there will be two zeros between each number.y=
% y=upsampleZero(x,P)
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
numVals=size(x,2);
y=zeros(N*P,numVals);
y(1:P:(N*P),:)=x;

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
