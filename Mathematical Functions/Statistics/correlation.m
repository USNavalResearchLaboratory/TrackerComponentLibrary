function corr=correlation(x,y)
%%CORRELATION Estimate the correlation between two discrete sequences x and
%             y. The strict definition of the autocorrelation is 
%             r(m)=sum_{n=-Inf}^Inf x(n+m)*y(n)' m=0,+/-1, +/-2, etc.
%             which is equivalent to 
%             r(m)=sum_{n=-Inf}^Inf x(n)*y(n+m)'
%             However, since the input sequences are finite, the
%             correlation returned by this function if x and y are length
%             N is r(m)=sum_{n=0}^(N-m-1} x(n+m)y(n)' for m>=0 and
%             r(-m)' for m<0 The values are arranged into a vector corr
%             such that corr(m)=r(m-N) for m=1,2,...(2*N-1)
%             Thus, r(0) is corr(N), r(-1) is corr(N-1) and r(1) is
%             corr(N+1).
%
%INPUTS: x A numDimX X 1 or a 1X numDimX sequence.
%        y If y is omitted, then an autocorrelation of x is computed (the
%          same as y=x. If y is provided, then is is a numDimY X1  or
%          1 X numDimY sequence. If numDImX~=numDimY, the shorted sequence
%          is zero-padded at the end to make tham have equal length.
%
%OUTPUTS: corr The (2*N-1) X 1 autocorrelation vector such that if
%              N=max(numDimX,numDimY), corr(N+k) is the correlation of an
%              offset of k in the sequence.
%
%As discussed in Chapter 2.9.3 of [1], discrete correlations can be easily
%expressed in terms of convolutions. As mentioned in Chapter 5.7 of [1],
%the convolution can be expressed in terms of discrete fourier transforms.
%Thus, that is how this is implemented as it is faster than explicit
%computation of all of the terms for long sums.
%
%EXAMPLE:
%Consider example 2.46 of [1].
% x=[1;3;-2;1;2;-1;4;4;2];
% y=[2;-1;4;1;-2;3];
% c=correlation(x,y);
% N=max([length(x),length(y)]);
% figure()
% stem((1-N):(N-1),c)
% %The stem command plots the correlation for the proper delays.
% %Now, we consider the autocorrelation of x
% c=correlation(x,x);
% N=length(x);
% figure()
% stem((1-N):(N-1),c)
%
%Note that in economics, the definition of cross correlation is typically
%applied to the sequences x-mean(x) and y-mean(y) instead of directly to x
%and y.
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach,
%    3rd ed. Boston: McGraw Hill, 2006.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    y=x;%Compute an autocorrelation.
end

numX=length(x);
numY=length(y);
if(numX<numY)
    x=[x(:);zeros(numY-numX,1)];
elseif(numY<numX)
    y=[y(:);zeros(numX-numY,1)];
end
numEls=max([numX,numY]);

%FFTs are used for efficiency with long sequences. This is similar to how
%FFTs are used with basic convolutions.
FFTLength=2^nextpow2(2*numEls-1);
r=ifft(fft(x(:),FFTLength).*conj(fft(y(:),FFTLength)));

%Keeps values corresponding to lags -(numEls-1):(numEls-1)
corr=[r((end-numEls+2):end);r(1:numEls)];

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
