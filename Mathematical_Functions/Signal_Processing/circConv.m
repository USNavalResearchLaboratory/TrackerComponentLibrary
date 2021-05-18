function convVal=circConv(x,y,period)
%%CIRCONV Take the circular convolution of sequences x and y. If x is the
%         same length as y and that length is l, the nth element of the
%         output is
%          convVal(n)=sum_{m=0}^l x(mod(m,period)+1)*y(mod(n-m,period)+1)
%         For sequences of different lengths and period not equaling the
%         length of any sequence, x and y are first aliased to be length
%         period (possibly with zero padding) and then the circular
%         convolution of the given period is performed. If
%         period=length(x)+length(y)-1, then the result is the same as the
%         conv function.
%
%INPUTS: x An N1XnumSeq set of numSeq real or complex vectors. Note that if
%          all of the x sequences are the same, then this can  just be a
%          single N1X1 vectior
%        y An N2XnumSeq set of numSeq real or complex vectors.
%   period The period of the circular convolution. If this parameter is
%          omitted or an empty matrix is passed, then period=min(N1,N2);
%
%OUTPUTS: convVal A periodXnumSeq vector that is the circular convolution
%                 of each vector in x with the corresponding one in y.
%
%Circular convolutions are discussed in Chapter 5.4.2 of [1]. The
%implementation here is with ffts. Such an implementation is mentioned in
%the chapter, but is left as a homework problem.
%
%EXAMPLE:
% x=[1;2;3;4];
% y=[12;-10;0;1];
% z=circConv(x,y)
%One will find that z=[-26;17;20;19].
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach,
%    3rd ed. Boston: McGraw Hill, 2006.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(period))
    period=min(length(x),length(y)); 
end

xAlias=aliasSequence(x,period);
yAlias=aliasSequence(y,period);

convVal=ifft(bsxfun(@times,fft(xAlias,period,1),fft(yAlias,period,1)),period,1);
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
