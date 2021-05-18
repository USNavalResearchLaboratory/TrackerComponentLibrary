function C=autocorrMats(x,maxLag,scaling)
%%AUTOCORRMATS Compute the sample autocorrelation matrices for a given
%           (assumed zero mean) vector sequence. The sample of a sequence
%           of N vector pairs for a lag of k is
%           C(:,:,k)=x(:,(k+1):N)*x(:,1:(N-k))';
%           is unnormalized and more commonly multiplied by (1/(N-j)) for
%           the unbiased estimate. This is the vector form of computing the
%           autocorrelation of a scalar x using
%           correlation(x,x,[],'unbiased'). However, whereas the
%           "correlation" function returns a symmetric sequence (the
%           negative delay values equal to the positive ones are given),
%           this function only provides delays >=0. In tracking, the
%           compution of such autocorrelation matrices arises in covariance
%           estimation of the system, such as in [1].
%
%INPUTS: x A xDimXN matrix of N xDimX1 vectors.
%   maxLag The output is given for delays of j=0 to maxLag. maxLag can
%          range from 0 to N-1. If this parameter is omitted, the default
%          of N-1 is used.
%  scaling This optionally specifies how the output is normalized. Possible
%          values are:
%          'none'     (The default if omitted or an empty matrix is passed)
%                     No scaling is performed.
%          'biased'   The output C is divided by N.
%          'unbiased' The kth matrix in C is divided by N-(k-1)
%
%OUTPUTS: C The xDimXxDimXmaxLag set of autocorrelation matrices. Matrix
%           C(:,:,k) gives the value for a delay of k-1 (k starting at 1
%           going up to maxLag+1).
%
%The expression for the sample autocorrelation is given in Equation 3.6 of
%[1].
%
%EXAMPLE 1:
%Here, one gets a set of 2D correlation matrices.
% x=[1,  2,   3,   4;
%   -9, -8, -16, -12];
% C=autocorrMats(x)
%
%EXAMPLE 2:
%Here, we show that for a complex scalar sequence, the result, after
%reshaping and mirroring, is the same as the output of the function
%"correlation" when considering a 1D sequence (differences being due to
%finite precision limitiations).
% x=rand(9,1)+1j*randn(9,1);
% c=correlation(x,x,[],'unbiased');
% c1=autocorrMats(x(:).');
% c1=c1(:);
% %One will see that the following are equal, except for finite precision
% %differences.
% c(:)
% [conj(c1(end:-1:2));c1(:)]
%
%REFERENCES:
%[1] J. Duník, O. Straka, O. Kost, and J. Havlík, "Noise covariance
%    matrices in state-space models: A survey and comparison of estimation
%    methods-part I," International Journal of Adaptive Control and Signal
%    Processing, vol. 31, no. 11, pp. 1505-1543, Nov. 2017.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(x,1);

if(nargin<3||isempty(scaling))
   scaling='unbiased'; 
end

N=size(x,2);

if(nargin<2||isempty(maxLag))
    maxLag=N-1;
end

%The indexation is different from Equation 3.6 of [1], but the result is
%the same.
C=zeros(xDim,xDim,maxLag+1);
for k=0:maxLag
    C(:,:,k+1)=x(:,(k+1):N)*x(:,1:(N-k))';
end

switch(scaling)
    case 'none'
    case 'biased'
        C=C/N;
    case 'unbiased'
        for k=0:maxLag
           C(:,:,k+1)=C(:,:,k+1)/(N-k);
        end
    otherwise
        error('Unknown scaling specified.')  
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
