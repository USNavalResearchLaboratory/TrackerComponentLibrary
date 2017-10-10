function [fracVal,E]=fracFourier(f,a,order)
%%FRACFOURIER Compute the discrete fractional Fourier transform.
%
%INPUTS: f The NXnumEls matrix of real or complex sequences over which the
%          fractional Fourier transform is desired. The transform is
%          performed over each column. If N is too large, overflow errors
%          will occur in the computation of the fractional Fourier
%          transform matrix.
%        a The order of the fractional Fourier transform.
%    order This parameter is only used if algorithm 0 is selected. This is
%          the order used when computing the Fourier transform matrix. It
%          can be from 0 to N-1. The default if omitted or an empty matrix
%          is passed is N/2. Alternatively the E output of this function
%          can be passed into order to speed up repeated calls to algorithm
%          0 for sequences x having the same length.
%
%OUTPUTS: fracVal The NXnumEls values of the fractional Fourier transform
%                 over each column of f. 
%               E This the set of eigenvectors from the function
%                 FourierTransMatEig. Passing this for multiple
%                 calls with f having the same length will speed up the
%                 algorithm.
%
%The fractional Fourier transform arises during the efficient computation
%of ambiguity functions and when performing a Wigner transform. This
%function implements the algorithm for the discrete fractional Fourier
%transform described in [1]. It is not the fastest algorithm, but it is
%generally more accurate than the fast approximations.
%
%We will note that for a column vector x,
%fft(x)==ifftshift(fracFourier(fftshift(x),1))*sqrt(length(x))
%ifft(x)==ifftshift(fracFourier(fftshift(x),3))/sqrt(length(x))
%and
%fracFourier(x,1)=fftshift(fft(ifftshift(x)))/sqrt(length(x))
%fracFourier(y,3)=fftshift(ifft(ifftshift(y)))*sqrt(length(y))
%to within finite precision limits.
%Additionally, f==fracFourier(fracFourier(f,a,order),-a,order) within
%finite precision limits. That is, the transform with -a is the inverse.
%
%EXAMPLE:
%Here is a simple example:
% x=0.0:0.01:(2*pi); x=x(:);
% N=length(x);
% f=sin(x);
% a=1.5;
% f0=fracFourier(f,a,N-1);
% figure(1)
% clf
% hold on
% plot(x,abs(f0),'-r','linewidth',4);
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] A. Bultheel and H. E. Martínez Sulbaran, "Computation of the
%    fractional Fourier transform," Applied and Computational Harmonic
%    Analysis, vol. 16, no. 3, pp. 182-202, May 2004.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(f))
   fracVal=[];
   E=[];
   return;
end

N=size(f,1);
numEls=size(f,2);

if(nargin<3||isempty(order))
    order=N/2;
end

fracVal=zeros(N,numEls);

[fracVal(:,1),E]=discFracFourier(f(:,1),a,order);

for curEl=2:numEls
    fracVal(:,curEl)=discFracFourier(f(:,curEl),a,E);
end
end

function [fracVal,E,lambda]=discFracFourier(f,a,order)
%%DISCFRACFOURIER This function implements the discrete fractional Fourier
%                 transformation described in Section 6 of [1].
%
%REFERENCES:
%[0] A. Bultheel and H. E. Martínez Sulbaran, "Computation of the
%    fractional Fourier transform," Applied and Computational Harmonic
%    Analysis, vol. 16, no. 3, pp. 182-202, May 2004.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    N=length(f);
    
    if(~isscalar(order))
        E=order;
    else 
        E=FourierTransMatEig(N,order);
    end
    
    isEven=(mod(N,2)==0);
    lambda=exp(-1j*pi/2*a*([0:(N-2),N-1+isEven])).';
    
    idx=mod((0:N-1)+fix(N/2),N)+1;
    f=f(:);
    fracVal(idx)=E*(lambda.*(E'*f(idx)));
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
