function z=real2ComplexSignal(x,M)
%%REAL2COMPLEXSIGNAL Given the real in-phase part of a signal, this
%          function returns a complex signal with in-phase and quadrature
%          parts. This function can also simultanoeusly decimate the signal
%          by 1/2 or interpolate the signal out by an integer multiple. The
%          signal length must be divisible by 2 to perform decimation.
%
%INPUTS: x The NX1  or 1XN signal. N must be divisble by 2 if M=1/2.
%        M A positive integer >=1. The length of the output sequence will
%          be M*N. However, one can also pass M=1/2 to decimate the signal
%          by 1/2. If this parameter is omitted or an empty matrix is
%          passed, a default value of 1 is used.
%
%OUTPUTS: z The (M*N)X1 complex signal.
%
%The algorithm is the "analytic" signal computation algorithm of [1]. The
%paper focusses on the case where N is even. However, the generalization
%for interpolation for N being odd is straightforward. This function can be
%used to take a real samples signal and return a "quadrature" component for
%signal processing purposes. 
%
%EXAMPLE:
% t=(1:100).';
% x=exp(1j*2*pi*t/100);
% real2ComplexSignal(real(x))-x
%One will get a result that is numerically zero for all of the elements.
%Note that using x=exp(-1j*2*pi*t/100), will, of course cause this function
%to return the same output, as the real part of the signal can't express
%the 180 degree phase shift. Thus, there is an ambiguity in the result.
%Note that it also works with an odd-length signal:
% t=(1:101).';
% x=exp(1j*2*pi*t/101);
% real2ComplexSignal(real(x))-x
%One also gets numerically zero results.
%
%REFERENCES:
%[1] S. L. Marple Jr., "Computing the discrete-time 'analytic' signal via
%    FFT," IEEE Transactions on Signal Processing, vol. 47, no. 9, pp.
%    2600-2603, Sep. 1999.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(M))
    M=1;
end

N=length(x);
X=fft(x(:),N);

if(M==1/2)%For this, N is a multiple of 2.
    z=(1/2)*ifft([X(1)+X(N/2+1);2*X(2:(N/2))],N/2);
else
    if(mod(N,2)==0)
        z=ifft([X(1);2*X(2:(N/2));X(N/2+1);zeros(M*N-N/2-1,1)],M*N);
    else
        z=ifft([X(1);2*X(2:((N+1)/2));X((N+1)/2+1);zeros(M*N-(N+1)/2-1,1)],M*N);
    end
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
