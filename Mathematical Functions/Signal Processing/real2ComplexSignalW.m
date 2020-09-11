function z=real2ComplexSignalW(s,returnPadded)
%%REAL2COMPLEXSIGNALW Given the real in-phase part of a signal, this
%          function returns a complex signal with in-phase and quadrature
%          parts. Unlike the function real2ComplexSignal, this function
%          uses the algorithm of [1] to produces an analytic signal that
%          has properties that help avoid alising issues when using the
%          Wigner-Ville transform.
%
%INPUTS: s The NX1 or 1XN real signal.
% returnPadded If true the analytic signal will be padded by N zeros.
%          Otherwise, just the length N analytic signal is returned. The
%          zero-padding is typically desired when computing the Wigner-Ville
%          transform. the default if omitted or an empty matrix is passed
%          is false.
%
%OUTPUTS: z The NX1 (or (2N)X1 is padded) complex analytic signal.The real
%           components of the first N elements will match (within finite-
%           precision bounds) the values of s.
%
%The algorithm is given in Section III of [1].
%
%EXAMPLE:
%Here, we demonstrate that the real components of the analytic signal match
%those of the original signal.
% x=randn(14,1);
% max(abs((real(real2ComplexSignalW(x))-x)./x))
%The returned value will be as close to 0 as one might expect within finite
%precision limits.
%
%REFERENCES:
%[1] J. M. O'Toole, M. Mesbah, and B. Boashash, "A new discrete analytic
%    signal for reducing aliasing in the discrete Wigner-Ville
%    distribution," IEEE Transactions on Signal Processing, vol. 56, no.
%    11, pp. 5427-5434, Nov. 2008.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(returnPadded))
    returnPadded=false; 
end

N=length(s);

%Step 1; zero pad to double the length.
sa=[s(:);zeros(N,1)];

%Step 2, take the DFT.
Sa=fft(sa);

%Steps 2 and 4: Filter and take the IDFT.
z=ifft([Sa(0+1);2*Sa(1+(1:(N-1)));Sa(N+1);zeros(N-1,1)]);

if(returnPadded)
    z((N+1):(2*N))=0; 
else
    z=z(1:N);
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
