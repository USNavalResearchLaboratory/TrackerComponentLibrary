function y=genDFT(x,k0,n0)
%%GENDFT Evaluate the generalized discrete Fourier transform. This is 
%        y(k)=sum_{n=0}^{N-1}x(n)*exp(-1j*(2*pi/N)*(k+k0)*(n+n0));
%        where the indexation on k and n starts at zero. This is the same
%        as the normal Fourier transform (the fft function), except for the
%        offsets a and b.
%
%INPUTS: x An NX1 or 1XN vector or an NXnumVec matrix of numVec vectors.
%          The values can be real or complex.
%       k0 A real offset in the Fourier domain. If this and the next input
%          are omitted, then k0=n0=-(N-1)/2 is used, which is the centered
%          DFT.
%       n0 A real offset in the time domain.
%
%OUTPUTS: y The generalized discrete Fourier transform of the vectors in x.
%           This has the same dimensions as x.
%
%The generalized DFT is introduced in [1]. Rather than using the fast
%algorithm described in [1] for implementation, the transformation is
%implemented by transforming the fft algorithms. As can be found in common
%texbooks, shifts in time and frequency can be achieved by multiplying the
%values in time and frequency by complex exponentials. Thus, this
%pre-multiplied x by a vextor of complex exponentials and then
%post-multiplied the output of the fft by a vector of complex exponentials.
%An explicit expression for this operation is given in Equation D.2 of [2]
%
%The generalized discrete Fourier transform can play a role in the
%computation of discrete sine and cosine transforms as can be seen in [3].
%
%REFERENCES:
%[1] G. Bongiovanni, P. Corsini, and G. Frosini, "One-dimensional and two-
%    dimensional generalized discrete Fourier transforms," IEEE
%    Transactions on Acoustics, Speech, and Signal Processing, vol. 24,
%    no. 1, pp. 97-99, Feb. 1976.
%[2] K. R. Rao, D. N. Kim, and J. J. Hwang, Fast Fourier Transform:
%    Algorithms and Applications. Springer, 2010.
%[3] S. Martucci, "Symmetric convolution and the discrete sine and cosine
%    transforms," IEEE Transactions on Signal Processing, vol. 42, no. 5,
%    pp. 1038-1051, May 1994.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isvector(x)&&size(x,1)<size(x,2))
    colVec=false;
    x=x.';
else
    colVec=true;
end

N=size(x,1);

if(nargin<2)%The default is the centered DFT
   k0=-(N-1)/2; 
   n0=k0;
end

k=(0:(N-1))';
%The weights in time.
wt=exp(-1j*2*pi*k0*k/N);
%The weights in frequency.
wf=exp(-1j*2*pi*(k+k0)*n0/N);

y=wf.*fft(wt.*x);

if(~colVec)
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
