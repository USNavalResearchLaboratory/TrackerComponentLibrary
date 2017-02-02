function C=discCosTrans(x,type)
%%DISCCOSTRANS Compute the 1D discrete cosine transformation of a given
%              real array. The weighting can be done in one of two ways.
%
%INPUTS: x  An (N+1)X1 or 1X(N+1) real vector whose discrete cosine
%           transform is desired.
%     type  A string specifying the type of discrete cosine transform
%           desired. If omitted, a unitary discrete cosine transform is
%           performed. Possible values are:
%           'II' Perform a type 2 discrete cosine transform. The item in
%                element k of the result is the sum from n=0 to N-1 of
%                2*x(n)*cos(pi*k*(n+1/2)/N), where all indexation is from
%                0.
%      'Unitary' The default if omitted. This is the same as type II,
%                except for scaling. The scaling is
%                CU(1)=CII(1)/(2*sqrt(N));
%                CU(2:end)=CII(2:end)*0.5*sqrt(2/N); where CII is the type
%                2 transform and CU is the unitary transform.
%
%OUTPUTS: C The (N+1)X1 or 1X(N+1) real discrete cosine transformation of
%           x. 
%
%The algorithm implemented is taken from [1].
%
%The inverse of this function is invDiscCosTrans.
%
%REFERENCES:
%[1] J. Makhoul, "A fast cosine transform in one and two dimensions," IEEE
%    Transactions on Acoustics, Speech, and Signal Processing, vol. ASSP-
%    28, no. 1, pp. 27-34, Feb. 1980.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    type='Unitary';
end

N=length(x);
if(size(x,1)==N)
    rowVec=true;
else
    rowVec=false;
end

x=x(:);
%Allocate space for the output.
v=zeros(N,1);
n1=0:1:floor((N-1)/2);
n2=floor((N+1)/2):1:(N-1);
v(n1+1)=x(2*n1+1);
v(n2+1)=x(2*N-2*n2-1+1);

k=(0:1:(N-1));

Wk4N=exp(-1j*2*pi*k/(4*N));
C=2*real(Wk4N(:).*fft(v));

switch(type)
    case 'II'
    case 'Unitary'
        %Modify the weighting on the type 2.
        C(1)=C(1)/(2*sqrt(N));
        C(2:end)=C(2:end)*0.5*sqrt(2/N);
    otherwise
        error('Discrete cosine transform type not suported.')
end

if(~rowVec)
    C=C';
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
