function C=discSinCosTrans(x,type)
%%DISCSINCOSTRANS Compute the 1D discrete sine or cosine transformation of
%              a given array. All 16 standard types of sine and cosine
%              transforms are available plus a unitary version of the type
%              II even cosine transform.
%
%INPUTS: x  An (N+1)X1 or 1X(N+1) real or complex vector whose discrete
%           sine or cosine transform is desired.
%     type  A string specifying the type of discrete sine or cosine
%           transform desired. If omitted or an empty matrix is passed, the
%           default is 'CIIe'. Possible values are:
%           'CIe' Perform a type I even discrete cosine transform. The
%                 item in element k of the result is
%                 x(0)+(-1)^k*X(N-1)+2*sum_{j=1}^{N-2}x(j)*cos(pi*j*k/(N-1))
%                 where all indexation is form 0.
%           'CIIe' Perform a type II even discrete cosine transform. The
%                 item in element k of the result is the sum from n=0 to
%                 N-1 of
%                 2*x(n)*cos(pi*k*(n+1/2)/N), where all indexation is from
%                 0.
%         'CIIeU' This is the same as type CIIe, except for scaling. The
%                 scaling is
%                 CU(1)=CII(1)/(2*sqrt(N));
%                 CU(2:end)=CII(2:end)*0.5*sqrt(2/N); where CIIe is the
%                 type 2 transform and CU is the unitary transform.
%         'CIIIe' Perform a type III even discrete cosine transform.
%          'CIVe' Perform a type IV even discrete cosine transform.
%           'CIo' Perform a type I odd discrete cosine transform.
%          'CIIo' Perform a type II odd discrete cosine transform.
%         'CIIIo' Perform a type III odd discrete cosine transform.
%          'CIVo' Perform a type IV odd discrete cosine transform.
%           'SIe' Perform a type I even discrete sine transform.
%          'SIIe' Perform a type II even discrete sine transform.
%         'SIIIe' Perform a type III even discrete sine transform.
%          'SIVe' Perform a type IV even discrete sine transform.
%         'SIIIe' Perform a type III even discrete sine transform.
%          'SIVe' Perform a type IV even discrete sine transform.
%           'SIo' Perform a type I odd discrete sine transform.
%          'SIIo' Perform a type II odd discrete sine transform.
%         'SIIIo' Perform a type III odd discrete sine transform.
%          'SIVo' Perform a type IV odd discrete sine transform.
%
%OUTPUTS: C The (N+1)X1 or 1X(N+1) real discrete cosine or sine
%           transformation of x. 
%
%The algorithms are implemented using the generalized discrete Fourier
%transform (GDFT) based upon the relations of the sequences and the GDFT
%that are given in [1]. The definition of the unitary version of the type
%II transformation is taken from [2]. 
%
%The inverse of this function is invDiscSinCosTrans.
%
%REFERENCES:
%[1] S. Martucci, "Symmetric convolution and the discrete sine and cosine
%    transforms," IEEE Transactions on Signal Processing, vol. 42, no. 5,
%    pp. 1038-1051, May 1994.
%[2] J. Makhoul, "A fast cosine transform in one and two dimensions," IEEE
%    Transactions on Acoustics, Speech, and Signal Processing, vol. ASSP-
%    28, no. 1, pp. 27-34, Feb. 1980.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    type='CIIe';
end

N=length(x);
if(size(x,1)==N)
    colVec=true;
else
    colVec=false;
end

x=x(:);
switch(type)
    case 'CIIeU'
        %Perform the CIIe transform
        C=genDFT([x;flip(x)],0,1/2);
        C=C(1:N);
        
        %Modify the weighting.
        C(1)=C(1)/(2*sqrt(N));
        C(2:end)=C(2:end)*0.5*sqrt(2/N);
    case 'CIe'
        C=fft([x;x((N-1):-1:2)]);
        C=C(1:N);
    case 'CIIe'
        C=genDFT([x;flip(x)],0,1/2);
        C=C(1:N);
    case 'CIIIe'
        C=genDFT([x;0;-flip(x(2:end))],1/2,0);
        C=C(1:N);
    case 'CIVe'
        C=genDFT([x;-flip(x)],1/2,1/2);
        C=C(1:N);
    case 'CIo'
        C=fft([x;flip(x(2:end))]);
        C=C(1:N);
    case 'CIIo'
        C=genDFT([x;flip(x(1:(end-1)))],0,1/2);
        C=C(1:N);
    case 'CIIIo'
        C=genDFT([x;-flip(x(2:end))],1/2,0);
        C=C(1:N);
    case 'CIVo'
        C=genDFT([x;0;-flip(x)],1/2,1/2);
        C=C(1:N);
    case 'SIe'    
        C=1j*fft([0;x;0;-flip(x)]);
        C=C(2:(N+1));
    case 'SIIe'
        C=1j*genDFT([x;-flip(x)],0,1/2);
        C=C(2:(N+1));
    case 'SIIIe'
        C=1j*genDFT([0;x(1:end);flip(x(1:(end-1)))],1/2,0);
        C=C(1:N);
    case 'SIVe'
        C=1j*genDFT([x;flip(x)],1/2,1/2);
        C=C(1:N);
    case 'SIo'
        C=1j*fft([0;x;-flip(x)]);
        C=C(2:(N+1));
    case 'SIIo'
        C=1j*genDFT([x;0;-flip(x)],0,1/2);
        C=C(2:(N+1));
    case 'SIIIo'
        C=1j*genDFT([0;x;flip(x)],1/2,0);
        C=C(1:N);
    case 'SIVo'
        C=1j*genDFT([x;flip(x(1:(end-1)))],1/2,1/2);
        C=C(1:N);
    otherwise
        error('Discrete cosine transform type not suported.')
end

%Get rid of finite precision errors. We know that if the purely input is
%real or imaginary, then the output must be the same.
if(isreal(x))
	C=real(C);
elseif(isImag(x))
    C=1j*imag(C);
end

if(~colVec)
    C=C.';
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
