function x=invDiscSinCosTrans(C,type)
%%INVDISCSINCOSTRANS Compute the inverse of a 1D discrete sine or cosine
%              transformation of a given array. All 16 standard types of
%              sine and cosine transforms are available plus a unitary
%              version of the type II even cosine transform.
%
%INPUTS: C  An (N+1)X1 or 1X(N+1) real or complex vector whose inverse
%           discrete cosine or sine transform is desired.
%     type  A string specifying the type of discrete inverse sine or cosine
%           transform desired. If omitted or an empty matrix is passed, the
%           default is 'CIIe'. Possible values are:
%           'CIe' Invert a type I even discrete cosine transform. The
%                 item in element k of the result is
%                 x(0)+(-1)^k*X(N-1)+2*sum_{j=1}^{N-2}x(j)*cos(pi*j*k/(N-1))
%                 where all indexation is form 0.
%           'CIIe' Invert a type II even discrete cosine transform. The
%                 item in element k of the result is the sum from n=0 to
%                 N-1 of
%                 2*x(n)*cos(pi*k*(n+1/2)/N), where all indexation is from
%                 0.
%         'CIIeU' This is the same as type CIIe, except for scaling. The
%                 scaling is
%                 CU(1)=CII(1)/(2*sqrt(N));
%                 CU(2:end)=CII(2:end)*0.5*sqrt(2/N); where CIIe is the
%                 type 2 transform and CU is the unitary transform.
%         'CIIIe' Invert a type III even discrete cosine transform.
%          'CIVe' Invert a type IV even discrete cosine transform.
%           'CIo' Invert a type I odd discrete cosine transform.
%          'CIIo' Invert a type II odd discrete cosine transform.
%         'CIIIo' Invert a type III odd discrete cosine transform.
%          'CIVo' Invert a type IV odd discrete cosine transform.
%           'SIe' Invert a type I even discrete sine transform.
%          'SIIe' Invert a type II even discrete sine transform.
%         'SIIIe' Invert a type III even discrete sine transform.
%          'SIVe' Invert a type IV even discrete sine transform.
%         'SIIIe' Invert a type III even discrete sine transform.
%          'SIVe' Invert a type IV even discrete sine transform.
%           'SIo' Invert a type I odd discrete sine transform.
%          'SIIo' Invert a type II odd discrete sine transform.
%         'SIIIo' Invert a type III odd discrete sine transform.
%          'SIVo' Invert a type IV odd discrete sine transform.
%
%OUTPUTS: x The inverse discrete cosine transformation of C.
%
%The algorithms are implemented using the relations between the forward and
%inverse transformation gien in Appendix A of [1]. The definition of the
%unitary version of the type II transformation  is taken from [2]. 
%
%This function is the inverse of discSinCosTrans.
%
%REFERENCES:
%[1] J. Makhoul, "A fast cosine transform in one and two dimensions," IEEE
%    Transactions on Acoustics, Speech, and Signal Processing, vol. ASSP-
%    28, no. 1, pp. 27-34, Feb. 1980.
%[2] S. Martucci, "Symmetric convolution and the discrete sine and cosine
%    transforms," IEEE Transactions on Signal Processing, vol. 42, no. 5,
%    pp. 1038-1051, May 1994.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(C);
if(nargin<2)
    type='CIIe';
end

switch(type)
    case 'CIe'
        x=discSinCosTrans(C,'CIe')/(2*(N-1));
    case 'CIIe'
        x=discSinCosTrans(C,'CIIIe')/(2*N);
    case 'CIIeU'
        %Adjust weighting to make it type CIIe.
        C(1)=C(1)*(2*sqrt(N));
        C(2:end)=C(2:end)/(0.5*sqrt(2/N));
        
        %Do the inverse for 'CIIe'
        x=discSinCosTrans(C,'CIIIe')/(2*N);
    case 'CIIIe'
        x=discSinCosTrans(C,'CIIe')/(2*N);
    case 'CIVe'
        x=discSinCosTrans(C,'CIVe')/(2*N);
    case 'CIo'
        x=discSinCosTrans(C,'CIo')/(2*N-1);
    case 'CIIo'
        x=discSinCosTrans(C,'CIIIo')/(2*N-1);
    case 'CIIIo'
        x=discSinCosTrans(C,'CIIo')/(2*N-1);
    case 'CIVo'
        x=discSinCosTrans(C,'CIVo')/(2*N+1);
    case 'SIe'
        x=discSinCosTrans(C,'SIe')/(2*(N+1));
    case 'SIIe'
        x=discSinCosTrans(C,'SIIIe')/(2*N);
    case 'SIIIe'
        x=discSinCosTrans(C,'SIIe')/(2*N);
    case 'SIVe'
        x=discSinCosTrans(C,'SIVe')/(2*N);
    case 'SIo'
        x=discSinCosTrans(C,'SIo')/(2*N+1);
    case 'SIIo'
        x=discSinCosTrans(C,'SIIIo')/(2*N+1);
    case 'SIIIo'
        x=discSinCosTrans(C,'SIIo')/(2*N+1);
    case 'SIVo'  
        x=discSinCosTrans(C,'SIVo')/(2*N-1);
    otherwise
        error('Discrete cosine transform type not suported.')
end

%Get rid of finite precision errors. We knw that if the purely input is
%real or imaginary, then the output must be the same.
if(isreal(C))
	x=real(x);
elseif(isImag(C))
    x=1j*imag(x);
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
