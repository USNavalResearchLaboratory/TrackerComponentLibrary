function x=invDiscCosTrans(C,type)
%%DISCCOSTRANS Compute the inverse 1D discrete cosine transformation of a
%              given real array.
%
%INPUTS: C  An (N+1)X1 or 1X(N+1) real vector whose inverse discrete cosine
%           transform is desired.
%     type  A string specifying the type of discrete cosine transform
%           desired. If omitted, the inverse of a unitary discrete cosine
%           transform is performed. Possible values are:
%           'II' Find the inverse of a type 2 discrete cosine transform.
%                The item in element k of the original transform is the sum
%                from n=0 to N-1 of 2*x(n)*cos(pi*k*(n+1/2)/N), where all
%                indexation is from 0.
%           'Unitary' The default if omitted. This is the same as type II,
%                     except for scaling. The scaling of the original
%                     transform whose inverse is desired is
%                     CU(1)=CII(1)/(2*sqrt(N));
%                     CU(2:end)=CII(2:end)*0.5*sqrt(2/N);
%                     where CII is the type 2 transform and CU is the
%                     unitary transform.
%
%OUTPUTS: x The inverse discrete cosine transformation of C.
%
%The algorithm implemented is taken from [1].
%
%This function is the inverse of discCosTrans.
%
%REFERENCES:
%[1] J. Makhoul, "A fast cosine transform in one and two dimensions," IEEE
%    Transactions on Acoustics, Speech, and Signal Processing, vol. ASSP-
%    28, no. 1, pp. 27-34, Feb. 1980.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(C);
if(size(C,1)==N)
    rowVec=true;
else
    rowVec=false;
end

if(nargin<2)
    type='Unitary';
end

switch(type)
    case 'II'
    case 'Unitary'
     %Adjust weighting to make it type II.
     C(1)=C(1)*(2*sqrt(N));
     C(2:end)=C(2:end)/(0.5*sqrt(2/N));
    otherwise
        error('Discrete cosine transform type not suported.')
end

C=[C(:);0];

k=0:1:(N-1);
V=0.5*exp(1j*2*pi*k(:)/(4*N)).*(C(k+1)-1j*C(N-k+1));

%The real command only deals with precision issues. The result should be
%real anyway.
v=real(ifft(V));

n1=0:1:floor((N-1)/2);
n2=floor((N+1)/2):1:(N-1);

%Allocate space.
x=zeros(N,1);
x(2*n1+1)=v(n1+1);
x(2*N-2*n2-1+1)=v(n2+1);

if(rowVec==false)
    x=x';
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
