function C=discSinCosMat(N,type)
%%DISCSINCOSMAT Obtain a matrix C such that when multiplied by a vector x,
%               C*x equals the value of a discrete sine or cosine transform
%               of x. All 16 discrete transforms (even and odd) are
%               implemented.
%
%INPUTS:   N The positive integer length of the transformation.
%       type A strong indicating the type of discrete cosine transformation
%            to perform. Possible values are:
%           'CIe'  (The default if omitted or an empty matrix is passed)
%                  The even type I cosine transform.
%           'CIIe' The even type II cosine transform.
%           'CIIIe' The even type III cosine transform.
%           'CIVe' The even type IV cosine transform.
%           'CIo' The odd type I cosine transform.
%           'CIIo' The odd type II cosine transform.
%           'CIIIo' The odd type III cosine transform.
%           'CIVo' The odd type IV cosine transform.
%           'SIe' The even type I sine transform.
%           'SIIe' The even type II sine transform.
%           'SIIIe' The even type III sine transform.
%           'SIVe' The even type IV sine transform.
%           'SIo' The odd type I sine transform.
%           'SIIo' The odd type II sine transform.
%           'SIIIo' The odd type III sine transform.
%           'SIVo' The odd type IV sine transform.
%
%OUTPUTS: C The NXN matrix to perform the selected discrete sine or cosine
%           transformation.
%
%Expressions for all of the transform matrices are given in the Appendix of
%[1]. There, the relations between the matrices for determining the inverse
%transformation matrices are given.
%
%REFERENCES:
%[1] S. Martucci, "Symmetric convolution and the discrete sine and cosine
%    transforms," IEEE Transactions on Signal Processing, vol. 42, no. 5,
%    pp. 1038-1051, May 1994.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(type))
   type='CIe'; 
end

switch(type)
    case 'CIe'
        N=N-1;
        
        m=0:N;
        C=2*cos(pi*bsxfun(@times,m(:),m)/N);
        C(:,1)=C(:,1)/2;
        C(:,end)=C(:,end)/2;
    case 'CIIe'
        m=0:(N-1);
        C=2*cos(pi*bsxfun(@times,m(:),m+1/2)/N);
    case 'CIIIe'
        m=0:(N-1);
        C=2*cos(pi*bsxfun(@times,m(:)+1/2,m)/N);
        C(:,1)=C(:,1)/2;
    case 'CIVe'
        m=0:(N-1);
        C=2*cos(pi*bsxfun(@times,m(:)+1/2,m+1/2)/N);
    case 'CIo'        
        m=0:(N-1);
        C=2*cos(2*pi*bsxfun(@times,m(:),m)/(2*N-1));
        C(:,1)=C(:,1)/2;
    case 'CIIo'
        m=0:(N-1);
        C=2*cos(2*pi*bsxfun(@times,m(:),m+1/2)/(2*N-1));
        C(:,end)=C(:,end)/2;
    case 'CIIIo'
        m=0:(N-1);
        C=2*cos(2*pi*bsxfun(@times,m(:)+1/2,m)/(2*N-1));
        C(:,1)=C(:,1)/2;
    case 'CIVo'
        N=N+1;
        m=0:(N-2);
        C=2*cos(2*pi*bsxfun(@times,m(:)+1/2,m+1/2)/(2*N-1));
    case 'SIe'
        N=N+1;
        m=1:(N-1);
        C=2*sin(pi*bsxfun(@times,m(:),m)/N);
    case 'SIIe'
        n=0:(N-1);
        C=2*sin(pi*bsxfun(@times,n(:)+1,n+1/2)/N);
    case 'SIIIe'
        m=0:(N-1);
        C=2*sin(pi*bsxfun(@times,m(:)+1/2,m+1)/N);
        C(:,end)=C(:,end)/2;
    case 'SIVe'
        m=0:(N-1);
        C=2*sin(pi*bsxfun(@times,m(:)+1/2,m+1/2)/N);
    case 'SIo'
        N=N+1;
        m=1:(N-1);
        C=2*sin(2*pi*bsxfun(@times,m(:),m)/(2*N-1));
    case 'SIIo'
        N=N+1;
        n=0:(N-2);
        C=2*sin(2*pi*bsxfun(@times,n(:)+1,n+1/2)/(2*N-1));
    case 'SIIIo'
        N=N+1;
        m=0:(N-2);
        C=2*sin(2*pi*bsxfun(@times,m(:)+1/2,m+1)/(2*N-1));
    case 'SIVo'
        m=0:(N-1);
        C=2*sin(2*pi*bsxfun(@times,m(:)+1/2,m+1/2)/(2*N-1));
        C(:,end)=C(:,end)/2;
    otherwise
        error('Unknown transform selected')
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
