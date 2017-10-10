function y=discHartleyTrans(x,N,dim)
%%DISCHARTLEYTRANS Compute the discrete Hartley transform. For a length N
%                  sequence, this is
%                  sum_{n=0}^{N-1}x(n+1)*(cos(2*pi/N*n*k)+sin(2*pi/N*n*k))
%                  This is similar to the discrete Fourier transform,
%                  except there is no imaginary value times the sine term
%                  and it is plus the sine term instead of minus. Note that
%                  x=(1/N)*discHartleyTrans(discHartleyTrans(x)). That is,
%                  it is its own inverse within a scale factor.
%
%INPUTS: x A vector or matrix over a particular dimensions of which (e.g.
%          columns), the discrete Hartely transform should be taken.
%        N The length of the transform to use. The default if this
%          parameter is omitted or an empty matrix is passed is to make the
%          transform the same length as the dimension of x over which the
%          transform is taken.
%      dim An optional parameter specifying over which dimension of x the
%          transform should be taken. If omitted or an empty matrix is
%          passed, the transform is taken over the entire vector, if a
%          vector, or the first non-unitary dimension if a matrix is
%          passed.
%
%OUTPUTS: y The discrete Hartley transform of the vectors in x. The
%           dimensions will be the same as x, though dimension dim will be
%           length N. 
%
%The Discrete Hartelty transform is defined in [1]. Here, we use the
%normalizing constant in a manner consistent with how Matlab uses if for
%the fft. Though it is possible to implement the discrete Hartely transform
%to be faster than an fft, we implement it here using two real ffts as it
%is much simpler.
%
%REFERENCES:
%[1] R. N. Bracewell, "Discrete Hartley transform," Journal of the Optical
%    Society of America, vol. 73, no. 12, pp. 1832-1835, Dec. 1983.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(dim))
    dim=find(size(x)>1,1);
    if(isempty(dim))
        dim=1;
    end
end

if(nargin<2||isempty(N))
    N=size(x,dim);
end

Xr=fft(real(x),N,dim);
Xi=fft(imag(x),N,dim);

y=(real(Xr)-imag(Xr))+1i*(real(Xi)+imag(Xi));

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
