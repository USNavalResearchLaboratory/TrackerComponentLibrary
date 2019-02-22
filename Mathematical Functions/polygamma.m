function val=polygamma(k,z)
%%POLYGAMMA Evaluate the polygamma function. This function is the same as
%           the psi function except it will work with negative values of z.
%           The polygamma function is the value of the kth derivative of
%           the derivative of gamma(z) divided by gamma(z). That is, it is
%           the kth derivative of the digamma function.
%
%INPUTS: This function can be used with 1 or two inputs. If two inputs are
%        given, then it is polygamma(k,z). If only one is given, then it is
%        polygamma(z) where k is taken to be zero. The inputs are:
%        k The number of derivatives to take, k>=0. The digamma function is
%          k=0.
%        z A matrix of points at which the polygamma function is to be
%          evaluated.
%
%OUTPUTS: val The values of the polygamma function. This is the same size
%             as z.
%
%For values of z>0, this function just calls psi(k,z). For values of z<0,
%this function uses the identity of [1] to deal with negative values. For
%example, for k=0, one uses the identity:
%psi(1-x)=psi(x)+pi/tan(pi*x)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Polygamma Function." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/PolygammaFunction.html
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==1)
    z=k;
    k=0;
end

numVals=numel(z);
val=zeros(size(z));

for curVal=1:numVals
    x=z(curVal);
    if(x>0)
        val(curVal)=psi(k,x);
    else
        x=1-x;
        val(curVal)=(-1)^k*(psi(k,x)+pi^(k+1)*cotDerivVal(pi*x,k));
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
