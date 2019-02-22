function val=HoelderMean(x,p,dim)
%%HOELDERMEAN Compute the Hölder mean (the generalized mean, also known as
%             the power mean) of a set of numbers. The Hölder mean of n
%             numbers equals
%             ((1/n)*sum_{=1}^n x_i)^(1/p)
%             Special cases are p=-Inf (minimum), p=Inf (maximum), p=-1
%             (harmonic mean), p=0 (geometric mean), p=1 (arithmetic mean)
%             and p=2 (root mean squared).
%
%INPUTS: x The set of samples over which the Hölder mean is to be found.
%          How the function handles vectors and matrices is determined by
%          the parameter dim.
%        p The positive or negative exponent of the mean as described
%          above. If this parameter is omitted or an empty matrix is
%          passed, p=1 is used.
%      dim An optional parameter specifying over which dimension of x the
%          Hölder mean is to be taken. If dim is omitted or an empty
%          matrix is passsed, then if x is a vector, the Hölder mean is
%          found over all elements of the vector. If x is a matrix, then
%          the mean is found over the columns of the matrix (resulting in a
%          row vector), and if x is an n-dimensional (n>2) matrix, then the
%          Hölder mean is found over the first non-singleton dimension
%          of the matrix.
%
%OUTPUT: val The Hölder mean of x taken over the appropriate dimension with
%            the exponent p.
%
%The Hölder mean is described in [1]. For the case of n=0, the function
%geometricMean is called.
%
%REFERENCES:
%[1] Cantrell, David W. and Weisstein, Eric W. "Power Mean." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/PowerMean.html
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(p))
   p=1; 
end

xDims=size(x);
if(nargin<3||isempty(dim))
    if(length(xDims)<3)%If x is a vector or matrix.
        if(xDims(1)==1||xDims(2)==1)%If x is a vector
            x=x(:);
            dim=1;
            n=length(x);
        else%x is a matrix, go over the columns.
            dim=2;
            n=xDims(2);
        end
    else%If x is an n-dimensional matrix.
        dim=find(xDims>1,1);
        if(isempty(dim))
            dim=1;
            n=1;
        else
            n=xDims(dim);
        end
    end
else
    n=xDims(dim);
end

switch(p)
    case -Inf
         val=min(x,[],dim);
    case -1
         val=harmonicMean(x,dim);
    case 0
         val=geometricMean(x,dim);
    case 1
         val=mean(x,dim);
    case Inf
         val=max(x,[],dim);
    otherwise
        val=((1/n)*sum(x.^p,dim)).^(1/p);
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
