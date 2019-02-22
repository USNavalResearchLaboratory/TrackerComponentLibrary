function symMean=symmetricMean(x,orders)
%%SYMMETRICMEAN Evaluate the system of polynomial functions of n numbers
%          known as a symmetric mean or the mean power product. For two
%          orders given, this is
%          (1/(n*(n-1))*sum_{i~=j}x(i)^orders(1)*x(j)^orders(2) 
%          For three orders, this continues to
%          (1/(n*(n-1)*(n-2))*sum_{i~=j~=k}x(i)^orders(1)*x(j)^orders(2)*x(k)^orders(3) 
%          and the pattern continues for higher orders. For just one order
%          given, it is the just mean(x.^orders(1)). Symmetric means arise
%          in the evlauation of k-statistics. That is, the computation of
%          unbiased central moments.
%
%INPUTS: x A length n vector of real points.
%   orders A length r vector of the orders of the exponents, n<=r.
%      
%OUTPUTS: symMean The value of the symmetric mean.
%
%This implements the formula given in [1]. For r=1, the mean value is
%explicitely found. For r>1, the problem is equivalent to computing a
%matrix permanent, so the perm algorithm is used rather than evaluating all
%binomial(n,r)*factorial(r) terms in the sum.
%
%EXAMPLE:
% x=[0, 8, -21, -13, 26, -19];
% orders=[1,3];
% symMean=symmetricMean(x,orders)
%One will get symMean=-27002.8.
%
%REFERENCES:
%[1] J. W. Tukey, "Keeping moment-like sampling computations simple," The
%    Annals of Mathematical Statistics, vol. 27, no. 1, pp. 37-54, 1956.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);
r=length(orders);

if(n<=r)
    error('It is required that n>r.') 
end

x=x(:);
orders=orders(:);

if(r==1)
    symMean=mean(x.^orders);
else
    X=bsxfun(@power,x,orders.');
    symMean=perm(X)/prod(n:-1:(n-r+1));
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
