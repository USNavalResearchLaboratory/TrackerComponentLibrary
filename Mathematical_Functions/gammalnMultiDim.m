function gammalnVal=gammalnMultiDim(a,p)
%%GAMMALNMULTIDIM Evaluate the natural logarithm of the multivariate gamma
%           function for a positive real argument. The multivariate gamma
%           function is defined as an integral starting on page 18 in [1]
%           as
%           Gamma_p(a)=int_{A>0} det(A)^{a-(1/2)*(p+1)*exp(trace(-A)) dA
%           where the integral is over all symmetric positive definite pXp
%           matrices (to which A>0 refers). It is assumed that 2*a>p-1.
%           This can also be expressed as
%           Gamma_p(a)=pi^((1/4)*p*(p-1))*prod_{i=1}^p gamma(a-(1/2)*(i-1))
%           The logarithm just changes that product to a sum and gamma to
%           gammaln.
%
%INPUTS: a A vector or matrix of positive values such that 2*a>p-1.
%        p A positive scalar integer>=1.
%
%OUTPUTS: gammalnVal The positive, scalar logarithm of multivariate gamma
%                    function values. The dimensionality of this is the
%                    same as that of a.
%
%This function can produce finite results in problems where gammaMultiDim
%would overflow.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. NAGAR Matrix Variate Distributions. Boca Raton,
%    Florida: Chapman & Hall/CRC Press, 1999 (Monographs and Surveys in
%    Pure and Applied Mathematics 104)
%
%March 2019 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(p<0||fix(p)~=p)
    error('p must be an integer >0.') 
end

if(~isreal(a)||~isreal(p))
    error('Both a and p must be real.')
end

if(2*a<=p-1)
    error('The multivariate gamma function is only defined for 2*a>p-1'); 
end

logSqrtPi=log(sqrt(pi));
gammalnVal=gammaln(a);
piCoeff=0;
for k=2:p
    piCoeff=piCoeff+logSqrtPi;
    gammalnVal=piCoeff+gammalnVal+gammaln(a+(1-k)/2);
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
