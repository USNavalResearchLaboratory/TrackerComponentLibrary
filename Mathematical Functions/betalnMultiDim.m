function betalnVal=betalnMultiDim(a,b,p)
%%BETALNMULTIDIM Evaluate the natural logarithm of the multivariate beta
%           function. This is defined on page 20 of [1] as the integral
%           beta_p(a,b)=int_{0<A<I_p}
%           det(A)^(a-(1/2)*(p+1))*det(I_p-A)^(b-(1/2)*(p+1)) dA
%           where the integration is over all pXp matrices such that I_p-A
%           is positive definite, where I_p is the pXp identity matrix.
%           Another representation is
%           beta_p(a,b)=gamma_p(a)*gamma_p(b)/gamma_p(a+b)
%           where gamma_p is the p-dimensional multivariate gamma
%           distribution.
%
%INPUTS: a,b Two positive, matrix of values, such that 2*a>p-1 and 2*b>p-1.
%            Both can have the same dimensions, or one is scalar and the
%            other a matrix.
%          p A positive scalar integer>=1.
%
%OUTPUTS: betalnVal The positive, scalar logarithm of multivariate beta
%                   function values. The dimensionality of this is
%                   consistent with that of a and b (whichever is a
%                   matrix).
%
%REFERENCES:
%[1] A. K. Gupta and D. K. NAGAR Matrix Variate Distributions. Boca Raton,
%    Florida: Chapman & Hall/CRC Press, 1999 (Monographs and Surveys in
%    Pure and Applied Mathematics 104)
%
%March 2019 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

betalnVal=gammalnMultiDim(a,p)+gammalnMultiDim(b,p)-gammalnMultiDim(a+b,p);

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
