function intVal=monomialIntExpAlt(alpha,UL)
%%MONOMIALINTEXPALT Evaluate the integral of w(x)*x^alpha where
%       w(x)=exp(-x) taken over the range of 0 to UL. This can also be used
%       in a multivariate context where w(x)=exp(-sum(x)) and the integral
%       is over w(x)*prod_{i=1}^N x(i)^alpha(i). That is, this function
%       gives the value of the desired monomial integral taken with an
%       exponential weighting function over the positive halfspace up to
%       UL.
%
%INPUTS: alpha An NX1 or 1XN vector of the exponents of the monomial term.
%              These values od not have to be integers. All elements must
%              be >-1.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The formula was explicitely solved in terms of the gamma function and the
%incomplete gamma function.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(UL))
    UL=Inf; 
end

if(UL<0||~isreal(UL)||isnan(UL))
    error('UL must be real, positive, and not a NaN.')
end

if(UL==0)
    intVal=0;
elseif(isfinite(UL))
    intVal=prod(gamma(alpha+1)-gamma(alpha+1)*gammainc(UL,alpha+1,'upper'));
else%If it is Inf
    intVal=prod(gamma(alpha+1));
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
