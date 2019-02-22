function vals=GegenbauerPolyChebMoments(lambda,k)
%%GEGENBAUERPOLYCHEBMOMENTS Evaluate the integral
%                int_{-1}^1 T_{2*r}(x)*(1-x^2)^(lambda-1/2) dx
%                for r=0:k, where T_n(x) is the nth order Chebyshev
%                polynomial of the first kind and lambda>-0.5. Gegenbauer
%                polynomials are orthogonal to the weight function
%                (1-x^2)^(lambda-1/2) when integrated over the region
%                (-1,1). This function only returns values with respect
%                to even Chuebyshev polynomials. The values with respect to
%                odd orders are all zero.
%
%INPUTS: lambda The real lambda term in the (1-x^2)^(lambda-1/2). lambda
%               >-1/2.
%             k A value >=0 such that 2*k is the maximum order of the
%               Chebyshev polynomial moment taken.
%
%OUTPUTS: vals A (k+1)X1 vector of moments where vals(i) is the moment with
%              respect to T_{2*(k-1))(x)
%
%The expressions for the integrals are given in [1] and play a role in the
%generation of quadrature points for integrals involving a Gegenbauer
%polynomial weighting function.
%
%REFERENCES:
%[1] D.B. Hunter, H.V. Smith, A quadrature formula of Clenshaw-Curtis type
%    for the Gegenbauer weight-function, Journal of Computational and
%    Applied Mathematics 177 (2005) 389-400.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%G is given in Equation 1.7;
j=(1:k).';
G=cumprod((j-1-lambda)./(j+lambda));

vals=zeros(k+1,1);
%The k=0 value
vals(1)=exp(gammaln(lambda+1/2)-gammaln(lambda+1))*sqrt(pi);

%Equation 1.6
vals(2:(k+1))=vals(1).*G;

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
