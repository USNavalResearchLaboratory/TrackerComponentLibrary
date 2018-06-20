function val=GegenbauerPolyChebMoment(lambda,k)
%%GEGENBAUERPOLYCHEBMOMENT Evaluate the integral
%                int_{-1}^1 T_{k}(x)*(1-x^2)^(lambda-1/2) dx, where T_n(x)
%                is the nth order Chebyshev polynomial of the first kind
%                and lambda>=-0.5. The expression (1-x^2)^(lambda-1/2) is
%                known as a Gegenbauer polynomial. 
%
%INPUTS: lambda A real scalar value for which the moment is desired.
%               lambda >-1/2.
%             k The scalar integer order of the Chebyshev polynomial moment
%               taken; k>=0.
%
%OUTPUTS: val The value of the moment. Note that all odd moments are zero.
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

if(mod(k,2))
    val=zeros(size(lambda));
    return;
end

k=k/2;
val=exp(gammaln(lambda+1/2)-gammaln(lambda+1))*sqrt(pi);

if(k>0)
    %G is given in Equation 1.7;
    j=(1:k).';
    G=prod((j-1-lambda)./(j+lambda));

    %Equation 1.6
    val=val.*G;
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
