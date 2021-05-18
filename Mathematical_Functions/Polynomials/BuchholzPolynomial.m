function p=BuchholzPolynomial(n,b,z)
%%BUCHHOLZPOLYNOMIAL Compute the value of a Buchhols polynomial of a given
%             order with the specified parameters. The Buchholz polynomial
%             p_n(b,z) arises when dealing with hypergeometric functions as
%             in [1].
%
%INPUTS: n The order of the polynomial. This is a non-negative integer.
%        b A number. This can be complex.
%        z A number. This can be complex.
%
%OUTPUTS: p The value of the Buchhols polynomial with the specified
%           parameters.
%
%The polynomials are implemented using Equations 5, 6, and 7 of [1].
%
%REFERENCES:
%[1] J. Abadand J. Sesma, "Buchholz polynomials: A family of polynomials
%    relating solutions of confluent hypergeometric and Bessel functions,"
%    Journal of Computational and Applied Mathematics, vol. 101, no. 1-2,
%    pp. 237-241, 15 Jan. 1999.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

maxFOrder=fix(n/2);
f=zeros(maxFOrder+1,1);

%Equation 6
f(0+1)=1;
for s=1:maxFOrder
    sumVal=0;
    for r=0:(s-1)
        sumVal=sumVal+binomial(2*s-1,2*r)*4^(s-r)*(abs(BernoulliNum(2*(s-r)))/(s-r))*f(r+1);
    end
    
    f(s+1)=-(b/2-1)*sumVal;
end

%Equation 7
g=zeros(n+1,1);
g(0+1)=1;
for m=1:n
    sumVal=0;
    for k=0:fix((m-1)/2)
        sumVal=sumVal+binomial(m-1,2*k)*4^(k+1)*abs(BernoulliNum(2*(k+1)))/(k+1)*g(m-2*k-1+1);
    end
    g(m+1)=-(1i*z/4)*sumVal;
end

%Equation 5
sumVal=0;
for s=0:maxFOrder
    sumVal=sumVal+binomial(n,2*s)*f(s+1)*g(n-2*s+1);
end

p=(1i*z)^n/factorial(n)*sumVal;

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
