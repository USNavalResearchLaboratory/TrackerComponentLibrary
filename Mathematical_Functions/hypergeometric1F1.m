function val=hypergeometric1F1(a,b,z)
%HYPERGEOMETRIC1F1 Evaluate the regular confluent hypergeometric function
%                  1F1(a;b;z) (the Kummer confluent hypergeometric
%                  function).
%
%INPUTS: a,b,z Three real or complex scalar numbers.
%
%OUTPUTS: val The scalar value of 1F1(a;b;z). 
%
%Hypergeometric functions arise in many applications, including
%combinatorics. This function is implemented using Buchholz polynomial as
%described in [1].
%
%REFERENCES:
%[1] J. Abad and J. Sesma, "Computation of the regular confluent
%    hypergeometric function," The Mathematica Journal, vol. 5, no. 4, pp.
%    74-76, 1995.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

nMax=50;

t=z*(a-b/2);

prevTerm=hypergeometric0F1(b,t);
sumVal=prevTerm;

didConverge=false;
for n=1:nMax
    curTerm=BuchholzPolynomial(n,b,z)*hypergeometric0F1(b+n,t)/(2^n*risingFactorial(b,n));
    sumVal=sumVal+curTerm;

    if(eps(sumVal)>abs(curTerm))
        didConverge=true;
        break
    end
end

if(didConverge==false)
   warning('Convergence of hypergeometric1F1 not achieved.') 
end
val=exp(z/2)*sumVal;

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
