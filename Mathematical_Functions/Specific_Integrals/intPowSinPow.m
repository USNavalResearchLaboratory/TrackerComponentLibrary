function val=intPowSinPow(u,n,m)
%%INTPOWSINPOW Evaluate the integral of u^n*sin(u)^m du.  A definite
%           integral can be evaluated, or an indefinite integral (with a
%           particular additive constant).
%
%INPUTS: u A 2XN (for definite integral) or a 1XN (for indefinite
%          integrals) set of N points. For definite integrals, u(1,:) are
%          the real lower bounds and u(2,:) are the real upper bounds. For
%          indefinite integrals, the integral is evaluated at the points in
%          u. The values in u should be real.
%        n The positive integer exponent of u.
%        m The positive integer exponent of the sine term.
%
%OUTPUTS: val The 1XN set of values of the integral of u^n*sin(u)^m.
%
%This function implements formulas 4 and 5 of Section 2.631 of [1].
%
%REFERENCES:
%[1] I. S. Gradshteyn and I. M. Ryzhik, Tables of Integrals, Series, and
%    Products, Corrected and Enlarged Edition. New York: Academic Press,
%    1980, translated from Russian.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(u,1);

if(isempty(u))
   val=[];
   return;
end

if(numDim==1)%An indefinite integral
    val=indefIntPowSinPow(u,n,m);
else%A definite integral
    val=indefIntPowSinPow(u(2,:),n,m)-indefIntPowSinPow(u(1,:),n,m);
end

end

function val=indefIntPowSinPow(u,n,m)

if(mod(m,2)==0)%Use formula 4.
    m=m/2;
    
    val=0;
    for k=0:(m-1)
        a=2*(m-k);
        x=a*u;
        val=val+(-1)^k*binomial(2*m,k)*(1/a)^(n+1)*intPowCos(x,n);
    end

    val=(-1)^m/(2^(2*m-1))*val;
    val=val+binomial(2*m,m)*u.^(n+1)/(2^(2*m)*(n+1));
else%Use formula 5.
    m=(m-1)/2;
    
    val=0;
    for k=0:m
        a=2*(m-k)+1;
        x=a*u;
        val=val+(-1)^k*binomial(2*m+1,k)*(1/a)^(n+1)*intPowSin(x,n);
    end
    
    val=val*(-1)^m/(2^(2*m));
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
