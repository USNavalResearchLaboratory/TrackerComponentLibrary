function intVal=symIntThirdKind(x,y,z,p)
%%SYMINTTHIRDKIND Compute the symmetric integral of the third kind. This is
%                 3/2 times the integral from 0 to infinity of
%                 ((t+x)*(t+y)*(t+z))^(-1/2)*(t+p)^(-1) dt where p is not
%                 zero and the components can be complex (see below).
%                 Sometimes, the integral is referred to as RJ.
%
%INPUTS: x,y,z,p The parameters of the integral. x,y, and z should have a
%                nonegative real part with at most one of them zero and
%                real(p)>0. Otherwise, let the magnitude of the complex
%                phase of p be less than pi and p~=0, x,y, and z can be
%                real and non-negative with at most one zero, or let two
%                of the variables x,y,z be nonzero and conjugate complex
%                with phases less in magnitude than pi and the third
%                variable be real and nonnegative.
%
%OUTPUTS: intVal The value of the symmetric integral of the third kind
%                given x, y, z, and p.
%
%The algorithm is taken from [1].
%
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Ignoring round-off errors, this should be the relative error of the
%solution.
r=eps;

x0=x;
y0=y;
z0=z;

A0=(x+y+z+2*p)/5;

delta=(p-x)*(p-y)*(p-z);

Q=max([abs(A0-x),abs(A0-y),abs(A0-z),abs(A0-p)])/nthroot(r/4,6);
A=A0;
n=0;

fourProd1=1;
fourProd3=1;
d=(sqrt(p)+sqrt(x))*(sqrt(p)+sqrt(y))*(sqrt(p)+sqrt(z));%d0
e=delta/(fourProd3*d^2);
sumVal=(6/(d*fourProd1))*symIntFirstKindDegen(1,1+e);

while(1)
    lambda=sqrt(x)*sqrt(y)+sqrt(x)*sqrt(z)+sqrt(y)*sqrt(z);
    A=(A+lambda)/4;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    z=(z+lambda)/4;
    p=(p+lambda)/4;

    n=n+1;
    fourProd1=fourProd1*4;%This is 4^n;
    fourProd3=fourProd3*4*4*4;%This is 4^(3*n);

    if(Q<fourProd1*abs(A))
        break;
    end

    if(n>500)
        warning('Convergence not achieved in 500 iterations.');
        break;
    end
    d=(sqrt(p)+sqrt(x))*(sqrt(p)+sqrt(y))*(sqrt(p)+sqrt(z));
    e=delta/(fourProd3*d^2);
    sumVal=sumVal+(6/(d*fourProd1))*symIntFirstKindDegen(1,1+e);
end

X=(A0-x0)/(fourProd1*A);
Y=(A0-y0)/(fourProd1*A);
Z=(A0-z0)/(fourProd1*A);

P=(-X-Y-Z)/2;

E2=X*Y+X*Z+Y*Z-3*P^2;
E3=X*Y*Z+2*E2*P+4*P^3;
E4=(2*X*Y*Z+E2*P+3*P^3)*P;
E5=X*Y*Z*P^2;

intVal=A^(-3/2)*(1-(3/14)*E2+(1/6)*E3+(9/88)*E2^2-(3/22)*E4-(9/52)*E2*E3+(3/26)*E5)/fourProd1+sumVal;

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
