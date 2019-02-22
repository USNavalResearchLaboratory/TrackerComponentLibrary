function intVal=symIntThirdKindDegen(x,y,z)
%%SYNINTTHIRDKINDDEGEN Compute the degenerate symmetric integral of the
%                      third kind. This is 3/2 times the integral from 0 to
%                      infinity of ((t+x)*(t+y))^(-1/2)*(t+z)^(-3/2) dt.
%                      This integral arises when numerically evaluating
%                      elliptic integrals. Sometimes, the integral is
%                      referred to as RD.
%
%INPUTS: x,y,z The three parameters of the function. z cannot be zero and
%              at most one of x and y can be zero. When complex, the values
%              are assumed to have complex phase angles less than pi
%              in magnitude.
%
%OUTPUTS: intVal The value of the degenerate symmetric integral of the
%                third kind given x, y, and z.
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

if(x==0&&y==z)
    intVal=(3*pi/4)*y^(-3/2);
else
    A0=(x+y+3*z)/5;
    Q=max([abs(A0-x),abs(A0-y),abs(A0-z)])/nthroot(r/4,6);
    
    sumTerm=0;
    n=0;
    A=A0;
    fourProd=1;
    while(1)
        lambda=sqrt(x)*sqrt(y)+sqrt(x)*sqrt(z)+sqrt(y)*sqrt(z);

        sumTerm=sumTerm+1/(fourProd*sqrt(z)*(z+lambda));

        x=(x+lambda)/4;
        y=(y+lambda)/4;
        z=(z+lambda)/4;
        A=(x+y+3*z)/5;

        n=n+1;

        fourProd=fourProd*4;%This is 4^n
        if(Q<fourProd*abs(A)&&n>7)
           break;
        end
        
        if(n>500)
            warning('Convergence not achieved in 500 iterations.');
            break;
        end
    end

    X=(A0-x)/(fourProd*A);
    Y=(A0-y)/(fourProd*A);
    Z=-(X+Y)/3;

    E2=X*Y-6*Z^2;
    E3=(3*X*Y-8*Z^2)*Z;
    E4=3*(X*Y-Z^2)*Z^2;
    E5=X*Y*Z^3;

    intVal=(1/fourProd)*A^(-3/2)*(1-(3/14)*E2+(1/6)*E3+(9/88)*E2^2-(3/22)*E4-(9/52)*E2*E3+(3/26)*E5)+3*sumTerm;
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
