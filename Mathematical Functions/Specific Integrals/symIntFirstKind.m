function intVal=symIntFirstKind(x,y,z)
%%SYMINTFIRSTKIND Evaluate the symmetric integral of the first kind. This
%                 is 0.5 times the integral from 0 to infinity of
%                 ((t+x)*(t+y)*(t+z))^(-1/2) dt, where the square root is
%                 taken real and positive if x,y,and z are positive. x, y,
%                 and z can be complex (see below). This integral arises
%                 when numerically evaluating elliptic integrals.
%                 Sometimes, the integral is referred to as RF.
%
%INPUTS: x,y,z The three parameters of the function. It is assumed that at
%              most one is zero. When complex, the nonzero ones are
%              assumed to have complex phase angles less than pi in
%              magnitude.
%
%OUTPUTS: intVal The value of the symmetric integral of the first kind
%                given x, y, and z.
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

%If z=0, then a simplified recursion can be used.
if(z==0)
    x=sqrt(x);
    y=sqrt(y);
    
    n=0;
    while(1)
        xPrev=x;
        yPrev=y;
        
        x=(xPrev+yPrev)/2;
        y=sqrt(xPrev*yPrev);
        
        if(abs(x-y)<2.7*sqrt(r)*abs(x))
            break;
        end
        
        n=n+1;
        if(n>500)
            warning('Convergence not achieved in 500 iterations.');
            break;
        end
    end
    
    intVal=pi/(x+y);
elseif(y==z)%The degenerate case
    intVal=symIntFirstKindDegen(x,y);
else
    A0=(x+y+z)/3;
    Q=max([abs(A0-x),abs(A0-y),abs(A0-z)])/nthroot(3*r,6);
    n=0;
    A=A0;
    fourProd=1;
    while(1)
        %One should not combine the terms in the square root products into
        %a single square root or else the function will give the wrong
        %answer when using complex variables.
        lambda=sqrt(x)*sqrt(y)+sqrt(x)*sqrt(z)+sqrt(y)*sqrt(z);
        A=(A+lambda)/4;
        x=(x+lambda)/4;
        y=(y+lambda)/4;
        z=(z+lambda)/4;

        %The termination condition.
        n=n+1;
        fourProd=fourProd*4;%This is 4^n
        
        if(Q<abs(A)*fourProd)
           break;
        end
        if(n>500)
            warning('Convergence not achieved in 500 iterations.');
            break;
        end
    end

    X=1-x/A;
    Y=1-y/A;
    Z=-X-Y;
    E2=X*Y-Z^2;
    E3=X*Y*Z;

    intVal=(1-(1/10)*E2+(1/14)*E3+(1/24)*E2^2-(3/44)*E2*E3)/sqrt(A);
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
