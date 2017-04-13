function intVal=symIntFirstKindDegen(x,y)
%%SYMINTFIRSTKINDDEGEN Compute the degenerate symmetric integral of the
%                      first kind. This is one half times the integral from
%                      0 to infinite of (t+x)^(-1/2)*(t+y)^(-1) dt. The
%                      parameters can be complex (see below). Sometimes,
%                      the integral is referred to as RC.
%
%INPUTS: x,y  The two parameters of the function. It is assumed that
%             y is not zero and that if the values are complex, the
%             magnitudes of their phases are less than pi.
%
%OUTPUTS: intVal The value of the degenerate symmetric integral of the
%                second kind given x and y.
%
%An algorthm for this function can be found in [1]. and is used for general
%values.
%
%However, if real(x)>0 and (real(y)>=0||imag(y)~=0), then the integral can
%be explicitly solved to get
%intVal=acos(sqrt(x)/sqrt(y))/(sqrt(y)*sqrt(1-x/y));
%Similarly, when x==0 and (imag(y)~=0||real(y)>=0) then
%intVal=pi/(2*sqrt(y));
%Thus, this simpler expressions are used when the necessary criteria are
%fulfilled. Though, care must be taken, because the simplification does not
%work when x=y.
%
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(real(x)>0&&(real(y)>=0||imag(y)~=0))
    intVal=acos(sqrt(x)/sqrt(y))/(sqrt(y)*sqrt(1-x/y));
    
    %If it is not finite, then x is probably too close to y and the other
    %algorithm should be used.
    if(isfinite(intVal))
        return;
    end
elseif(x==0&&(real(y)>=0||imag(y)~=0))
    intVal=pi/(2*sqrt(y));
    return;
end

%Ignoring round-off errors, this should be the relative error of the
%solution.
r=eps;

%If y is negative, a transformation must be made.
if(y<0)
    y=-y;
    coeffVal=sqrt(x/(x+y));
    x=x+y;
else
    coeffVal=1;
end

A0=(x+2*y)/3;
Q=(3*r)^(-1/8)*abs(A0-x);

n=0;
fourProd=1;
while(1)
    lambda=2*sqrt(x)*sqrt(y)+y;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    A=(x+2*y)/3;

    n=n+1;
    fourProd=fourProd*4;%This is 4^n;
    if(Q<fourProd*abs(A)&&n>7)
        break;
    end

    if(n>500)
        warning('Convergence not achieved in 500 iterations.');
        break;
    end
end

s=(y-A0)/(fourProd*A);
p=[(9/8);(159/208);(9/22);(3/8);(1/7);(3/10);0;1];

intVal=coeffVal*polyval(p,s)/sqrt(A);

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
