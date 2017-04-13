function wz=Faddeeva(z)
%%FADDEEVA Compute the Faddeeva function (also known as the Kramp function)
%          of the possibly complex argument z. This is a scaled
%          complementary error function. It equals exp(-z^2)*(1-erf(z))
%          However, Matlab's built-in erf function cannot handle complex
%          numbers. This function implements the algorithm described in
%          [1].
%
%INPUTS: z A scalar, vector or matrix of values at which one wishes to
%          evaluate the Faddeeva function.
%
%OUTPUTS: wz A matrix having the same dimensions as z in which the values
%            of the Faddeeva function are computed for the values in z.
%
%The Faddeeva arises in a number of physics problems involving waves. It is
%also related to the complex error function, Fresnel integrals. The
%algorithm of [1] is used. It is a modification of the algorithm of [2].

%REFERENCES:
%[1] G. P. M. Poppe and C. M. J. Wijers, "More efficient computation of the
%    complex error function," ACM Transactions on Mathematical Software,
%    vol. 16, no. 1, pp. 38-46, Mar. 1990.
%[2] W. Gautschi, "Efficient computation of the complex error function,"
%    SIAM Journal on Numerical Analysis, vol. 7, no. 1, pp. 187-198, Mar.
%    1970.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVals=numel(z);
wz=zeros(size(z));

zList=z;

for curZ=1:numVals
    z=zList(curZ);

    xOrig=real(z);
    yOrig=imag(z);

    %Put the values into the first quadrant.
    x=abs(xOrig);
    y=abs(yOrig);
    z=x+1i*y;

    %Evaluates exp(-z^2)*erfc(-1i*z);

    %The choice of v, x0,and y0 as in Section 2.2
    x0=6.3;
    y0=4.4;

    %The rho term in Equation 2.7
    rho=sqrt((x/x0)^2+(y/y0)^2);

    if(rho>0.292)
        %Use the Laplace continued fraction.
        if(rho>1)%Outer Region Q.
            %Equation 2.13
            nu=round(3+1442/(26*rho+77));%The minimum number of terms needed.
            N=0;
            h=0;
        else%%Inner Region R
            %Equation 2.11
            s=(1-y/y0)*sqrt(1-rho^2);

            %Equation 2.14;
            nu=round(16+26*s);
            N=round(7+34*s);
            h=1.88*s;
        end

        r=0;
        s=0;

        %Equation 2.5
        for n=nu:-1:0
            r=(1/2)/(h-1i*z+(n+1)*r);
            if(n<=N&&h~=0)
               s=r*((2*h)^n+s);
            end
        end

        if(h>0)
            w=(2/sqrt(pi))*s;
        else
            w=(2/sqrt(pi))*r;
        end
    else
        %Inner Region S
        %This is evaluated using Equation 2.16

        %Equation 2.18
        sPrime=(1-0.85*y/y0)*rho;

        %Equation 2.17
        N=round(6+72*sPrime);

        %Do the sum in Equation 2.16
        prodVal=z*2*1i/sqrt(pi);
        %The first term for n=0
        sumVal=prodVal;
        for n=1:N
            prodVal=prodVal*z^2*(2*n-1)/(n*(2*n+1));
            sumVal=sumVal+prodVal;
        end

        w=exp(-z^2)*(1+sumVal);
    end

    %Use the transformation of Section 3 to account for the proper quadrant of
    %the estimate.
    if(yOrig<0)
        w=2*exp(-z^2)-w;
        if(xOrig>0)
            w=conj(w);
        end
    elseif(xOrig<0)
        w=conj(w); 
    end
    
    wz(curZ)=w;
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
