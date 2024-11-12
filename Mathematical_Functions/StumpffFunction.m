function vals=StumpffFunction(xIn,k)
%%STUMPFFFUNCTION Evaluate the Stumpff function, supporting orders k=0 to
%            k=21. This is just a scaled version of Battin's universal
%            function from Chapter 4.5 of [1]. The function evaluates the
%            series C_k(x)=sum_{n=0}^Inf (-1)^n*x^n/factorial(k+2*n)
%            The Stumpff function is described in [2] and arises in
%            Keplerian celestial mechanics.
%
%INPUTS: x A real scalar or matrix
%        k The scalar degree of the Stumpff function. This function
%          supports integers degrees from 0 to 21.
%
%OUTPUTS: vals A matrix the same size as x holding the Stumpff function of
%              degree k evalauted at x.
%
%Explicit expressions for the infinite series were obtained using
%Mathematica up to order 21. For a magnitude of x>200, those expressions
%are used. For a magnitude of x<200, the series sum is evaluated, since the
%explicit expressions are less accurate for small values of x and the sum
%is less accurate for large values of x -due to finite precision
%limitiations.
%
%REFERENCES:
%[1] R. H. Battin, An introduction to the Mathematics and Methods of
%    Astrodynamics. Reston, Virginia: American Institute of Aeronautics
%    and Astronautics, 1999.
%[2] E. W. Weisstein, "Stumpff Function", Eris Weisstein's World of
%    Physics, Online: https://scienceworld.wolfram.com/physics/StumpffFunction.html
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

vals=zeros(size(xIn));
N=numel(xIn);

for curEl=1:N
    x=xIn(curEl);

    if(k==fix(k)&&k>=0&&k<=21)
        if(abs(x)<=200)
            val=0;
            for n=0:100
                term=(-1)^n*(x^n/factorial(k+2*n));
                %This loop should always terminate long before n=100.
                if(abs(term)<eps(val))
                    break;
                end
                val=val+term;
            end
            vals(curEl)=val;
            continue
        end
    
        switch(k)
            case 0
                val=cos(sqrt(x));
            case 1
                val=sin(sqrt(x))/sqrt(x);
            case 2
                val=(1-cos(sqrt(x)))/x;
            case 3
                val=1/x-sin(sqrt(x))/x^(3/2);
            case 4
                val=(-2+x+2*cos(sqrt(x)))/(2*x^2);
            case 5
                val=(-6+x)/(6*x^2)+sin(sqrt(x))/x^(5/2);
            case 6
                val=(24+(-12+x)*x-24*cos(sqrt(x)))/(24*x^3);
            case 7
                val=1/x^3+(-20+x)/(120*x^2)-sin(sqrt(x))/x^(7/2);
            case 8
               val=(-720+x*(360+(-30+x)*x)+720*cos(sqrt(x)))/(720*x^4);
            case 9
                val=1/x^4+(840+(-42+x)*x)/(5040*x^3)+sin(sqrt(x))/x^(9/2);
            case 10
                val=1/x^5+(-20160+x*(1680+(-56+x)*x))/(40320*x^4)-cos(sqrt(x))/x^5;
            case 11
                val=1/x^5+(-60480+x*(3024+(-72+x)*x))/(362880*x^4)-sin(sqrt(x))/x^(11/2);
            case 12
                val=-(1/x^6)+1/(2*x^5)+(-151200+x*(5040+(-90+x)*x))/(3628800*x^4)+cos(sqrt(x))/x^6;
            case 13
                val=-(1/x^6)+1/(6*x^5)+(-332640+x*(7920+(-110+x)*x))/(39916800*x^4)+sin(sqrt(x))/x^(13/2);
            case 14
                val=1/x^7-1/(2*x^6)+1/(24*x^5)+(-665280+x*(11880+(-132+x)*x))/(479001600*x^4)-cos(sqrt(x))/x^7;
            case 15
                val=1/x^7-1/(6*x^6)+1/(120*x^5)+(-1235520+x*(17160+(-156+x)*x))/(6227020800*x^4)-sin(sqrt(x))/x^(15/2);
            case 16
                val=-(1/x^8)+1/(2*x^7)-1/(24*x^6)+1/(720*x^5)-1/(40320*x^4)+1/(3628800*x^3)-1/(479001600*x^2)+1/(87178291200*x)+cos(sqrt(x))/x^8;
            case 17
                val=-(1/x^8)+1/(6*x^7)-1/(120*x^6)+1/(5040*x^5)+(-3603600+x*(32760+(-210+x)*x))/(1307674368000*x^4)+sin(sqrt(x))/x^(17/2);
            case 18
                val=1/x^9-1/(2*x^8)+1/(24*x^7)-1/(720*x^6)+1/(40320*x^5)-1/(3628800*x^4)+1/(479001600*x^3)-1/(87178291200*x^2)+1/(20922789888000*x)-cos(sqrt(x))/x^9;
            case 19
                val=1/x^9-1/(6*x^8)+1/(120*x^7)-1/(5040*x^6)+1/(362880*x^5)-1/(39916800*x^4)+1/(6227020800*x^3)-1/(1307674368000*x^2)+1/(355687428096000*x)-sin(sqrt(x))/x^(19/2);
           case 20
                val=-(1/x^10)+1/(2*x^9)-1/(24*x^8)+1/(720*x^7)-1/(40320*x^6)+1/(3628800*x^5)-1/(479001600*x^4)+1/(87178291200*x^3)-1/(20922789888000*x^2)+1/(6402373705728000*x)+cos(sqrt(x))/x^10;
            otherwise%k==21
                val=-(1/x^10)+1/(6*x^9)-1/(120*x^8)+1/(5040*x^7)-1/(362880*x^6)+1/(39916800*x^5)-1/(6227020800*x^4)+1/(1307674368000*x^3)-1/(355687428096000*x^2)+1/(121645100408832000*x)+sin(sqrt(x))/x^(21/2);
        end
        vals(curEl)=val;
    else
        error('This function only supports integer values of k from 0 to 21.')
    end
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
