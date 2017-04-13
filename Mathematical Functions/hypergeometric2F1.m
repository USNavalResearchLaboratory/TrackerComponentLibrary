function val=hypergeometric2F1(a,b,c,z)
%%HYPERGEOMETRIC2F1 Evaluate the Gaussian hypergeometric function
%            2F2(a,b,c,z) for real arguments with z in the region -1<=z<=1.
%            The Gaussian hypergeometric function arises in a number of
%            combinatorial problems.
%
%INPUTS: a,b,c Real scalar values.
%            z A real value such that z<=1.
%
%OUTPUTS: val The value of the Gaussian hypergeometric function
%             2F1(a,b,c,z).
%
%Implementing a numerically stable version of the Gaussian hypergeometric
%function that converges in a reasonable amount of time is difficult. Here,
%the algorithm of [1] for z<1 is implemented along with the asymptotic
%solution for z=1. For -1<z<0.5, a Taylor series solution is used, with a
%transformation needed for negative values of z. For z>0.5, convergence of
%the series would be slow. Thus, a transformation is used. However, values
%of a, b, and c lead to numerical instabilities if c-a-b is almost zero.
%Thus, two special case fixes from [1] are implemented to deal with those
%scenarios. Additionally, a different transformation from [1] is used for
%z<-1. That transformation has its own numerical instabilities if a-b is
%almost an integer. Thus, a special fix is also implemented to deal with
%that case.
%
%Unlike in [1], the gamma function built into Matlab is used rather than
%the Chebyshev approximation listed. However, the forward differencing
%terms given in the paper cannot be directly evaluated due to finite
%precision limitations. These are terms like
%(1/gamma(z+epsVal)-1/gamma(z))/epsVal . If epsVal is large, things are
%fine, but when small, problems arise. Thus, special techniques described
%in the paper utilizing Chebyshev polynomial approximations are used. The
%paper does not provide the coefficients for the polynomials. They were,
%however, obtained directly from the author of [1].
%
%EXAMPLES:
%As noted in [2], exact solutions are available for a few points. For
%example:
% v1=hypergeometric2F1(1/3,2/3,5/6,27/32)
% v2=hypergeometric2F1(1/4,1/2,3/4,80/81)
% v3=hypergeometric2F1(1/8,3/8,1/2,2400/2401)
% v4=hypergeometric2F1(1/6,1/2,2/3,125/128)
% v5=hypergeometric2F1(1/12,5/12,1/2,1323/1331)
% v6=hypergeometric2F1(1/12,5/12,1/2,121/125)
% %Should respectively equal
% v1True=8/5
% v2True=9/5
% v3True=2/3*sqrt(7)
% v4True=(4/3)*2^(1/6)
% v5True=(3/4)*(11)^(1/4)
% v6True=2^(1/6)*15^(1/4)/(4*sqrt(pi))*(gamma(1/3)^3/gamma(1/4)^2)*(1+sqrt(3))
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.
%[2] Weisstein, Eric W. "Hypergeometric Function." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/HypergeometricFunction.html
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %If it does not converge in this many terms, then it probably will not
    %converge.
    maxNumTerms=1000;

    %Special cases
    if(z==1)
        val=gamma(c)*gamma(c-a-b)/(gamma(c-a)*gamma(c-b));
        return;
    end

    %For the implemented cases in Table I, we transform z to a new variable
    %that lies between 0 and 1/2. This also requires transforming the other
    %parameters. For case I, we do not implement the special case. rather,
    %we just have an error.
    if(z<-1)
        a1=a;
        b1=c-b;
        c1=a-b+1;
        a2=b;
        b2=c-a;
        c2=b-a+1;
        w=1/(1-z);

        k=round(a-b);
        epsVal=a-b-k;
        
        if(abs(epsVal)<=0.05)%The test tolerance
            if(a<b)
                temp=a; 
                a=b;
                b=temp;
            end
            k=round(a-b);
            epsVal=a-b-k;
            
            val=caseI(a,b,c,w,k,epsVal,maxNumTerms);
        else%a-b is not close to being an integer. We can use Equation 2
            %directly.
            f1=hypGeom2F1Taylor(a1,b1,c1,w,maxNumTerms);
            f2=hypGeom2F1Taylor(a2,b2,c2,w,maxNumTerms);

            %Equation 2
            val=w^a*(gamma(c)*gamma(b-a)/(gamma(b)*gamma(c-a)))*f1+...
                w^b*(gamma(c)*gamma(a-b)/(gamma(a)*gamma(c-b)))*f2;
        end
    elseif((-1<=z) && (z<0))%Case II in Table I
        %We are using the transformation from Equation 3.
        a1=a;
        b1=c-b;
        c1=c;
        
        w=z/(z-1);
        
        %-1<=z<0, there are no special cases, so the sum in Equation
        %1 can be used directly via the hypGeom2F1Taylor function and
        %the result transformed back.
        f=hypGeom2F1Taylor(a1,b1,c1,w,maxNumTerms);
        val=((1-w)^a)*f;
    elseif((0<=z)&&(z<=0.5))%Case III in Table I
        %No transformation is needed.
        %0<=z<=0.5, so the sum in Equation 1 can be used directly.
        val=hypGeom2F1Taylor(a,b,c,z,maxNumTerms);
    elseif((0.5<z)&&(z<=1))%Case IV in Table I
        %We are using the transformation from Equation 4. We have two sets
        %of coefficients, one for each hypergeometric function. However, if
        %the solution is too close to a region with singularities, we will
        %use the integral solution.
        a1=a;
        b1=b;
        c1=a+b-c+1;
        
        a2=c-a;
        b2=c-b;
        c2=c-a-b+1;
        
        w=1-z;
        
        %If c-a-b is approximately equal to an integer, then a fix is
        %needed, because each of the two terms in Equation 4 is infinite,
        %but their sum remains finite.
        
        k=round(c-a-b);
        epsVal=c-a-b-k;
        %If c-a-b is too close to being an integer.
        if(abs(epsVal)<0.05)
            if(k>=0)
                val=caseIVa(a,b,c,w,k,epsVal,maxNumTerms);
            else
                k=-k;
                val=caseIVb(a,b,c,w,k,epsVal,maxNumTerms);
            end
        else%c-a-b is not close to being an integer. We can use
            %Equation 4 directly.
            val1=hypGeom2F1Taylor(a1,b1,c1,w,maxNumTerms);
            val2=hypGeom2F1Taylor(a2,b2,c2,w,maxNumTerms);

            %Equation 4
            val=(gamma(c)*gamma(c-a-b)/(gamma(c-a)*gamma(c-b)))*val1+w^(c-a-b)*(gamma(c)*gamma(a+b-c)/(gamma(a)*gamma(b)))*val2;
        end
    else
        error('This function requires that abs(z)<=1 or z<-1 with abs(a-b)>0.05.')
    end
end

function val=hypGeom2F1Taylor(a,b,c,z,maxNumTerms)
%%HYPGEOM2F1TAYLOR This is the Gausian hypergeometric function implemented
%                 for the case where 0<=z<=0.5, implemented using its power
%                 series. Care is taken to avoid overflows when computing
%                 the ratios of the falling factorial terms. Enough terms
%                 are used to reach machine precision accuracy.
%
%Section III of [1] describes a criterion used for terminating the sum.
%Here, the criterion used is that the sum terminates when the new terms
%would not change the current value as they are less than the mantissa of
%the current floating point value.
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.

%Evaluate the hypergeometric function using the sum of Equation 1 computed
%using the recursion of Equation 26.
term=1;
val=1;
didSucceed=false;
for k=1:maxNumTerms
    term=term*(a+k-1)*(b+k-1)*z/((c+k-1)*k);
    
    %If the term to add is so small it would not change the result.
    if(abs(term)<eps(val))
        didSucceed=true;
        break;
    end
    
    val=val+term;
end

if(didSucceed==false)
     warning('The algorithm failed to converge.')
end

end

function val=caseI(a,b,c,w,k,epsVal,maxNumTerms)
%%CASEI This implements special case I of [1], which occurs when a-b is
%%close to being an integer.
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.

%First, calculate the finite sum in Equation 16.
termVal=w^b*gamma(k+epsVal);%The n=0 term.
sumVal=termVal;
for n=1:(k-1)
    termVal=-termVal*((-1+b+n)*(-1-a+c+n)*w)/((epsVal+k-n)*n);
    sumVal=sumVal+termVal;
end
val=gamma(c)/(gamma(a)*gamma(c-b))*sumVal;

gammaA=gamma(a);
gammaCA=gamma(c-a);
%Next, compute the infinite sum in Equation 16

%g values that do not depend on n. g6 and g7 are given by Equation 15.
%However, to reduce finite precision effects, the complicated method of
%Section III for computing such quantities is used via the
%invGammaForwardDiff function.
g6=-invGammaForwardDiff(a-k,-epsVal);
g7=invGammaForwardDiff(c-a+k,epsVal);

%f values for the n=0 case.
%These definitions of f are from Appendix A.
f1e=risingFactorial(a-k-epsVal,k)/gammaA;
f2=1/gammaCA;
f2e=1/gammaCA;

%Equation B3
f3=-gamma(1);
f3e=-gamma(1-epsVal);
%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B1
    f3=f3/(k+n);
    f3e=f3e/(k+n+epsVal);
end

%Equation B4
f4=gamma(1);
f4e=gamma(1+epsVal);

f5=w^(a);
f5e=w^(a-epsVal);

%From appendix B
g1=0;
%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B5
    g1=(a-epsVal+n-1)*g1-risingFactorial(a-k,n+k-1)/gammaA;
end

%From appendix B
g2=0;

%As described in Section III of the text, finite differences of gamma
%functions can be difficult to compute to a high precision. Thus, instead
%of using 
%g3=-(gamma(1-epsVal)-gamma(1))/epsVal;
%from Appendix B or the limit g3=gamma(1)*polygamma(1); as epsVal->0,
%a solution utilizing Chebyshev polynomials is used. Note that said
%solution must be incremented up to time n=0 using Equation B16. For g4, a
%similar problem exists and again. Problems with using
%g4=(gamma(1+epsVal)-1)/epsVal exist with small values of g4, with the
%limit as epsVal->0 being g4=gamma(1)*polygamma(1). Thus Chebyshev
%interpolation is again used.
[g3,g4]=computeg3g4(epsVal);

%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B7
    g3=g3/(k+n+epsVal)+1/((k+n+epsVal)*factorial(k+n));
end

%A special approximation for the forward difference is needed to compute
%g5.
g5=-powForwardDiff(a,w,-epsVal);

sumVal=0;
n=0;
didSucceed=false;
outerCoeff=(-1)^k*gamma(c);
while(1)
    aPoch=risingFactorial(a,n);
    caPoch=risingFactorial(c-a,k+n);
    
    term2Add=aPoch*f2e*f3e*f4*f5*g6...
             -f1e*caPoch*f3*f4e*f5e*g7...
             -g1*f2e*f3e*f4*f5...
             +f1e*g2*f3e*f4*f5...
             +f1e*f2*g3*f4*f5...
             -f1e*f2*f3*g4*f5...
             -f1e*f2*f3*f4e*g5;
    term2Add=outerCoeff*term2Add;
    
    if(eps(sumVal)>abs(term2Add))
        didSucceed=true;
        break;
    end
    sumVal=sumVal+term2Add;
    
	n=n+1;
    
    if(n>maxNumTerms||~isfinite(sumVal))
        break;
    end
    
    %These definitions of f are from Appendix A.
    %f1e=risingFactorial(a-k-epsVal,n+k)/gammaA;
    f1e=f1e*(a+n-1-epsVal);
    
    %f2=risingFactorial(c-a+k,n)/gammaCA;
    f2=f2*(c-a+k+n-1);
    
    %f2e=risingFactorial(c-a+k+epsVal,n)/gammaCA;
    f2e=f2e*(c-a+epsVal+k+n-1);
    
    %Equation B1
    f3=f3/(k+n);
    f3e=f3e/(k+n+epsVal);
    
    %Equation B2
    f4=f4/(n);
    f4e=f4e/(n-epsVal);

    %From Appendix A
    f5=f5*w;
    f5e=f5e*w;
    
    %Equation B5
    g1=(a-epsVal+n-1)*g1-risingFactorial(a-k,n+k-1)/gammaA;
    %Equation B6
    g2=(c-a+k+epsVal+n-1)*g2+risingFactorial(c-a+k,n-1)/gammaCA;
    %Equation B7
    g3=g3/(k+n+epsVal)+1/((k+n+epsVal)*factorial(k+n));
    %Equation B8
    g4=g4/(n-epsVal)+1/((n-epsVal)*factorial(n));
    %Equation B9
    g5=w*g5;
end

if(didSucceed==false)
    warning('The algorithm failed to converge.')
end
val=val+sumVal;

end


function val=caseIVa(a,b,c,w,k,epsVal,maxNumTerms)
%%CASEIVA This implements the special case IVa of [1], which occurs when
%         c-a-b is close to being zero.
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.

%First, calculate the finite sum in Equation 17.

if(k>0)
    termVal=gamma(k+epsVal);%The n=0 term.
    sumVal=termVal;
    for n=1:(k-1)
        termVal=termVal*(-(((-1+a+n)*(-1+b+n)*w)/(n*(k-n+epsVal))));
        sumVal=sumVal+termVal;
    end
    val=(gamma(c)/(gamma(c-a)*gamma(c-b)))*sumVal;
else
	val=0; 
end

gammaA=gamma(a);
gammaB=gamma(b);

%Next, compute the infinite sum in Equation 17.

%g values that do not depend on n. g6 and g7 are given by Equation 15.
%However, to reduce finite precision effects, the complicated method of
%Section III for computing such quantities is used via the
%invGammaForwardDiff function.
g6=invGammaForwardDiff(a+k,epsVal);
g7=invGammaForwardDiff(b+k,epsVal);

%f values for the n=0 case.
%From Appendix A
f1=1/gammaA;
%f1e=1/gammaA;%The f1 value offset by epsilon is not needed.
f2=1/gammaB;
f2e=1/gammaB;
%Equation B12
f3=-gamma(1);
f3e=-gamma(1-epsVal);
%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B10
    f3=f3/(k+n);
    f3e=f3e/(k+n+epsVal);
end

%Equation B13
f4=gamma(1);
f4e=gamma(1+epsVal);
%From Appendix A
f5=w^k;
f5e=w^(k+epsVal);

%These g values are from Appendix B
g1=0;
g2=0;

%As described in Section III of the text, finite differences of gamma
%functions can be difficult to compute to a high precision. Thus, instead
%of using 
%g3=-(gamma(1-epsVal)-gamma(1))/epsVal;
%from Appendix B or the limit g3=gamma(1)*polygamma(1); as epsVal->0,
%a solution utilizing Chebyshev polynomials is used. Note that said
%solution must be incremented up to time n=0 using Equation B16. For g4, a
%similar problem exists and again. Problems with using
%g4=(gamma(1+epsVal)-1)/epsVal exist with small values of g4, with the
%limit as epsVal->0 being g4=gamma(1)*polygamma(1). Thus Chebyshev
%interpolation is again used.
[g3,g4]=computeg3g4(epsVal);

%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B16
    g3=g3/(k+n+epsVal)+1/((k+n+epsVal)*factorial(k+n));
end

%A special approximation for the forward difference is needed to compute
%g5.
g5=powForwardDiff(k,w,epsVal);

sumVal=0;
n=0;
didSucceed=false;
outerCoeff=(-1)^k*gamma(c);
while(1)
    aPoch=risingFactorial(a,k+n);
    bPoch=risingFactorial(b,k+n);
    
    term2Add=-aPoch*f2*f3*f4e*f5*g6...
             -f1*bPoch*f3*f4e*f5*g7...
             -aPoch*bPoch*epsVal*f3*f4e*f5*g6*g7...
             +g1*f2e*f3e*f4*f5e...
             +f1*g2*f3e*f4*f5e...
             +f1*f2*g3*f4*f5e...
             -f1*f2*f3*g4*f5e...
             +f1*f2*f3*f4e*g5;
    term2Add=outerCoeff*term2Add;
    
    if(eps(sumVal)>abs(term2Add))
        didSucceed=true;
        break;
    end
    sumVal=sumVal+term2Add;
    
    n=n+1;
    
    if(n>maxNumTerms||~isfinite(sumVal))
        break;
    end
    
    %These definitions of f are from Appendix A.
    %f1=risingFactorial(a+k,n)/gammaA;
    f1=f1*(a+k+n-1);
    %The f1 value offset by epsilon is not needed.
    %f1e=risingFactorial(a+k+epsVal,n)/gammaA;
    %f1e=f1e*(a+k+epsVal+n-1);
    %f2=risingFactorial(b+k,n)/gammaB;
    f2=f2*(b+k+n-1);
    %f2e=risingFactorial(b+k+epsVal,n)/gammaB;
    f2e=f2e*(b+k+epsVal+n-1);
    %Equation B10
    f3=f3/(k+n);
    f3e=f3e/(k+n+epsVal);
    
    %Equation B11
    f4=f4/n;
    f4e=f4e/(n-epsVal);
    %From Appendix A
    f5=f5*w;
    f5e=f5e*w;
    
    %Equation B14
    g1=(a+k+epsVal+n-1)*g1+risingFactorial(a+k,n-1)/gammaA;
    %Equation B15
    g2=(b+k+epsVal+n-1)*g2+risingFactorial(b+k,n-1)/gammaB;
    %Equation B16
    g3=g3/(k+n+epsVal)+1/((k+n+epsVal)*factorial(k+n));
    %Equation B17
    g4=g4/(n-epsVal)+1/((n-epsVal)*factorial(n));
    %Equation B18
    g5=w*g5;
end

if(didSucceed==false)
    warning('The algorithm failed to converge.')
end
val=val+sumVal;

end

function val=caseIVb(a,b,c,w,k,epsVal,maxNumTerms)
%%CASEIVA This implements the special case IVb of [1], which occurs when
%         c-a-b is close to being zero.
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.

gammaA=gamma(a);
gammaB=gamma(b);
gammaC=gamma(c);

%First, calculate the finite sum in Equation 18.
if(k>0)
    termVal=w^(epsVal-k)*gamma(k-epsVal);%The n=0 term.
    sumVal=termVal;
    for n=1:(k-1)
        termVal=termVal*((-1-a+c+n)*(-1-b+c+n)*w)/(n*(-k+n+epsVal));
        sumVal=sumVal+termVal;
    end
    val=(gammaC/(gammaA*gammaB))*sumVal;
else
   val=0; 
end

%Next, compute the infinite sum in Equation 17.

%g values that do not depend on n. g6 and g7 are given by Equation 15.
%However, to reduce finite precision effects, the complicated method of
%Section III for computing such quantities is used via the
%invGammaForwardDiff function.
g6=invGammaForwardDiff(a-k,epsVal);
g7=invGammaForwardDiff(b-k,epsVal);

%f values for the n=0 case.
%From Appendix A
f1=risingFactorial(a-k,k)/gammaA;
%The f1 value offset by epsilon is not needed.
%f1e=risingFactorial(a-k+epsVal,k)/gammaA;
f2=risingFactorial(b-k,k)/gammaB;
f2e=risingFactorial(b-k+epsVal,k)/gammaB;

%Equation B21
f3=-gamma(1);
f3e=-gamma(1-epsVal);

%Equation B22
f4=gamma(1);
f4e=gamma(1+epsVal);
%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B20
    f4=f4/(k+n);
    f4e=f4e/(k+n-epsVal);
end
%From Appendix A
f5=1;
f5e=w^epsVal;

%From Appendix B
g1=0;
g2=0;
%We now have to loop to get the g1 and g2 terms up to n=0
for n=(-k+1):1:0
    %Equation B23
    g1=(a+epsVal+n-1)*g1+risingFactorial(a-k,k+n-1)/gammaA;
    %Equation B24
    g2=(b+epsVal+n-1)*g2+risingFactorial(b-k,k+n-1)/gammaB;
end

%As described in Section III of the text, finite differences of gamma
%functions can be difficult to compute to a high precision. Thus, instead
%of using 
%g3=-(gamma(1-epsVal)-gamma(1))/epsVal;
%from Appendix B or the limit g3=gamma(1)*polygamma(1); as epsVal->0,
%a solution utilizing Chebyshev polynomials is used. Problems with using
%g4=(gamma(1+epsVal)-1)/epsVal exist with small values of g4, with the
%limit as epsVal->0 being g4=gamma(1)*polygamma(1). Thus Chebyshev
%interpolation is again used. Note that equation B26 must eb used to
%increment g4 to the proper time.
[g3,g4]=computeg3g4(epsVal);

%We now have to loop to get the term up to n=0
for n=(-k+1):1:0
    %Equation B26
    g4=g4/(k+n-epsVal)+1/((k+n-epsVal)*factorial(k+n));
end

%A special approximation for the forward difference is needed to compute
%g5.
g5=powForwardDiff(k,w,epsVal);

sumVal=val;
n=0;
didSucceed=false;
outerCoeff=(-1)^k*gammaC;
while(1)
    aPoch=risingFactorial(a,n);
    bPoch=risingFactorial(b,n);
    
    term2Add=-aPoch*f2*f3*f4e*f5*g6...
             -f1*bPoch*f3*f4e*f5*g7...
             -aPoch*bPoch*epsVal*f3*f4e*f5*g6*g7...
             +g1*f2e*f3e*f4*f5e...
             +f1*g2*f3e*f4*f5e...
             +f1*f2*g3*f4*f5e...
             -f1*f2*f3*g4*f5e...
             +f1*f2*f3*f4e*g5;
    term2Add=outerCoeff*term2Add;
  
    if(eps(sumVal)>abs(term2Add))
        didSucceed=true;
        break;
    end
    sumVal=sumVal+term2Add;
    
    n=n+1;
    
    if(n>maxNumTerms||~isfinite(sumVal))
        break;
    end
    
    %These definitions of f are from Appendix A
    f1=risingFactorial(a-k,n+k)/gammaA;
    %The f1 value offset by epsilon is not needed.
    %f1e=risingFactorial(a-k+epsVal,n+k)/gammaA;
    f2=risingFactorial(b-k,n+k)/gammaB;
    f2e=risingFactorial(b-k+epsVal,n+k)/gammaB;
    
    %Equation B19
    f3=f3/n;
    f3e=f3e/(n+epsVal);
    
    %Equation B20
    f4=f4/(k+n);
    f4e=f4e/(k+n-epsVal);
    %From Appendix A
    f5=f5*w;
    f5e=f5e*w;
    
    %Equation B23
    g1=(a+epsVal+n-1)*g1+risingFactorial(a-k,k+n-1)/gammaA;
    %Equation B24
    g2=(b+epsVal+n-1)*g2+risingFactorial(b-k,k+n-1)/gammaB;
    %Equation B25
    g3=g3/(n+epsVal)+1/((n+epsVal)*factorial(n));
    %Equation B26
    g4=g4/(k+n-epsVal)+1/((k+n-epsVal)*factorial(k+n));
    %Equation B27
    g5=w*g5;
end

if(didSucceed==false)
	warning('The algorithm failed to converge.')
end
val=sumVal;

end

function [g3,g4]=computeg3g4(epsVal)
%%COMPUTEG3G4 This function implements the algorithm of Section III of [1]
%             for computing g4 and g4. The coefficients for the Chebychev
%             polynomials are not given in [1]. They were provided by the
%             author of [1].
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.

%The CII coefficients of Section III. The first element is not
%used.
c= [0.11528686913857579339872890819003657e1;
   -0.39836641427188668813550502856567435;
    0.16381491849746834445969671065563396;
   -0.41349972584595838242416447164595642e-1;
    0.11739888104509743948748485834561229e-1;
   -0.31509159742825717845846783104528302e-2;
    0.85084809366682540330028115184077086e-3;
   -0.22845443192182297253614554810213881e-3;
    0.61296656896858907270916323759970391e-4;
   -0.16433766723011959082591541534833589e-4;
    0.44046701847148520660258125028242579e-5;
   -0.11803851479587223345492859134791582e-5;
    0.31630339312403588488305625683201151e-6;
   -0.84755796666686117564957022251013564e-7;
    0.22710572677209079780536954678987573e-7;
   -0.60853209609268373214751556259951644e-8;
    0.16305620921375867864482570008163625e-8;
   -0.43690846345047718022878883179027790e-9;
    0.11706935476739890379554689241357534e-9;
   -0.31368649843198552351255033209421610e-10;
    0.84052057618382692960217222664957228e-11;
   -0.22521682699590609081199019088965996e-11;
    0.60346669123807723976181127096882828e-12;
   -0.16169841538137032176079290114309245e-12;
    0.43326960175123609635570088625382667e-13;
   -0.11609424034675431553315176322024985e-13;
    0.31107358004300087572452155428660087e-14;
   -0.83351914632193111475558815401948979e-15;
    0.22334078222557889355389486422061460e-15;
   -0.59843982246058550382747881611851515e-16;
    0.16035146716190080240936859943115090e-16;
   -0.42966046133076898235808019603294715e-17;
    0.11512717363557431988678458870224873e-17;
   -0.30848233202835882015258583966299712e-18;
    0.82657591746540727258216017499064442e-19;
   -0.22148034956862123422799663231945171e-19;
    0.59345480806145642339133686333296721e-20;
   -0.15901573656881585725893714030807897e-20;
    0.42608138203898096080539369435375448e-21;
   -0.11416816226321087557458906349840213e-21;
    0.30591266842950015571055286508657438e-22;
   -0.81969053674548061989664444282339330e-23;
    0.21963543471485197662543467891802004e-23;
   -0.58851140572211577956963471197095354e-24;
    0.15769121438531798083082131134888596e-24;
   -0.42253211944581570323425035302537635e-25;
    0.11321706791574145306428072576766804e-25;
   -0.30335842761477973373797446515125892e-26;
    0.81281383350578045680446098123885346e-27;
   -0.21782407988772728568103833180457024e-27;
    0.58395544064782062129754390403734767e-28;
   -0.15729062977489325257494410942884130e-28;
    0.42390612257722955199550993363196147e-29;
   -0.11242203351086692027388616387423238e-29;
    0.27892280419588143241883200553486195e-30;
   -0.75766427928255356179910217971637866e-31];

numCoeffs=length(c);

%Space for the Chebyshev polynomials
T3=zeros(numCoeffs,1);
T4=zeros(numCoeffs,1);

%Space for the polynomial differences of Section III.
F3=zeros(numCoeffs,1);
F4=zeros(numCoeffs,1);

x3=0;
T3(1)=1;
T3(2)=2*(x3-epsVal);
F3(1)=0;
F3(2)=-2;
g3=c(2)*F3(2);

x4=0;
T4(1)=1;
T4(2)=2*(x4+epsVal);
F4(1)=0;
F4(2)=2;
g4=c(2)*F4(2);

for i=2:(numCoeffs-1)
    T3(i+1)=4*(x3-epsVal)*T3(i)-T3(i-1);
    T4(i+1)=4*(x4+epsVal)*T4(i)-T4(i-1);
    F3(i+1)=-4*T3(i)+4*x3*F3(i)-F3(i-1);
    F4(i+1)=4*T4(i)+4*x4*F4(i)-F4(i-1);
    g3=g3+c(i+1)*F3(i+1);
    g4=g4+c(i+1)*F4(i+1);
end

g3=-g3;
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
