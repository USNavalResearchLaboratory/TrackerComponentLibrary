function Q=MarcumQ(mu,alpha,beta)
%%MARCUMQ  Evaluate the generalized Marcum Q function of a desired order,
%          Q_M(alpha,beta) including when the order is not an integer.
%
%INPUTS: mu The scalar order (mu>0) of the Marcum Q function. The function
%           is most stable if mu is an integer or if alpha is zero if mu is
%           not an integer.
%     alpha The scalar first parameter of the Q function (alpha>=0).
%      beta The scalar second parameter of the Q function (beta>=0).
%
%OUTPUTS: Q The value of the generalized Marcum Q function with the given
%           parameters.
%
%DEPENDENCIES: GammaD.m
%
%Note that the parameter order is different from that used in the marcumq
%function that is part of Matlab's statistics toolbox. Also note that
%Matlab's function cannot handle non-integer orders at all.
%
%This function tries to combine the best implementations of the MarcumQ
%function for a wide range of parameters. The most reliable results are
%obtained if mu is an integer or alpha is zero. The algorithms used are
%specifically,
%- When alpha=0, the MarcumQ function becomes a regularized gamma function.
%  Thus, the build-in gammainc function can be used for a very wide range
%  of (non-integer) values of mu and beta.
%- When alpha is not zero, but mu is a multiple of 0.5, then
%  the algorithm of [1] is used.
%- If alpha is not zero and mu is not a multiple of 0.5, then the
%  equivalency of the Marcum Q function with the noncentral chi square CDF
%  and therefore the noncentral gamma CDF is used to compute the result.
%
%Note that a result from the [2] was used to modify the implementation of
%Benton's algorithm in [1] for large input values.
%
%Also, for those looking to understand the difficulty of computing the
%MarcumQ function, the following algorithms are implemented as
%separate functions in this file but are not used, because they are
%inferior to the algorithms chosen. The algorithms are not
%used either because they are not as stable or they are slower. The
%algorithms are taken from [3], [4], [5], [6], and [7].
%
%REFERENCES:
%[1] D. Benton and K. Krishnamoorthy, "Computing Discrete Mixtures of
%    Continuous Distributions: Noncentral Chisquare, Noncentral t and the
%    Distribution of the Square of the Sample Multiple Correlation
%    Coefficient," Computational Statistics & Data Analysis, vol. 43, no.
%    2, pp.249-26, 28 Jun. 2003.
%[2] D. A. Shnidman, "Note on "The Calculation of the Probability of
%    Detection and the Generalized Marcum Q-Function"," IEEE Transactions
%    on Information Theory, vol. 37, no. 4. pg. 1233, Jul. 1991.
%[3] A. Annamalai Jr. and C. Tellambura, "A Simple Exponential Integral
%    Representation of the Generalized Marcum Q-Function Q_M(a,b) for
%    Real-Order M with Applications," in Proceedings of the IEEE Military
%    Communications Conference, San Diego, CA, 16-19 Nov. 2008.
%[4] A. Annamalai Jr., A., C. Tellambura, and J. Matyjas, "A New Twist on
%    the Generalized Marcum Q-Function Q_M(a,b) with Fractional-Order M and
%    Its Applications," in Proceedings of the IEEE Consumer Communications
%    and Networking Conference, Las Vegas, NV, 10-13 Jan. 2009.
%[5] C. G. Ding, "Algorithm AS 275: Computing the Non-Central Chi^2
%    Distribution Function," Journal of the Royal Statistical Society
%    Series C (Applied Statistics) vol. 41. no. 2. pp. 478-482, 1992.
%[6] A. H. Ross, "Algorithm for Calculating the Noncentral Chi-Square
%    Distribution," IEEE Transactions on Information Theory, vol. 45, no. 4,
%    pp. 1327-1333, May 1999.
%[7] D. A. Shnidman, "The Calculation of the Probability of Detection and
%    the Generalized Marcum Q-Function," IEEE Transactions on Information 
%    Theory, vol. 35, no. 2, pp. 389-400, Mar. 1989.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%December 2014 David A. Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(mu<=0||alpha<0||beta<0||~isreal(mu)||~isreal(alpha)||~isreal(beta))
    error('Invalid Inputs');
end

%Check for the special case of a zero alpha. This will work for 
if(alpha==0)
%If alpha is zero, the MarcumQ function reduces to a regularized gamma
%function. This identity is programmed into Mathematica.
    Q=1-gammainc(beta^2/2,mu);
    return
end

%Use Ross' algorithm if mu is of the form integer+0.5 and the precision
%bounds are satisfied.
diff=mu-fix(mu);
if((diff==0.5)||(diff==0))
    Q=AlgorithmBenton(mu,alpha,beta);
    return
end

%If mu is not a multiple of 0.5, then use equivalency of the Marcum Q
% function with the noncentral chi square CDF and therefore the noncentral
% gamma CDF to compute the result.
    Q=1-GammaD.CDF(beta.^2,mu,2,alpha^2);
end

function Q=AlgorithmSchnidman(mu,alpha,beta)
%This implements the algorithm of [1] but with the improved exponential
%term of [2].

%To determine whether to compute Q or 1-Q, the decision criterion used in
%[3] is employed here.
%
%REFERENCES:
%[1] D. A. Shnidman, "The Calculation of the Probability of Detection and
%    the Generalized Marcum Q-Function," IEEE Transactions on Information
%    Theory, vol. 35, no. 2, pp. 389-400, Mar 1989.
%[2] D. A. Shnidman, "Note on "The Calculation of the Probability of
%    Detection and the Generalized Marcum Q-Function"," IEEE Transactions
%    on Information Theory, vol. 37, no. 4. pg. 1233, Jul. 1991.
%[3] A. H. Ross, "Algorithm for Calculating the Noncentral Chi-Square
%    Distribution," IEEE Transactions on Information Theory, vol. 45, no. 4,
%    pp. 1327-1333, May 1999.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    N=mu;
    %X and Y are related to alpha and beta as in Equation 7.
    X=alpha^2/(2*N);
    Y=beta^2/2;
    NX=N*X;

    %Compute the Chernoff bound. This is necessary to determine which
    %algorithm to use.
    %Equation 22
    lambda=1-N/Y-sqrt((N/(2*Y))^2+NX/Y);
    %Equation 23, the Chernoff bound.
    CB=exp(-lambda*Y+NX*lambda/(1-lambda)-N*log(1-lambda));
    
    %The paper describes a way of avoiding computation if the Chernoff
    %bound is below some minimum value (i.e. the maximum amount of error we
    %want to allow). Here, the minimum value is chosen to be eps.
    if(CB<eps)
        if(Y>NX+N)
            Q=0;
        else
            Q=1;
        end
        return;
    end
    
    %Ross's criterion for determining whether to find Q or 1-Q.
    if(alpha^2+mu<beta^2)
        %Use the series to find the Q value directly.
        %This is form 2 in the paper. Form two consists of two sums. The
        %first sum is fairly simply to evaluate. One need just determine 
        %the first value in the summation that should be used.
        mMin=sumStart(Y);
        
        m=mMin:(N-1);
        if(~isempty(m))
           Q=sum(expFunc(Y,m));
        else
            Q=0;
        end
        
        %Now, the double-sum term must be evaluated.
        kMin=mMin;
        mMin=max(mMin,N);%The outer sum is from mMin to infinity.
        
        kInit=kMin:((mMin-1)-N);
        if(~isempty(kInit))
            innerSum=sum(expFunc(NX,kInit));
        else
            innerSum=0;
        end
        
        %Sum until the term no longer chnges the sum.
        m=mMin;
        deltaQ=Inf;
        while(deltaQ>eps(Q))
            innerSum=innerSum+expFunc(NX,m-N);
            
            deltaQ=expFunc(Y,m)*(1-innerSum);
            Q=Q+deltaQ;
            
            m=m+1;
        end
    else
        %Use the  series to find 1-Q.
        %This is form 4 in the paper. Form 4 is similar to just the second
        %summation term in form 2. Thus, a similar method of evaluating the
        %sums is employed here.
        
        kMin=sumStart(Y);
        mMin=kMin;
    
        mInit=mMin:(N-1+kMin-1);
        if(~isempty(mInit))
            innerSum=sum(expFunc(Y,mInit));
        else
            innerSum=0;
        end
        
        P=0;
        %Sum until the term no longer chnges the sum.
        k=kMin;
        deltaP=Inf;
        while(deltaP>eps(P))
            innerSum=innerSum+expFunc(Y,N-1+k);
            
            deltaP=expFunc(NX,k)*(1-innerSum);
            P=P+deltaP;
            
            k=k+1;
        end
        Q=1-P;
    end
    
    %Section 5 of the paper describes how to determine when the sums can
    %start. If L(m,z)=sum_{n=0}^m exp(-z)*z^n/n!, we want to throw out the
    %first M-1 terms that are so small they do not contribute any
    %meaningful amount.
    function M=sumStart(z)
        %An arbitrary, small amount. The paper provides some guidelines for
        %choosing this quantity, but they are rather vague.
        eM=10^(-30);
        
        %Equation 41
        G=-4*log(4*eM*(1-eM))/5;
        sqrtArg=G*(2*z-G);
        
        %If the square root would have been imaginary, then the sum should
        %start at the first element.
        if(sqrtArg<0)
            M=0;
        else
            %Equation 40
            M=fix(z+0.5-sqrt(sqrtArg));
        end
    end
end

function Q=AlgorithmDing(mu,alpha,beta)
%This is the algorithm of [1]. This will work for small values of alpha and
%values of mu that are either integers or an integer+0.5. However, if mu is
%large, then this algorithm can be slow or suffer from overflow/underflow
%errors.
%
%REFERENCES:
%[1] C. G. Ding, "Algorithm AS 275: Computing the Non-Central Chi^2
%    Distribution Function," Journal of the Royal Statistical Society
%    Series C (Applied Statistics) vol. 41. no. 2. pp. 478-482, 1992.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    x=beta^2;
    nu=2*mu;
    lambda=alpha^2;

    s=exp(-lambda/2);
    u=s;
    
    t=(1/gamma(nu/2+1))*(x/2)^(nu/2)*exp(-x/2);
    P=s*t;
    i=1;
    truncErrBound=Inf;
    while(truncErrBound>eps(P))
        %This is from Equation 7 in Ding's paper as the error bound for
        %truncating the series at this term. As in the paper, this will
        %take at least 2*nu iterations, as nu+2*i<x otherwise and the
        %denominator is negative. Because of the large minimum number of
        %iterations for big nu, a different algorithm is preferable when nu
        %becomes large.
        truncErrBound=t*x/(nu+2*i-x);
        if(truncErrBound<0)
            truncErrBound=Inf;
        end
        
        t=t*x/(nu+2*i);
        u=u*lambda/(2*i);
        s=s+u;
        P=P+s*t;
        
        i=i+1;
    end
    
    Q=1-P;
end

function Q=AlgorithmBenton(mu,alpha,beta)
%This function implements the algorithm of [1], whose code is given in
%Section 7.3 as Algorithm 7.3. Some of the comments
%are taken from the paper. The algorithm will evaluate the MarcumQ
%function at integer values of mu as well as values of mu that are integers
%+0.5. Changes from the original code are the substitution of expFunc(y,M)
%from [2] for the code that computed terms of the form exp(-y)*y^M/M! to
%avoid overflow/underflow errors. Also, the error tolerance is set to
%eps(sumVal) and some code for the fact that we want Q and not 1-Q is
%added as well as the translation of the input parameters
%into the form for the noncentral Chi-Squared distribution.
%
%REFERENCES:
%[1] D. Benton and K. Krishnamoorthy, "Computing Discrete Mixtures of
%    Continuous Distributions: Noncentral Chisquare, Noncentral t and the
%    Distribution of the Square of the Sample Multiple Correlation
%    Coefficient," Computational Statistics & Data Analysis, vol. 43, no.
%    2, pp.249-26, 28 Jun. 2003.
%[2] D. A. Shnidman, "Note on "The Calculation of the Probability of
%    Detection and the Generalized Marcum Q-Function"," IEEE Transactions
%    on Information Theory, vol. 37, no. 4. pg. 1233, Jul. 1991.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

maxIter=5000;

n=2*mu;
lambda=alpha^2;
y=beta^2;

x=y/2;
del=lambda/2;
k=fix(del);
a=n/2+k;
%Compute the gamma distribution function using (4.5) at (x; a), and assign
%it to "gamkf" and "gamkb" so that they can be called laterforfor ward as
%well as backward computations:
gamkf=gammainc(x,a);
gamkb=gamkf;

if(lambda==0)
   P=gamkf;
   Q=1-P;
   return
end

%Compute the Poisson probability at (k; del) and assign it to "poikf" and 
%"poikb" so that they can be used as initial values for forward and
%backward recursions:
poikf=expFunc(del,k);
poikb=poikf;
%"xtermf" is an initialization to compute the second term in (4.3)
%recursively:
xtermf=expFunc(x,a-1);
%"xtermb" is an initialization to compute the second term in (4.4)
%recursively:
xtermb=xtermf*x/a;
sumVal=poikf*gamkf;
remain=1-poikf;
i=1;

while(1)
    xtermf=xtermf*x/(a+i-1);
    gamkf=gamkf-xtermf;
    poikf=poikf*del/(k+i);
    sumVal=sumVal+poikf*gamkf;
    errorVal=remain*gamkf;
    remain=remain-poikf;
    
%Do forward and backward computations k times or until convergence:
    if(i>k)
        if(errorVal<=eps(sumVal))
            break;
        elseif(i>maxIter)
            warning('Maximum number of iterations Exceeded. Unable to meet the error tolerance.')
            break;
        else
           i=i+1; 
        end
    else
        xtermb=xtermb*(a-i+1)/x;
        gamkb=gamkb+xtermb;
        poikb=poikb*(k-i+1)/del;
        sumVal=sumVal+gamkb*poikb;
        remain=remain-poikb;
        if(remain<=eps(sumVal))
            break;
        elseif(i>maxIter)
            warning('Maximum number of iterations Exceeded. Unable to meet the error tolerance.')
            break;
        else
           i=i+1; 
        end
    end
end

P=sumVal;
%The max deals with negative values within precision bounds.
Q=max(1-P,0);
%The min deals with positive values within precision bounds.
Q=min(Q,1);
end

function Q=AlgorithmAnnamalai(mu,alpha,beta)
%This function implements the algorithm of [1]., which will work with
%non-integer values of mu, but which will produce bad results if
%(alpha/beta)^(1-mu) is large. The same algorithm is also presented in [2].
%
%The algorithm is just an implementation of equations in the paper.
%However, a slight ad-hoc modification to deal with numerical precision
%limitations near singular points has been added to the integral function.
%
%REFERENCES:
%[1] A. Annamalai Jr. and C. Tellambura, "A Simple Exponential Integral
%    Representation of the Generalized Marcum Q-Function Q_M(a,b) for
%    Real-Order M with Applications," in Proceedings of the IEEE Military
%    Communications Conference, San Diego, CA, 16-19 Nov. 2008.
%[2] A. Annamalai Jr., A., C. Tellambura, and J. Matyjas, "A New Twist on
%    the Generalized Marcum Q-Function Q_M(a,b) with Fractional-Order M and
%    Its Applications," in Proceedings of the IEEE Consumer Communications
%    and Networking Conference, Las Vegas, NV, 10-13 Jan. 2009.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    zeta=alpha/beta;
    %THis is from Equation
    intFun=@(theta)intFuncQ(theta);
    H=(zeta^(1-mu)/pi)*integral(intFun,0,pi);

    %This is Equation 6.
    intFunC=@(x)(x./zeta).^(mu-1)./(zeta+x).*exp(-0.5*(alpha^2+beta^2+alpha*beta*(x+1./x)));
    C=(sin(mu*pi)/pi)*integral(intFunC,0,1);
    
    %The cases are from Equation 5.
    if(alpha<beta)
        Q=H+C;
    elseif(alpha==beta)
        Q=0.5+H+C;
    else
        Q=1+H+C;
    end
    
    if(Q<0)
        warning('Precision limitations have made the result unreliable');
        Q=0;
    end
    
    if(Q>1||~isfinite(Q))
        warning('Precision limitations have made the result unreliable');
        Q=1;
    end
    
    function val=intFuncQ(theta)
    val=exp(-0.5*(alpha^2+beta^2-2*alpha*beta*cos(theta))).*(cos((mu-1)*theta)-zeta*cos(mu*theta))./(1+zeta^(2)-2*zeta*cos(theta));
    
    %An ad-hoc fix for endpoints of integration that might run into numeric
    %precision problems.
    thetaSingular=theta(~isfinite(val));
    val(~isfinite(val))=exp(-0.5*(alpha^2+beta^2-2*alpha*beta*cos(thetaSingular))).*(cos((mu-1)*thetaSingular)-zeta*cos(mu*thetaSingular))./(1e-10);
    end
end


function Q=AlgorithmRoss(mu,alpha,beta)
%This function implements the algorithm of [1] to solve the MarcumQ
%function for values of mu>=1.5 that are integers or integers +0.5. Due to
%the precision limitations, the algorithm should not be used if
%log(beta/alpha) < (log(realmax())-log(6/eps))/(mu-1), as is documented in
%the paper.
%
%The implementation comes directly from implementing all of the equations
%in Table II in the paper to get Q or 1-Q.
%
%REFERENCES:
%[1] A. H. Ross, "Algorithm for Calculating the Noncentral Chi-Square
%    Distribution," IEEE Transactions on Information Theory, vol. 45, no. 4,
%    pp. 1327-1333, May 1999.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    m=mu*2;
    a=alpha;
    b=beta;
    
    if(mod(m,2)==0)
        if(a^2+m<b^2)
            r=a/b;
            fv=1;
        else
            r=b/a;
            fv=0;
        end

        tv=1;
        fPrev=0;
        q=2;
        v=1;
        s=exp(-(a-b)^2/2);

        gv=1;
        gPrev=0;
        DenTest=(3/eps)*(1+sqrt(pi*a*b/2));
        NumTest1=r*DenTest;
        NumTest2=(6/eps)*(b/a)^(m/2-1);

        if(a^2+m<b^2)
            while(v<=m/2||fv<=tv*NumTest1||fv<=NumTest2||gv<=DenTest)
                tv=r*tv;
                if(v<m/2)
                    pv=tv+1/tv; 
                else
                    pv=tv;
                end

                temp=fv;
                fv=pv+(2*v/(a*b))*fv+fPrev;
                fPrev=temp;

                temp=gv;
                gv=q+(2*v/(a*b))*gv+gPrev;
                gPrev=temp;

                v=v+1;
            end

            Q=s*fv/gv;
        else
            while(v<=m/2||fv<=tv*NumTest1||fv<=NumTest2||gv<=DenTest)
                tv=r*tv;
                if(v<m/2)
                    pv=0; 
                else
                    pv=tv;
                end

                temp=fv;
                fv=pv+(2*v/(a*b))*fv+fPrev;
                fPrev=temp;

                temp=gv;
                gv=q+(2*v/(a*b))*gv+gPrev;
                gPrev=temp;
                v=v+1;
            end

            Q=1-s*fv/gv;
        end
    else
        if(a^2+m<b^2)
            r=a/b;
            fv=sqrt(r);
        else
            r=b/a;
            fv=0;
        end

        tv=sqrt(r);
        fPrev=0;
        q=0;
        v=3/2;
        s=(1/sqrt(2*pi*a*b))*(1-exp(-2*a*b))*exp(-(a-b)^2/2);

        gv=1;
        gPrev=0;
        DenTest=(3/eps)*(1+sqrt(pi*a*b/2));
        NumTest1=r*DenTest;
        NumTest2=(6/eps)*(b/a)^(m/2-1);

        if(a^2+m<b^2)
            while(v<=m/2||fv<=tv*NumTest1||fv<=NumTest2||gv<=DenTest)
                tv=r*tv;
                if(v<m/2)
                    pv=1/tv; 
                else
                    pv=0;
                end

                temp=fv;
                fv=pv+(2*v/(a*b))*fv+fPrev;
                fPrev=temp;

                temp=gv;
                gv=q+(2*v/(a*b))*gv+gPrev;
                gPrev=temp;

                v=v+1;
            end

            Q=(0.5-0.5*erf((b+a)/sqrt(2)))+(0.5-0.5*erf((b-a)/sqrt(2)))+s*fv/gv;
        else
            while(v<=m/2||fv<=tv*NumTest1||fv<=NumTest2||gv<=DenTest)
                tv=r*tv;
                if(v<m/2)
                    pv=0; 
                else
                    pv=tv;
                end

                temp=fv;
                fv=pv+(2*v/(a*b))*fv+fPrev;
                fPrev=temp;

                temp=gv;
                gv=q+(2*v/(a*b))*gv+gPrev;
                gPrev=temp;

                v=v+1;
            end

            Q=1-s*fv/gv;
        end
    end
end


function val=expFunc(y,M)
%This function evaluates exp(-y)*y^M/M! while trying to avoid
%underflows. This function is used by [1] when imp[lementing the MarcumQ
%function as well as in the function AlgorithmBenton in this file.
%
%REFERENCES:
%[1] D. A. Shnidman, "Note on "The Calculation of the Probability of
%    Detection and the Generalized Marcum Q-Function"," IEEE Transactions
%    on Information Theory, vol. 37, no. 4. pg. 1233, Jul. 1991.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(y>0)
        if(y<10000)
            val=exp(-y+M*log(y)-gammaln(M+1));
        else
            %This is the expression from the 1991 Shidman paper. It is
            %supposed to be particularly accurate when y is very large.
            %Note that in the expression for A in the paper, the terms
            %z+1/2 and 1+1/(2*z) have the wrong sign in front of the z's.
            %It is not always used, because it is not as precise when y is
            %not large.
            z=M+1;
            A=(z-0.5)*((1-y/z)/(1-1/(2*z))+log(y/z))-0.5*log(2*pi*y)-J(z);
            val=exp(A);
        end
    else
        %Deal with the zero exponent case.
        val=0+(M==0);
    end
    
    function val=J(z)
    %This is the Binet series from the Shidman papers. In the event
    %that one moves up to higher precision arithmetic, it might be
    %necessary to add more terms.
    val=1./(12*z+2./(5*z+53./(42*z+1170./(53*z+53./z)))); 
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


