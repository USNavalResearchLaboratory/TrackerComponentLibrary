classdef VonMisesD
%%VONMISESD Functions to handle the von Mises distribution, which is a
%    common circular distribution. The distribution is given in Chapter
%    3.5.4 of [1].
%Implemented methods are: circVar, PDF, CDF,trigMoment, normProdDist,
%                         convDistApprox, params4TrigMoment, rand, entropy
%
%REFERENCES:
%[1] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester: John
%    Wiley and Sons, 2000.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=circVar(kappa)
%%CIRCVAR Obtain the circular variance of a particular von Mises
%         distribution.
%
%INPUTS: kappa The concentration of the distribution 0<=kappa<inf.
%
%OUTPUTS: val The circular variance of the distribution.
%
%The circular variance is defined as 1-abs(r), where r is the mean
%resultant value of the second moment. The mean resultant values for all
%moment of the von Mises distribution are given in Chapter 2.2.4 of [1].
%
%EXAMPLE:
%Here, we verify that the result agrees with the value obtained by
%sampling. 
% numRuns=1e5;
% kappa=100;
% mu=0.25;
% x=VonMisesD.rand([numRuns,1],mu,kappa);
% circVarEst=findCircVar(x);
% abs(VonMisesD.circVar(kappa)-circVarEst)./abs(VonMisesD.circVar(kappa))
%One will see that the relative error is on the order of 1e-3.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=1-BesseliRatio(1,kappa);

end

function val=PDF(x,mu,kappa)
%%PDF Evaluate a von Mises PDF at a given point with the specified
%     parameters.
%
%INPUTS: x The point(s) at which the von Mises distribution is to be
%          evaluated.
%       mu The mean direction of the distribution. This is normally between
%          -pi and pi.
%    kappa The concentration of the distribution 0<=kappa<inf.
%
%OUTPUTS: val The value(s) of the von Mises PDF with the specified
%             parameters evaluated at x.
%
%The von Mises distribution is used to model errors on angular quantities.
%It is described in [1] and [2].
%
%REFERENCES:
%[1] N. Fisher, Statistical Analysis of Circular Data. Cambridge, UK:
%    Cambridge University Press, 1993.
%[2] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester,
%    England: John Wiley & Sons, 2000.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=exp(kappa*cos(x-mu))/(2*pi*besseli(0,kappa));
end

function vals=CDF(x,mu,kappa,startAng)
%%CDF Evaluate  the cumulative distribution function (CDF) of the von Mises
%     distribution. Since this is a circular distribution, the CDF has to
%     be defined in terms of a specific starting angle. Increading density
%     goes in the direction of increasing angle from the starting point
%     (around the circle).
%
%INPUTS: z The point(s) at which the CDF should be evaluated.
%       mu The mean direction of the distribution. This is normally between
%          -pi and pi.
%    kappa The concentration of the distribution 0<=kappa<inf.
%  startAng The starting angle. This is where the CDF is zero. The default
%           if this parameter is omitted or an empty matrix is passed is
%           zero.
%
%OUTPUTS: vals The scalar value(s) of the von Mises CDF with the given PDF
%              evaluated at the points in z.
%
%The function Besseli0Inc gives the CDF of the von Mises distribution
%starting from -pi. Thus, by taking the difference of two evaluations of
%Besseli0Inc (and determining whether we are going the long or the short
%way around the circle) one can evaluate the CDF.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(startAng))
    startAng=-pi;
end

vals=Besseli0Inc(x-mu,kappa)-Besseli0Inc(startAng-mu,kappa);

sel=vals<0;
vals(sel)=vals(sel)+1;
end

function [rho,theta,R]=trigMoment(n,mu,kappa)
%%TRIGMOMENT Compute the nth trigonometric moment of a von Mises
%            distribution having the given parameters. A direction theta on
%            the unit circle can be modeled as a complex quantity having
%            unit magnitude as exp(1j*theta). The nth raw complex
%            trigonometric moment is the expected value of exp(1j*n*theta),
%            where the integral for the expected value is taken over any
%            interval of length 2*pi.
%
%INPUTS: n The order of the moment desired. This is >=1.
%       mu The mean direction of the distribution. This is normally between
%          -pi and pi.
%    kappa The concentration of the distribution 0<=kappa<inf. For
%          kappa>500, the results might be inaccurate.
%
%OUTPUTS: rho The (complex) mean resultant value for the nth moment. Note
%             that rho=R*exp(1j*theta).
%       theta The (real) trigonometric mean angle in radians for the nth
%             moment. This is between -pi and pi.
%           R The (real) mean resultant length for the nth moment.
%
%The moments of a von Mises distribution are given in Chapter 2.2.4 of [1]
%and Chapter 3.5.4 of [2].
%
%EXAMPLE:
%Here, we verify the first trigonometric moment based on samples.
% numRuns=1e5;
% kappa=100;
% mu=0.25;
% x=VonMisesD.rand([numRuns,1],mu,kappa);
% [rhoS,thetaS,RS]=findTrigMomentFromSamp(1,x);
% [rho,theta,R]=VonMisesD.trigMoment(1,mu,kappa);
% abs(rhoS-rho)
% abs(thetaS-theta)
% abs(RS-R)
%Here, we will see that the errors tend to be on the order of 1e-4,
%dominated by the angular error.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%[2] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester: John
%    Wiley and Sons, 2000.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

%For n=1 or 2, there should be no need to use anything but the default
%maximum number of iterations in BesseliRatio.
if(n==1)
    rho=BesseliRatio(1,kappa)*exp(1j*n*mu);
elseif(n==2)
    rho=BesseliRatio(2,kappa)*BesseliRatio(1,kappa)*exp(1j*n*mu);
else
    rho=besseli(n,kappa,1)/besseli(0,kappa,1)*exp(1j*n*mu);
end

if(nargout>1)
    theta=angle(rho);
    R=abs(rho);
end
end

function [mu,kappa]=normProdDist(mu1,kappa1,mu2,kappa2)
%%NORMPRODDIST The product of two von Mises distributions is an
%              unnormalized von Mises distribution. This finds the
%              parameters of the product distribution.
%
%INPUTS: mu1, kappa1 The mean direction and concentration of the first
%                    distribution. Note that 0<=kappa<inf.
%        mu2, kappa2 The mean direction and concentration of the second
%                    distribution.
%
%OUTPUTS: mu,kappa The mean direction and concentration of the product
%                  distribution -pi<=mu<pi.
%
%The product of two von Mises distributions is described in [1].
%
%REFERENCES:
%[1] M. Azmani, S. Reboul, J.-B. Choquel, and M. Benkelloun, "A recursive
%    fusion filter for angular data," in Proceedings of the 2009 IEEE
%    International Conference on Robotics and Biomimetics, Guilin, China,
%    19-23 Dec. 2009, pp. 882-887.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

C=kappa1*cos(mu1)+kappa2*cos(mu2);
S=kappa1*sin(mu1)+kappa2*sin(mu2);
mu=atan2(S,C);
kappa=sqrt(C^2+S^2);
end

function [mu,kappa]=convDistApprox(mu1,kappa1,mu2,kappa2)
%%CONVDISTAPPROX The convolution of two von Mises distributions is not a
%          von Mises distribution. However, it can be approximated by a von
%          Mises distribution that matches the moments, as is done in
%          Equation 7 of [1]. This function returns the paramters of the
%          approximation.
%  
%INPUTS: mu1, kappa1 The mean direction and concentration of the first
%                    distribution. Note that 0<=kappa<inf.
%        mu2, kappa2 The mean direction and concentration of the second
%                    distribution.
%
%OUTPUTS: mu, kappa The mean direction and concentration of a von Mises
%                   approximation to the product distribution. -pi<=mu<pi.
%
%REFERENCES:
%[1] I. Markovic and I. Petrovic, "Bearings-only tracking with a mixture of
%    von Mises distributions," in IEEE/RSJ International Conference on
%    Intelligent Robots and Systems, Vilamoura, Algarve, Portugal, 7-12
%    Oct. 2012, pp. 707-712.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
    mu=wrapRange(mu1+mu2,-pi,pi);

    AkProd=BesseliRatio(1,kappa1)*BesseliRatio(1,kappa2);
    kappa=BesseliRatioInv(1,AkProd);
end


function [mu,kappa]=params4TrigMoment(rho)
%%PARAMS4TRIGMOMENT Obtain the parameters (mu and kappa) of the von Mises
%           distribution such that the von Mises distribution matches the
%           given complex mean resultant length for the first moment. 
%
%INPUTS: rho The (complex) mean resultant value for the first moment. Note
%            that abs(rho)<=1.
%
%OUTPUTS: mu, kappa The mean direction and the concentration parameter of
%                   the von Mises distribution with the specified mean
%                   resultant value.
%
%The moments of a von Mises distribution are given in Chapter 2.2.4 of [1].
%This function just implements an inverse.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    mu=atan2(imag(rho),real(rho));
    kappa= BesseliRatioInv(1,abs(rho));
end

function vals=rand(N,mu,kappa)
%%RAND Generate von Mises random variables with the given parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       mu The mean direction of the distribution. This is normally
%          between -pi and pi.
%    kappa The concentration of the distribution 0<=kappa<inf.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated von Mises random variables.
%
%The algorithm if that of Section 4 of [1], which is a rejection sampling
%method based on using a wrapped Cauchy distribution as a proposal density.
%However, the range of the random variable u1 in [1] does not appear to be
%correct. The correct formulation is based on Chapter 7.3 of [2] using the
%errata available at http://luc.devroye.org/errors.pdf so u1 is randomly
%generated between -1 and 1 as opposed to between 0 and 1.
%
%REFERENCES:
%[1] D. J. Best and N. I. Fisher, "Efficient simulation of the von Mises
%    distribution," Journal of the Royal Statistical Society. Seriec C
%    (Applied Statistics), vol. 28, no. 2, pp. 152-157, 1979.
%[2] L. Devroye, Non-Uniform random Variate Generation. New York:
%    Springer-Verlag, 1986.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end

    vals=zeros(dims);%Allocate space
    numSamp=prod(dims);

    %Step 0
    tau=1+sqrt(1+4*kappa^2);
    rho=(tau-sqrt(2*tau))/(2*kappa);
    r=(1+rho^2)/(2*rho);
    
    for curSamp=1:numSamp
        while(1)
            %Step 1
            u1=2*rand(1)-1;%In range [-1,1]
            z=cos(pi*u1);
            f=(1+r*z)/(r+z);
            c=kappa*(r-f);

            %Step 2
            u2=rand(1);%in range [0,1]
            if(c*(2-c)-u2>0)
               break;
            end

            %Step 3
            if(log(c/u2)+1-c>=0)
                break;
            end
        end
        %Step 4
        u3=2*rand(1)-1;%In range [-1,1]
        vals(curSamp)=wrapRange(sign(u3)*acos(f)+mu,-pi,pi);
    end
end
    
function entropyVal=entropy(kappa)
%%ENTROPY Obtain the differential entropy of the von Mises distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x)) dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: kappa The concentration of the distribution 0<=kappa<inf.
%
%OUTPUTS: entropyVal The value of the differential entropy in nats.
%
%Differential entropy is defined in Chapter 8 of [1].
%
%If kappa is very large, besseli(0,kappa) will overflow and a warning will
%be issued. In such an instance, an approximation based on the first term
%of the series on page 289 of [2] is used. The approximation becomes
%increasingly bad as kappa becomes larger.
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%[2] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    besseliVal=besseli(0,kappa);
    
    if(isfinite(besseliVal))
        entropyVal=-kappa*BesseliRatio(1,kappa)+log(2*pi*besseliVal);
    else
        %Use the approximation exp(kappa)/sqrt(2*pi*kappa) for values of
        %kappa so large that besseli overflows
        warning('besseli(0,kappa) overflowed; using an approximation.')
        entropyVal=-kappa*BesseliRatio(1,kappa)+(1/2)*(4000+log(pi/1000));
    end
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
