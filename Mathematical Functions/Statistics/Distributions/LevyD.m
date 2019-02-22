classdef LevyD
%%LEVYD Functions to handle the Lévy distribution. This is a one-sided
%       distribution with a very long tail. It arises when considering the
%       amount of time it takes a Brownian motion drift to hit a
%       particular point.
%Implemented methods are: PDF, CDF, invCDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=PDF(x,mu,c)
%%PDF Evaluate the probability density function (PDF) of the Lévy
%     distribution.
%
%INPUTS: x The point or a matrix of points at which the PDF of the Lévy
%          distribution should be evaluated. x>=mu for nonzero values.
%       mu The location parameter of the distribution.
%        c The scale parameter of the distribution.
%
%OUTPUTS: val The value(s) of the PDF evaluated at the point(s) in x.
%
%The Lévy distribution is presented in Section 2.1 of [1]. The version
%implemented here allows the zero point to be shifted with the location
%parameter mu.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=4;
% c=1/10;
% numSamples=1e6;
% figure(1)
% clf
% histogram(LevyD.rand([numSamples,1],mu,c),'Normalization','pdf','BinLimits',[4,10])
% hold on
% numPoints=1000;
% x=linspace(4,20,numPoints);
% vals=LevyD.PDF(x,mu,c);
% plot(x,vals,'linewidth',2)
% axis([4,10,0,4])
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] W. A. Woyczynski, "Lévy Processes in the Physical Sciences," in Levy
%    Processes: Theory and Applications, O. E. Barndorff-Nielson, T. 
%    Mikosch, and S. I. Resnick, Eds. Boston: Birkhäuser, 2001.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=sqrt(c/(2*pi))*exp(-c./(2*x-mu))./(x-mu).^(3/2);%(sqrt(x-mu).^3);
    val(x<mu)=0;
end
    
function prob=CDF(x,mu,c)
%%CDF Evaluate the cumulative distribution function (CDF) of the Lévy
%     distribution.
%
%INPUTS: x The point or a matrix of points at which the CDF of the Lévy
%          distribution should be evaluated.x>mu for nonzero values.
%       mu The location parameter of the distribution.
%        c The scale parameter of the distribution.
%
%OUTPUTS: prob The value(s) of the CDF of the Lévy distribution evaluated
%              at x.
%
%The Lévy distribution is presented in Section 2.1 of [1]. The version
%implemented here allows the zero point to be shifted with the location
%parameter mu. The CDF was just obtained from the integral of the PDF.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=4;
% c=1/10;
% x=5;
% numSamples=1e6;
% prob=LevyD.CDF(x,mu,c)
% probSamp=mean(LevyD.rand([numSamples,1],mu,c)<=x)
%One will see that both values are about 0.7518.
%
%REFERENCES:
%[1] W. A. Woyczynski, "Lévy Processes in the Physical Sciences," in Levy
%    Processes: Theory and Applications, O. E. Barndorff-Nielson, T. 
%    Mikosch, and S. I. Resnick, Eds. Boston: Birkhäuser, 2001.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=erfc(sqrt(c./(2*(x-mu))));
    prob(x<mu)=0;
end

function x=invCDF(prob,mu,c)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the Lévy distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          mu The location parameter of the distribution.
%           c The scale parameter of the distribution.
%
%OUTPUTS: x The argument(s) of the CDF that would give the probability or
%           probabilities in prob.
%
%The Lévy distribution is presented in Section 2.1 of [1]. The version
%implemented here allows the zero point to be shifted with the location
%parameter mu. The CDF was just obtained from the integral of the PDF and
%is in terms of an error function. Thus, the inverse error function
%provides the inverse CDF value.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% mu=4;
% c=1/10;
% x=5;
% xBack=LevyD.invCDF(LevyD.CDF(x,mu,c),mu,c)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] W. A. Woyczynski, "Lévy Processes in the Physical Sciences," in Levy
%    Processes: Theory and Applications, O. E. Barndorff-Nielson, T. 
%    Mikosch, and S. I. Resnick, Eds. Boston: Birkhäuser, 2001.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

   x=mu+c./(2*erfcinv(prob).^2);
end

function x=rand(N,mu,c)
%%RAND Generate Lévy distributed random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       mu The location parameter of the distribution.
%        c The scale parameter of the distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Lévy random variables.
%
%This function implements the inverse transformation method of Chapter 5.1
%of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end
    
    U=rand(dims);
    x=LevyD.invCDF(U,mu,c);
end

function entropyVal=entropy(c)
%%ENTROPY Obtain the differential entropy of the Lévy distribution given in
%         nats. The differential entropy of a continuous distribution is
%         entropy=-int_x p(x)*log(p(x)) dx where the integral is over all
%         values of x. Units of nats mean that the natural logarithm is
%         used in the definition. Unlike the Shannon entropy for discrete
%         variables, the differential entropy of continuous variables can
%         be both positive and negative.
%
%INPUTS: c The scale parameter of the distribution.
%
%OUTPUTS: entropyVal The value of the differential entropy in nats.
%
%Differential entropy is defined in Chapter 8 of [1].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    gammaVal=-psi(1);
    entropyVal=(1+3*gammaVal+log(16*pi*c^2))/2;
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
