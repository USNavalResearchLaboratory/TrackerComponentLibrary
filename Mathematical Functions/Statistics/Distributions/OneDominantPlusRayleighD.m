classdef OneDominantPlusRayleighD
%%ONEDOMINANTPLUSRAYLIEGHD Functions to handle the one-dominant-plus-
%    Rayleigh distribution. This distribution arises in the Swerling III
%    and IV models in Chapter 11 of [1].
%
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.   
    
methods(Static)
function val=mean(A0)
%%MEAN Obtain the mean of the one-dominant-plus-Rayleigh distribution.
%
%INPUTS: A0 The parameter of the distribution. This is the most probable
%           value of the distribution.
%
%OUTPUTS: val The mean of the one-dominant-plus-Rayleigh distribution.
%
%The distribution is defined in Equation 11.4-1 of [1]. The mean was found
%by integration.
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% A0=2;
% meanVal=OneDominantPlusRayleighD.mean(A0)
% numSamp=1e6;
% sampMeanVal=mean(OneDominantPlusRayleighD.rand([numSamp,1],A0))
%One will find both values are about 2.170.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(1/2)*A0*sqrt((3*pi)/2);
end

function val=var(A0)
%%VAR Obtain the variance of the one-dominant-plus-Rayleigh distribution.
%
%INPUTS: A0 The parameter of the distribution. This is the most probable
%           value of the distribution.
%
%OUTPUTS: val The variance of the one-dominant-plus-Rayleigh distribution.
%
%The distribution is defined in Equation 11.4-1 of [1]. The variance was
%found using the identity
%variance=E{x^2}-E{x}^2
%and the necessary expected values were computed by integration.
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% A0=2;
% varVal=OneDominantPlusRayleighD.var(A0)
% numSamp=1e6;
% sampVarVal=var(OneDominantPlusRayleighD.rand([numSamp,1],A0))
%One will find both values are about 0.6209.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    val=(1/24)*A0^2*(32-9*pi);
end

function vals=PDF(x,A0)
%%%PDF Evaluate the one-dominant-plus-Rayleigh probability distribution
%      function (PDF) at one or more desired points.
%
%INPUTS: x The point(s) at which the PDF is to be  evaluated. x>=0 for
%          nonzero values.
%       A0 The parameter of the distribution. This is the most probable
%          value of the distribution.
%
%OUTPUTS: The value(s) of the PDF.
%
%The distribution is defined in Equation 11.4-1 of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% A0=2;
% numSamples=1e5;
% figure(1)
% clf
% histogram(OneDominantPlusRayleighD.rand([numSamples,1],A0),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,6,numPoints);
% vals=OneDominantPlusRayleighD.PDF(x,A0);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    vals=zeros(size(x));
    sel=(x>=0);

    vals(sel)=9*x(sel).^3/(2*A0^4).*exp(-3*x(sel).^2/(2*A0^2));
end
    
function val=CDF(x,A0)
%%CDF Evaluate the cumulative distribution function (CDF) of the one-
%     dominant-plus-Rayleigh distribution at one or more desired points.
%
%INPUTS: x The point(s) at which the CDF is to be evaluated.
%       A0 The parameter of the distribution. This is the most probable
%          value of the distribution.
%
%OUTPUTS: val The value(s) of the CDF at the given points.
%
%The distribution is defined in Equation 11.4-1 of [1]. The CDF was found
%via integration.
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% A0=2;
% x=3;
% numSamples=1e5;
% prob=OneDominantPlusRayleighD.CDF(x,A0)
% probSamp=mean(OneDominantPlusRayleighD.rand([numSamples,1],A0)<=x)
%One will see that both values are about 0.850.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=1-exp(-((3*x.^2)/(2*A0^2))).*(1+(3*x.^2)/(2*A0^2));
end

function x=invCDF(prob,A0)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the one-dominant-plus-Rayleigh distribution.
%    
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          A0 The parameter of the distribution. This is the most probable
%             value of the distribution.
%
%OUTPUTS: x The argument(s) of the CDF that would give the probability or
%           probabilities in prob.
%
%The distribution is defined in Equation 11.4-1 of [1]. The CDF was found
%via integration. The inverse of the CDF involves the Lambert-W function.
%The Lambert W function can have up to 2 real solutions. The -1 solution
%branch was chosen, because it is the one guaranteeitng x to be positive.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% A0=2;
% x=3;
% xBack=OneDominantPlusRayleighD.invCDF(OneDominantPlusRayleighD.CDF(x,A0),A0)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    x=-(2/3)*A0^2*(1+LambW((prob-1)/exp(1),-1));
    x(prob<=0)=0;
    %Deal with possible finite precision errors with the abs command.
    x=sqrt(abs(x));
end

function vals=rand(N,A0)
%%RAND Generate one-dominant-plus-Rayleigh random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of  random variables.
%       A0 The parameter of the distribution. This is the most probable
%          value of the distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated one-dominant-plus-Rayleigh random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
%
%REFERENCES:
%[1] S. H. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end
    
    U=rand(dims);

    vals=OneDominantPlusRayleighD.invCDF(U,A0);
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
