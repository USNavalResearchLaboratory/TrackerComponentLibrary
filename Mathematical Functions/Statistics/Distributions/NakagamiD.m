classdef NakagamiD
%%NAKAGAMID Functions to handle the Nakagami-m distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(m,Omega)
%%MEAN  Obtain the mean of the Nakagami-m distribution for given shape
%       and spread parameters.
%
%INPUTS:    m     The shape parameter of the distribution.
%           Omega The spread parameter of the distribution.
%
%OUTPUTS: val  The mean of the Nakagami-m distribution.
%
%The logarithms are used to avoid overflows in the gamma function ratio 
%when m is large.
%
%The mean of the Nakagami-m distribution is given on the inside of the
%front cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=gammaln(m+1/2)-gammaln(m)+log(Omega/m)/2;

    val=exp(val);
end

function val=var(m,Omega)
%%VAR   Obtain the variance of the Nakagami-m distribution for given
%       shape and spread parameters.
%
%INPUTS:    m     The shape parameter of the distribution.
%           Omega The spread parameter of the distribution.
%
%OUTPUTS: val  The variance of the Nakagami-m distribution.
%
%The logarithms are used to avoid overflows in the gamma function ratio 
%when m is large.
%
%The variance of the Nakagami-m distribution is given on the inside of the
%front cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    gammaRat=exp(gammaln(m+1/2)-gammaln(m));
    
    val=Omega*(1-gammaRat^2/m);
end  
    
function val=PDF(x,m,Omega)
%%PDF          Evaluate the probability density function (PDF) of the
%              Nakagami-m  distribution at one or more desired points.
%
%INPUTS:    x   The point(s) at which the Nakagami-m PDF is to be 
%               evaluated.
%           m   The shape parameter of the distribution.
%           Omega The spread parameter of the distribution.
%
%OUTPUTS:   val The value(s) of the Nakagami-m PDF with shape parameter m
%               and spread parameter Omega evaluated at x.
%
%The PDF of the Nakagami-m distribution is given on the inside of the
%front cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(x>=0);

    val(sel)=(2/gamma(m))*((m/Omega)^m)*x(sel).^(2*m-1).*exp(-(m/Omega)*x(sel).^2);
end 

function prob=CDF(x,m,Omega)
%%CDF       Evaluate the cumulative distribution function (CDF) of the
%           Nakagami-m distribution at one or more desired points.
%
%INPUTS:    x   The point(s) at which the Nakagami-m CDF is to be 
%               evaluated.
%           m   The shape parameter of the distribution.
%           Omega The spread parameter of the distribution.
%
%OUTPUTS:   prob The value(s) of the Nakagami-m CDF with shape parameter m
%                and spread parameter Omega evaluated at x.
%
%The CDF of the Nakagami-m distribution can be expressed in terms of the
%incomplete gamma function as described in [1]. The incomplete gamma
%function is built into Matlab without the use of any toolboxes.
%
%REFERENCES:
%[1] %M. Pätzold, Mobile Radio Channels, Ed. 2, Wiley, 2012, pg. 30.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=zeros(size(x));
    sel=(x>=0);

    prob(sel)=gammainc((m/Omega)*x(sel).^2,m);
end

function val=invCDF(prob,m,Omega)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the Nakagami-m distribution.
%
%INPUTS:    prob The probability or probabilities (0<=prob<=1) at which the 
%                argument of the CDF is desired.
%           m    The shape parameter of the distribution.
%           Omega The spread parameter of the distribution.
%
%OUTPUTS:   val  The argument(s) of the CDF that would give the probability
%                or probabilities in prob.
%
%The CDF of the Nakagami-m distribution can be expressed in terms of the
%incomplete gamma function as described in [1]. The inverse incomplete
%gamma function is built into Matlab without the use of any toolboxes and
%thus is called here.
%
%REFERENCES:
%[1] %M. Pätzold, Mobile Radio Channels, Ed. 2, Wiley, 2012, pg. 30.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=sqrt((Omega/m)*gammaincinv(prob,m));
end

function vals=rand(N,m,Omega)
%%RAND          Generate Nakagami-m distributed random variables with the
%               given parameters.
%
%INPUTS:    N      If N is a scalar, then rand returns an NXN 
%                  matrix of random variables. If N=[M,N1] is a two-element  
%                  row vector, then rand returns an MXN1 matrix of random
%                  variables.
%           m      The shape parameter of the distribution.
%           Omega  The spread parameter of the distribution.
%
%OUTPUTS:   vals   A matrix whose dimensions are determined by N of the
%                  generated Nakagami-m random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    U=rand(dims);

    vals=NakagamiD.invCDF(U,m,Omega);
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
