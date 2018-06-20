classdef InverseGammaD
%%INVERSEGAMMAD Functions to handle the inverse gamma distribution.
%Implemented methods are: mean, var, PDF, CDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(k, beta)
%%MEAN Obtain the mean of the inverse gamma distribution for given shape
%      and scale parameters.
%
%INPUTS: k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
%
%OUTPUTS: val The mean of the inverse gamma distribution.
%
%If k<1, then the distribution has no mean and a NaN is returned.
%
%The distribution arises when estimating radar cross sections in Ch.
%10.3.2 of [1].
%
%REFERENCES:
%[1] W. Koch, Tracking and Sensor Fusion: Methodological Framework and
%    Selected Applications. Berlin, Germany: Springer, 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(k<1)
        val=NaN;
    else
        val=beta/(k-1);
    end
end

function val=var(k,beta)
%%VAR Obtain the variance of the inverse gamma distribution for given shape
%     and scale parameters.
%
%INPUTS: k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
%
%OUTPUTS: val  The variance of the inverse gamma distribution.
%
%If k<2, then the distribution has no variance and a NaN is returned.
%
%The distribution arises when estimating radar cross sections in Ch.
%10.3.2 of [1].
%
%REFERENCES:
%[1] W. Koch, Tracking and Sensor Fusion: Methodological Framework and
%    Selected Applications. Berlin, Germany: Springer, 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(k<2)
        val=NaN;
    else
        val=beta^2/((k-1)^2*(k-2));
    end
end    
    
function val=PDF(x,k,beta)
%%PDF Evaluate the probability density function (PDF) of the inverse gamma
%     distribution at the desired points.
%
%INPUTS: x The point(s) at which the inverse gamma PDF is to be evaluated.
%          Note that x>=0.
%        k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
%
%OUTPUTS: val The value(s) of the PDF of the inverse gamma distribution
%             with parameters k and beta evaluated at x.
%
%The inverse gamma distribution is the distribution of the inverse of a
%gamma random variable. Thus, with a little bit of algebra, it can be
%expressed in terms of the gamma PDF.
%
%The distribution arises when estimating radar cross sections in Ch.
%10.3.2 of [1].
%
%REFERENCES:
%[1] W. Koch, Tracking and Sensor Fusion: Methodological Framework and
%    Selected Applications. Berlin, Germany: Springer, 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(beta^k/gamma(k)).*x.^(-k-1).*exp(-beta./x);
end
    
function val=CDF(x,k,beta)
%%CDF Evaluate the cumulative distribution function (CDF) of the inverse
%     gamma distribution at desired points.
%
%INPUTS: x The point(s) at which the inverse gamma CDF is to be  evaluated.
%          Note that x>=0.
%        k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
%
%OUTPUTS: prob The value(s) of the CDF of the inverse gamma distribution
%              with parameters k and beta evaluated at x.
%
%The inverse gamma distribution is the distribution of the inverse of a
%gamma random variable. Thus, with a little bit of algebra, it can be
%expressed in terms of the gamma PDF.
%
%The distribution arises when estimating radar cross sections in Ch.
%10.3.2 of [1].
%
%REFERENCES:
%[1] W. Koch, Tracking and Sensor Fusion: Methodological Framework and
%    Selected Applications. Berlin, Germany: Springer, 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=1-GammaD.CDF(1./x,k,1/beta);
end

function vals=rand(N,k,beta)
%%RAND Generate inverse gamma distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row  vector, then rand
%          returns an MXN1 matrix of random variables.
%        k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated inverse gamma random variables.
%
%This just takes the inverse of gamma distributed random variables.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    vals=1./randGamma(N,k,1/beta);
end

function entropyVal=entropy(k,beta)
%%ENTROPY Obtain the differential entropy of the inverse gamma
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: k The shape parameter of the distribution; k>0.
%     beta The scale parameter of the distribution; beta>0.
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

    entropyVal=k+log(beta)+gammaln(k)-(1+k)*psi(k);
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
