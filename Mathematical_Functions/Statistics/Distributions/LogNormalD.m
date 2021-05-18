classdef LogNormalD
%%LOGNORMALD Functions to handle the scalar and multivariate lognormal
%   distribution. A lognormal random variable is obtained by evaluating
%   e^x, where x is a normal random variable. If x is a vector, then the
%   exponent is performed element-by-element.
%Implemented methods are: mean,cov, PDF, CDF, rand, entropy (only for
%                                                    scalar distributions)
%
%An expression for the PDF of a multivariate lognormal distribution is
%given in [1], where it is shown that a joint multivariate log-normal and
%normal PDF can be formed.
%
%Note that as derived in [2], the lognormal distribution cannot be uniquely
%identified via its moments.
%
%REFERENCES:
%[1] S. J. Fletcher and M. Zupanski, "A hybrid multivariate normal and
%    lognormal distribution for data assimilation," Atmospheric Science
%    Letters, vol. 7, no. 2, pp. 43-46, Apr./Jun. 2006.
%[2] C. C. Heyde, "On a property of the lognormal distribution," The
%    Journal of the Royal Statistical Society Series B (Methodological),
%    vol. 25, no. 2, pp. 392-393, 1963.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(mu,Sigma)
%%MEAN Obtain the mean of the lognormal distribution.
%
%INPUTS: mu The mean of the normal distribution that is transformed to
%           obtain the lognormal distribution. If the PDF is multivariate,
%           then this is a column vector.
%     Sigma The variance (if scalar) or covariance matrix (if
%           multidimensional) of the Gaussian PDF that is transformed to
%           get the lognormal distribution. The variance cannot be zero
%           and the covariance matrix can not be singular.
%
%OUTPUTS: val The mean of the lognormal distribution.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=exp(mu+0.5*diag(Sigma));

end

function val=cov(mu,Sigma)
%%VAR Obtain the covariance matrix of the lognormal distribution (the
%     variance if scalar).
%
%INPUTS: mu The mean of the normal distribution that is transformed to
%           obtain the lognormal distribution. If the PDF is multivariate,
%           then this is a column vector.
%     Sigma The variance (if scalar) or covariance matrix (if
%           multidimensional) of the Gaussian PDF that is transformed to
%           get the lognormal distribution. The variance cannot be zero and
%           the covariance matrix can not be singular.
%
%OUTPUTS: val The variance of the lognormal distribution.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    muDim=size(mu,1);

    val=zeros(muDim,muDim);
    
    for row=1:muDim
        for col=1:muDim
            val(row,col)=exp(mu(row)+mu(col)+0.5*(Sigma(row,row)+Sigma(col,col)));
        end
    end
    
    val=val.*(exp(Sigma)-1);
end

function vals=PDF(z,mu,Sigma)
%%PDF Evaluate a scalar or multivariate lognormal PDF at a certain point
%     given the mean and the covariance matrix.
%
%INPUTS: z The point at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluation at
%          multiple points are desired, then this is a matrix with each
%          column being the a point (a vector). This must be non-negative.
%       mu The mean of the normal distribution that is transformed to
%          obtain the lognormal distribution. If the PDF is multivariate,
%          then this is a column vector.
%    Sigma The variance (if scalar) or covariance matrix (if
%          multidimensional) of the Gaussian PDF that is transformed to get
%          the lognormal distribution. The variance cannot be zero and the
%          covariance matrix can not be singular.
%
%OUTPUTS: vals The scalar value of the lognormal PDF evaluated at z. If z
%              is a matrix, then vals will be a column vector.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

numPoints=size(z,2);
vals=zeros(numPoints,1);
normFac=1/sqrt(det(Sigma*2*pi));

    for curPoint=1:numPoints
        diff=log(z(:,curPoint))-mu;
        vals(curPoint) = normFac*prod(1./z(:,curPoint))*exp(-0.5*invSymQuadForm(diff,Sigma));
    end
    
    %Deal with values where z=0 is passed.
    sel=~isfinite(vals);
    vals(sel)=0;
end

function vals=CDF(z,mu,var)
%%CDF Evaluate a scalar lognormal CDF at a certain point given the mean and
%     the variance of normal distribution transformed to get the lognormal
%     distribution. 
%
%INPUTS: z The point(s) at which the CDF should be evaluated.
%       mu The mean of the scalar normal distribution that is transformed
%          to obtain the lognormal distribution.
%      var The variance of the Gaussian PDF that is transformed to get the
%          lognormal distribution.
%
%OUTPUTS: vals The scalar value(s) of the lognormal CDF with mean mu and
%              variance var evaluated at the point z.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numPoints=length(z);

    vals=zeros(numPoints,1);
    for curPoint=1:numPoints
       vals(curPoint)=(1+erf((log(z)-mu)/sqrt(2*var)))/2;
    end
end

function x=rand(N,mu,Sigma)
%%RAND Generate multivariate lognormal random variables with a given mean
%      vector and covariance matrix.
%
%INPUTS: N The number of random variables to generate.
%       mu The mean of the normal distribution that is transformed to
%          obtain the lognormal distribution.
%    Sigma The variance or covariance matrix of the Gaussian PDF that is
%          transformed to get the lognormal distribution.
%
%OUTPUT: x An xDimXN matrix of random instances of the multivariate
%          lognormal distribution.
%
%The lognormal distribution just comes from raising e to a normal random
%variable.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    x=exp(GaussianD.rand(N,mu,Sigma));
end

function entropyVal=entropy(mu, Sigma)
%%ENTROPY Obtain the differential entropy of the scalar log-normal
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: mu The mean of the normal distribution that is transformed to
%           obtain the lognormal distribution.
%     Sigma The variance of the Gaussian PDF that is transformed to get the
%           lognormal distribution.
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

    entropyVal=mu+(1/2)*log(2*pi*exp(1)*Sigma);
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
