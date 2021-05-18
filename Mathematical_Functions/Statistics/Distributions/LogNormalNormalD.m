classdef LogNormalNormalD
%%LOGNORMALNORMALD Functions to handle a hybrid normal and lognormal
%   distribution. Given jointly normal random variables with mean mu and
%   covariance matrix Sigma, a hybrid normal lognormal random variable y is
%   one where y(1:NLog)=exp(x(1:NLog) and the rest of the element of y are
%   just the corresponding normal elements of x. Such a joint normal-
%   lognormal distribution might be useful, for example, for range
%   measurements that are modeled as log-normal and Doppler measurments
%   that are modeled as normal.
%Implemented methods are: PDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function vals=PDF(z,mu,Sigma,NLog)
%%PDF          Evaluate the PDF of a joint normal-lognormal random vector
%              at certain points given the mean and covariance matrix of
%              the underlying Gaussian noise.
%
%INPUTS:    z   The point at which the PDF should be evaluated. If the PDF
%               is multivariate, then this is a column vector. If
%               evaluation at multiple points are desired, then this is a
%               matrix with each column being the a point (a vector). This
%               must be non-negative. The first NLog components of each
%               vector must be the lognormal ones.
%           mu  The mean of the normal distribution that drives the normal
%               and lognormal components of the noise. This is an xDimX1
%               column vector.
%         Sigma The covariance matrix of the Gaussian PDF that is
%               transformed to drive the joint lognormal-normal
%               distribution. This cannot be singular.
%          NLog The number of components of the vectors in z that are the
%               lognormal ones. The first NLog vectors are lognormal.
%
%OUTPUTS: vals  The scalar value of the joint lognormal-normal PDF
%               evaluated at z. If z is a matrix, then vals will be a
%               row vector.
%
%The implementations of the PDF of this class is taken from [1].
%
%REFERENCES:
%[1] S. J. Fletcher and M. Zupanski, "A hybrid multivariate normal and
%    lognormal distribution for data assimilation," Atmospheric Science
%    Letters, vol. 7, no. 2, pp. 43-46, Apr./Jun. 2006.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numPoints=size(z,2);
    normFac=1/sqrt(det(Sigma*2*pi));
    vals=zeros(1,numPoints);
    
    for curPoint=1:numPoints
        prodVal=prod(1./z(1:NLog,curPoint));
        zTilde=[log(z(1:NLog));z((NLog+1):end)];
        diff=zTilde-mu;
        vals(curPoint)=normFac*prodVal*exp(-0.5*invSymQuadForm(diff,Sigma));
    end 
end


function x=rand(N,mu,Sigma,NLog)
%%RAND   Generate hybird lognormal-normal random variables with a given
%        mean vector and covariance matrix.
%
%INPUTS: mu     The mean of the normal distribution that drives both the
%               lognormal component as well as the normal component.
%        Sigma  The variance or covariance matrix of the Gaussian PDF that
%               drives the lognromal and normal components.
%         NLog  The number of elements of the random vector that are
%               distributed log-normal.
%
%OUTPUT: x   An xDimXN matrix of random instances of the lognormal-normal
%            distribution. The first NLog components of each random
%            instance are the lognormal components.
%
%The lognormal distribution part of the vector just comes from raising e to
%a normal random variable.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    x=GaussianD.rand(N,mu,Sigma);
    x(1:NLog)=exp(x(1:NLog));
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
