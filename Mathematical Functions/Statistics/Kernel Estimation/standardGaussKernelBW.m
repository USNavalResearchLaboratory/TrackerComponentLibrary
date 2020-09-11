function [H,mu,SR]=standardGaussKernelBW(param1,method,spherData,N)
%%STANDARDGAUSSKERNELBW Estimate the (univariate or multivariate) kernel
%             bandwidth when approximating a PDF using a set of samples
%             with a Gaussian kernel. This function provides the bandwidth
%             estimate based on the optimal solution if the underlying
%             density of the samples is normal. If the density is not
%             normal, then this function is just an approximation. The
%             kernel bandwidth estimate can be used with the points xi in
%             the function kernelApprox to approximate the continuous
%             density.
%
%INPUTS: param1 If N is an empty matrix or is omitted, then this is:
%              xi A numDimXnumSamples set of samples of the distribution
%                 from which a Gaussian kernel estimator should interpolate
%                 a PDF. For multidimensional estimation, there must be
%                 >=numDim samples. Specifically, the sample covariance
%                 matrix must be positive definite.
%           If N is given, then this is:
%              SR A root covariance matrix such that SR*SR' equals the 
%                 (exact or approximated) covariance matrix of the PDF
%                 being approximated.
%    method An optional parameter specifying the method that should be used
%           to estimate the bandwidth. Possible values are
%           0 Use the method implied by Equation 3.30 and 3.28 in [1],
%             which uses the lesser of the sample standard deviation or a
%             scaled interquartile range. This approach is supposed to work
%             well when the underlying distribution of the samples is not
%             actually normal and might be bimodal. This method cannot be
%             used if N is given.
%           1 (The default if omitted or an empty matrix is passed) Use the
%             method of Equation 3.28 in [1], based on the sample standard
%             deviation.
%           2 Use the method of Equation 3.29 in [1], based on the
%             interquartile range of the samples. This method cannot be
%             used if N is given.
% spherData When the underlying density is truly normal, then the
%           approximation here is only valid if the cross terms of the
%           correlation are zero. To make that the case, if this parameter
%           is true, then the data is "sphered" --premultiplied by the
%           inverse of the square root (lower-triangular Cholesky
%           decomposition) of the sample covariance matrix, and then the
%           final result is multipled by the square root of the sample
%           covariance matrix to scale everything back. However, when this
%           function is used for approximations, sphering can sometimes
%           make thing worse. Thus, if this is false, then only the
%           diagonal terms of the sample covariance matrix are used and the
%           data is not squared. The default if this parameter is omitted
%           or an empty matrix is passed is true.
%         N If this parameter is provided, then param1 is assumed to be SR
%           and N is the number of samples that are being smoothed.
%           Otherwise, if this is omitted or an empty matrix is passed,
%           param1 is xi.
%           
%OUTPUTS: H The kernel bandwidth matrix that can be used in the function
%           kernelApprox. This will be a lower-triangular matrix.
%         mu If xi is given, then this is the sample mean. Otherwise, this
%            is an empty matrix.
%         SR The root covariance matrix from the input or computed from xi.
%
%The methods in Chapter 3 of [1] are for a scalar distribution. We are able
%to apply them to a multivariate distribution by pre-multiplying the
%multivariate distribution by the inverse of the sample standard deviation
%of the points. If the samples were from a normal distribution (the
%approximation used here), then we can replace the scalar bandwidth
%formula with the solution for numDim-dimensions given in Chapter 2.4.2 of
%[2]. The pre-multiplication to decorrelate the dimensions is necessary for
%the use of that formula. Ultimately, we then have then rotate the final
%result back. The techniques from Chapter 3 of 1 that are used here only
%differ in how they try to obtain a robust estimate of the standard
%deviation values.
%
%REFERENCES:
%[1] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%[2] A. D. Bowman and A. Azzalini, Applied Smoothing Techniques for Data
%    Analysis: The Kernel Approach with S-Plus Illustrations. Oxford:
%    Clarendon Press, 2004.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    N=[]; 
end

if(nargin<3||isempty(spherData))
    spherData=true;
end

if(nargin<2||isempty(method))
    method=1;
end
numDims=size(param1,1);

if(isempty(N))%xi is given
    xi=param1;
    [mu,S]=calcMixtureMoments(xi);
    SR=cholSemiDef(S,'lower',1);%Not fully triangular. but S=SR*SR'
    if(method~=1&&spherData)
        xi=SR\xi;
    end
    N=size(xi,2);%Number of points
else%S and N are given.
    if(method~=1)
        error('If N is given, then method=1 must be used.')
    end
    
    mu=[];
    SR=param1;
end

if(numDims>1&&spherData)
    sigma1=ones(numDims,1);
else
    sigma1=diag(SR);
end

switch(method)
    case 0%Use Equation 3.30 and 3.28 in [1].
        sigma2=interquartileRange(xi,2)/(GaussianD.invCDF(0.75)-GaussianD.invCDF(0.25));
        sigma=min([sigma1,sigma2],[],2);
    case 1%Use Equation 3.28 in [1].
        sigma=sigma1;
    case 2%Use the interquartile range (scaled appropriately) instead of
          %the standard deviation in [1]. This is Equation 3.29.
        sigma=interquartileRange(xi,2)/(GaussianD.invCDF(0.75)-GaussianD.invCDF(0.25));
    otherwise
        error('Unknown method specified')
end

h=(4/(N*(numDims+2)))^(1/(numDims+4))*sigma;

%Undo any rotation/scaling.
if(numDims>1&&spherData)
    %H is actually aking to the square root of a covariance matrix. h is a
    %list of scalars--the diagonls of the rotated matrix. This builds that
    %diagonal matrix and then scales it.
    H=SR*diag(h);
else
    H=diag(h);
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
