function H=standardGaussKernelBW(xi,method,spherData)
%%STANDARDGAUSSKERNELBW Estimate the (univariate of multivariate) kernel
%             bandwidth when approximating a PDF using a set of samples
%             with a Gaussian kernel. This function provides the bandwidth
%             estimate based on the optimal solution if the underlying
%             density being samples is normal. If the density is not
%             normal, then this function is just an approximation. The
%             kernel bandwidth estimate can be used with the points xi in
%             the function kernelApprox to approximate the continuous
%             density.
%
%INPUTS: xi A numDimXnumSamples set of samples of the distribution from
%           which a Gaussian kernel estimator should interpolate a PDF. For
%           multidimensional estimation, there must be >=numDim samples.
%           Specifically, the sample covariance matrix must be positive
%           definite.
%    method An optional parameter specifying the method that should be used
%           to estimate the bandwidth. Possible values are
%           0 (The default if omitted or an empty matrix is passed) Use the
%             method implied by Equation 3.30 and 3.28 in [1], which uses
%             the lesser of the sample standard deviation or a scaled
%             interquartile range. This approach is supposed to work well
%             when the underlying distribution of the samples is not
%             actually normal and might be bimodal.
%           1 Use the method of Equation 3.28 in [1], based on the sample
%             standard deviation.
%           2 Use the method of Equation 3.29 in [1], based on the
%             interquartile range of the samples.
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
%           
%OUTPUTS: H The kernel bandwidth matrix that can be used in the function
%           kernelApprox. This will be a lower-triangular matrix.
%
%The methods in Chapter 3 of [1] are for a scalar distribution. We are able
%to apply them to a multivariate distribution by pre-multiplying the
%multivaraite distribution by the inverse of the sample standard deviation
%of the points. If the samples were from a normal distribution (the
%approximation used here), then the distributions for all of the dimensions
%would have been made independent by the rotation and we can thus use the
%scalar method on each dimension and then rotate the final result back.
%
%REFERENCES:
%[1] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(spherData))
    spherData=true;
end

if(nargin<2||isempty(method))
    method=0;
end
numDims=size(xi,1);

%Make the components independent if a multivariate distribution is given.
if(numDims>1)
    [~,S]=calcMixtureMoments(xi);
    
    if(spherData)
        %Normalize the estimates.
        SR=chol(S,'lower');
        xi=SR\xi;
        sigma1=ones(numDims,1);
    else
        sigma1=sqrt(diag(S));
    end
else
    sigma1=std(xi,[],2);
end

n=size(xi,2);%Number of points

switch(method)
    case 0%Use Equation 3.30 and 3.28 in [1].
        sigma2=interquartileRange(xi,[],2)/(GaussianD.invCDF(0.75)-GaussianD.invCDF(0.25));
        sigma=min([sigma1,sigma2],[],2);
        h=(4/3)^(1/5)*sigma*n^(-1/5);
    case 1%Use Equation 3.28 in [1].
        h=(4/3)^(1/5)*sigma1*n^(-1/5);
    case 2%Use the interquartile range (scaled appropriately) instead of
          %the standard deviation in [1]. This is Equation 3.29.
        sigma=interquartileRange(xi,[],2)/(GaussianD.invCDF(0.75)-GaussianD.invCDF(0.25));
        h=(4/3)^(1/5)*sigma*n^(-1/5);
    otherwise
        error('Unknown method specified')
end

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
