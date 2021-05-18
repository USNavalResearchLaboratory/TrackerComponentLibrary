function kernelEst=kernelApprox(x,xi,H,K)
%%KERNELAPPROX Perform multivariate kernel density estimation given a set
%           of random uniformely weighted probability distribution function
%           (PDF) samples, a bandwidth estimate and a kernel. One wants to
%           interpolate values of the PDF  between the points. Thus, one
%           puts some type of a kernel (typically a Gaussian PDF) around
%           each point. However, the scaling of the kernel (the covariance
%           matrix in a Gaussian PDF), determines how peaky the results
%           will be. A function like getKernelBandwidthMatrix can estimate
%           a good bandwidth matrix to use. This function allows one to use
%           the points, kernel, and bandwidth matrix to interpolate values
%           between the points.
%
%INPUTS: x A numDimXnumPoints matrix of numPoints points where interpolated
%          PDFvalues are desired.
%       xi A numDimXnumNodes set of samples of the PDF.
%        H The numDimXnumDim square bandwidth of the kernel estimator.
%          This is often a lower-triangular matrix. See below for details.
%        K This selects the kernel to use. Possible inputs are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Gaussian kernel.
%          1 Use the Epanechnikov kernel.
%          One can also pass a function handle to the kernel. A kernel
%          should be some type of PDF. K(x) should be able to return a
%          1XnumNodes set of values when given a numDimXnumNodes set of x
%          values.
%
%OUTPUTS: kernelEst A numDimXnumPoints matrix of the estimated value of the
%               PDF implied by the samples of the PDF and the kernel at
%               each of the points in x.
%
%As described in Chapter 6 of [1], a method of interpolating a PDF given
%points is to use the formula
%f=(1/n)*sum_{i=1}^n(1/det(H))*K(H\(x-xi(:,i))
%That is what this function implements.
%
%REFERENCES:
%[1] D. W. Scott, Multivariate Density Estimation: Theory, Practice, and
%    Visualization, 2nd ed. Hoboken, NJ: Wiley, 2015.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(x,2);

n=size(xi,2);

if(nargin<4||isempty(K))
    %Use a standard Gaussian kernel if none is specified.
    K=0;
end

if(~isa(K,'function_handle'))
    switch(K)
        case 0%Gaussian kernel 
            K=@(x)GaussianKernel(x);
        case 1%Epanechnikov kernel
            K=@(x)EpanechnikovKernel(x);
        otherwise
            error('Unknown kernel specified.')
    end
end

sumVal=zeros(numPoints,1);
for curPoint=1:numPoints
    sumVal(curPoint)=sum(K(H\bsxfun(@minus,x(:,curPoint),xi)));
end

kernelEst=(1/n)*(1/abs(det(H)))*sumVal;
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
