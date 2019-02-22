classdef ComplexGaussianD
%%COMPLEXGAUSSIAND Functions to hand scalar and multivariate circularly-
%                  symmetric complex Gaussian distributions.
%Implemented methods are: mean,cov, PDF, PDFI, PDFS, rand, randS
%
%When dealing with circularly-symmetric complex normal distributions, one
%must pay attention to factors of 2 and sqrt(2) that appear here and then
%and are not in real Gaussian distributions. For example, if z is a
%zero-mean scalar circularly-symmetric complex normally distributed random
%variable with variance sigma^2, then norm(z) is Rayleigh distributed
%random variable with parameter sigma/sqrt(2). This constrasts with
%sqrt(x1^2+x2^2), where x1 and x2 are independent zero mean real normal
%random variables with standard deviations sigma^2. In that instance,
%norm(z) is Rayleigh distributed random variable with parameter
%sigma. 
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(mu)
%%MEAN Obtain the mean of the circularly-symmetric complex Gaussian
%      distribution.
%
%INPUTS: mu The mean of the PDF. If the PDF is multivariate, then this ia a
%           column vector.
%
%OUTPUTS: val The mean of the Gaussian distribution.
%
%The circularly-symmetric complex Gaussian distribution is parameterized by
%its mean and covariance matrix. Thus, this function just returns the mean
%it is given.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end

function val=cov(Sigma)
%%COV Obtain the covariance matrix of the Gaussian distribution (the
%     variance if scalar).
%
%INPUTS: Sigma The variance (if scalar) or Hermitian covariance matrix (if
%              multidimensional) of the PDF. The variance cannot be zero
%              and the covariance matrix cannot be singular.
%
%OUTPUTS: val The covariance matrix of the Gaussian distribution.
%
%The circularly-symmetric complex Gaussian distribution is parameterized by
%its mean and covariance matrix. Thus, this function just returns the
%covariance matrix
%it is given.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=Sigma;
end

function vals=PDF(z,mu,Sigma)
%%PDF Evaluate a scalar or multivariate circularly-symmetric complex
%     Gaussian (normal) PDF at a certain point given the mean and the
%     covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If values at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDImX1 column vector. If omitted or an empty matrix is passed,
%          a zero mean is used.
%    Sigma The variance (if scalar) or numDimXnumDim Hermitian covariance
%          matrix (if multidimensional) of the PDF. The variance cannot be
%          zero and the covariance matrix cannot be singular. If omitted or
%          an empty matrix is passed, the identity matrix is used as the
%          covariance matrix.
%
%OUTPUTS: vals The scalar values of the complex normal PDF with mean mu
%              and covariance matrix Sigma evaluated at the points in z. If
%              multiple points are passed (z is a matrix), then val is a
%              row vector.
%
%The circularly-symmetric complex normal distribution is presented with
%respect to the real normal distribution in [1].
%
%REFERENCES:
%[1] N. R. Goodman, "Statistical analysis based on a certain multivariate
%    complex Gaussian distribution (an introduction)," Annals of
%    Mathematical Statistics, vol. 34, no. 1, pp. 152-177, Mar. 1963.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numDim=size(z,1);
    numPoints=size(z,2);
    vals=zeros(1,numPoints);
    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(Sigma))
       Sigma=eye(numDim,numDim); 
    end
    
    
    for curPoint=1:numPoints
        diffVal=z(:,curPoint)-mu;
        %The determinant of a Hermitian matrix is real. Also, the quadratic
        %argument of the exponent function is real. The real commands just
        %get rid of issues due to finite precision errors.
        vals(curPoint)=(real(det(pi*Sigma)))^(-1)*exp(-real(diffVal'*inv(Sigma)*diffVal));
    end
end

function vals=PDFI(z,mu,SigmaInv)
%%PDFI Evaluate a scalar or multivariate circularly-symmetric complex
%      Gaussian (normal) PDF at a certain points given the mean and the
%      inverse of the covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If values at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDimX1 column vector. If omitted or an empty matrix is passed,
%          a zero mean is used.
% SigmaInv The inverse variance (if scalar) or numDimXnumDim inverse
%          covariance matrix (if multidimensional) of the PDF. SigmaInv can
%          be singular. If omitted or an empty matrix is passed, the
%          identity matrix is used as the covariance matrix.
%
%OUTPUTS: val The scalar value of the normal PDF with mean mu and inverse
%             covariance matrix SigmaInv evaluated at the point z. If
%             multiple points are passed (z is a matrix), then val is a row
%             vector.
%
%The circularly-symmetric complex normal distribution is presented with
%respect to the real normal distribution in [1].
%
%REFERENCES:
%[1] N. R. Goodman, "Statistical analysis based on a certain multivariate
%    complex Gaussian distribution (an introduction)," Annals of
%    Mathematical Statistics, vol. 34, no. 1, pp. 152-177, Mar. 1963.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numPoints=size(z,2);
    vals=zeros(1,numPoints);
    n=size(z,1);
    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(SigmaInv))
       SigmaInv=eye(numDim,numDim); 
    end
    
    for curPoint=1:numPoints
        %Note that det(A^(-1))=1/det(A) and that det(a*A)=a^n*det(A), where
        %a is a scalar and A is an nXn matrix.
        
        diffVal=z(:,curPoint)-mu;
        %The determinant of a Hermitian matrix is real. Also, the quadratic
        %argument of the exponent function is real. The real commands just
        %get rid of issues due to finite precision errors.
        vals(curPoint)=(pi)^(-n)*real(det(SigmaInv))*exp(-real(diffVal'*SigmaInv*diffVal));
    end
end

function val=PDFS(z,mu,S)
%%GAUSSIANPDFS Evaluate a scalar or multivariate circularly-symmetric
%              complex Gaussian (normal) PDF at a certain point given the
%              mean and the square root of the covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If values at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDimX1 column vector.
%        S The square root of the variance (if scalar) or the numDimXnumDim
%          lower-triangular square root of the covariance matrix (if
%          multidimensional) of the PDF such that S*S'=Sigma, where Sigma
%          is the covariance matrix. S cannot be a singular matrix. If
%          omitted or an empty matrix is passed, the identity matrix is
%          used.
%
%OUTPUTS: val The scalar value(s) of the normal PDF with mean mu and square
%             root covariance matrix S evaluated at the points in z. I
%             multiple points are passed (z is a matrix), then val is a
%             row vector.
%
%The circularly-symmetric complex normal distribution is presented with
%respect to the real normal distribution in [1].
%
%REFERENCES:
%[1] N. R. Goodman, "Statistical analysis based on a certain multivariate
%    complex Gaussian distribution (an introduction)," Annals of
%    Mathematical Statistics, vol. 34, no. 1, pp. 152-177, Mar. 1963.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(S))
       S=eye(numDim,numDim); 
    end

%Note that (S*S')^(-1)=(S')^(-1)*S^(-1)
    diff=S\bsxfun(@minus,z,mu);
%Note that det(S*S')=det(S)*det(S') and that det(S)=det(S') so
%det(S*S')=det(S)^2. Also, det(a*S)=a^ndet(S), where a is a scalar and S is
%an nXn matrix. Thus,
%det(2*pi*S*S')=det(sqrt(2*pi)*S)^2=(2*pi)^n*det(S)^2
    n=size(z,1);
    %The abs in the determinant is necessary if the main diagonal of S has
    %negative terms. S can still be such that S*S'=Sigma as the sign of
    %those terms is not unique due to the squaring. Also, the quadratic
    %argument of the exponent function is real. The real commands just
    %get rid of issues due to finite precision errors.
    
    val = (1/((pi)^(n)*abs(det(S)^2)))*exp(-real(sum(diff.*conj(diff),1))); 
end


function x=rand(N,mu,P)
%%RAND Generate multivariate circularly-symmetric complex Gaussian random
%      variables with a given mean vector and covariance matrix.
%
%INPUTS: N The number of random variables to generate.
%       mu The xDim X1 mean of the multivariate Gaussian to generate. If
%          omitted or an empty matrix is passed, a mean of zero is used.
%          and xDim is taken to be 1.
%        P The xDim X xDim positive definite hermitian covariance matrix
%          of the multivariate Gaussian to generate. If this parameter is
%          omitted or an empty matrix is passed, then the identity matrix
%          will be used.
%
%OUTPUT: x An xDimXN matrix of random instances of the multivariate
%          Gaussian distribution.
%
%The circularly-symmetric complex normal distribution is presented with
%respect to the real normal distribution in [1]. The variables are
%generated as the sum of two real normal distributed random vectors, with
%with covariance matrices P/2, one multiplied by sqrt(-1).
%
%REFERENCES:
%[1] N. R. Goodman, "Statistical analysis based on a certain multivariate
%    complex Gaussian distribution (an introduction)," Annals of
%    Mathematical Statistics, vol. 34, no. 1, pp. 152-177, Mar. 1963.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(mu))
        xDim=1;
        mu=zeros(xDim,1);
    else
        xDim=size(mu,1);
    end

    if(nargin<3||isempty(P))
        P=eye(xDim,xDim);
    end

    xDim=size(mu,1);
    
    S=chol(P,'lower');
    x=bsxfun(@plus,mu,(1/sqrt(2))*(S*randn(xDim,N)+1j*S*randn(xDim,N)));
end


function x=randS(N,mu,S)
%%RANDS Generate multivariate Gaussian random variable with a given mean
%       vector and lower-triangular square root covariance matrix.
%
%INPUTS: N The number of random variables to generate.
%       mu The xDim X1 mean of the multivariate Gaussian to generate.
%        S The xDim X xDim lower triangular square root covariance matrix
%          of the multivariate Gaussian to generate. If this parameter is
%          omitted or an empty matrix is passed, then the identity matrix
%          is used.
%
%OUTPUT: x An xDimXN matrix of random instances of the multivariate
%          Gaussian distribution.
%
%The circularly-symmetric complex normal distribution is presented with
%respect to the real normal distribution in [1]. The variables are
%generated as the sum of two real normal distributed random vectors, with
%with covariance matrices P/2, one multiplied by sqrt(-1).
%
%REFERENCES:
%[1] N. R. Goodman, "Statistical analysis based on a certain multivariate
%    complex Gaussian distribution (an introduction)," Annals of
%    Mathematical Statistics, vol. 34, no. 1, pp. 152-177, Mar. 1963.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDim=size(mu,1);
    
    if(nargin<3||isempty(S))
        S=eye(xDim,xDim);
    end
    
    x=bsxfun(@plus,mu,(1/sqrt(2))*(S*randn(xDim,N)+1j*S*randn(xDim,N)));
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

