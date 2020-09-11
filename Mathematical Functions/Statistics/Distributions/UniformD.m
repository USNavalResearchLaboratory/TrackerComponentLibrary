classdef UniformD
%%UNIFORMD Functions to handle the scalar and multivariate (hyper-
%     rectangularly bounded) uniform distributions.
%Implemented methods are: mean, cov, PDF, CDF (scalar only), rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(bounds)
%%MEAN Obtain the mean of the uniform distribution
%
%INPUTS: bounds A 2XnumDim matrix of the minimum and maximum bounds of each
%               of the numDim dimensions of the uniform distribution.
%               bounds(1,:) is the set of minimum bounds and bounds(2,:) is
%               the set of maximum bounds.
%
%OUTPUTS: val  The numDimX1 mean od the uniform distribution.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(bounds(2,:)+bounds(1,:))'/2;
end
  
function val=cov(bounds)
%%COV Obtain the covariance matrix of the uniform distribution (the
%     variance if scalar).
%
%INPUTS: bounds A 2XnumDim matrix of the minimum and maximum bounds of each
%               of the numDim dimensions of the uniform distribution.
%               bounds(1,:) is the set of minimum bounds and bounds(2,:) is
%               the set of maximum bounds.
%
%OUTPUTS: val  The covariance matrix of the uniform distribution. This is a
%              diagonal matrix.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=diag((1/12)*(bounds(2,:)-bounds(1,:)).^2);
end

function vals=PDF(z,bounds)
%%PDF Evaluate a scalar or multivariate uniform  PDF at a certain point
%     given the bounds of the PDF.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluation at
%          multiple points are desired, then this is a matrix with each
%          column being the a point (a vector).
%   bounds A 2XnumDim matrix of the minimum and maximum bounds of each of
%          the numDim dimensions of the uniform distribution. bounds(1,:)
%          is the set of minimum bounds and bounds(2,:) is the set of
%          maximum bounds.
%
%OUTPUTS: vals The scalar value of the uniform PDF with the given bounds.
%              If multiple points are passed (z is a matrix), then val is a
%              column vector.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numDim=size(z,1);
    numPoints=size(z,2);
    val=1/prod(bounds(2,:)-bounds(1,:));
    vals=repmat(val,[numPoints,1]);
    
    sel=false(1,numPoints);
    for curDim=1:numDim
        sel=sel|z(curDim,:)<bounds(1,curDim)|z(curDim,:)>bounds(2,curDim);
    end
    vals(sel)=0;
end

function val=CDF(z,bounds)
%%CDF Evaluate a scalar uniform cumulative distribution function (CDF) at a
%     certain point given the bounds of the distribution.
%
%INPUTS: z The scalar point(s) at which the CDF should be evaluated.
%   bounds A 2X1 vector of the minimum and maximum bounds of the scalar
%          uniform distribution. bounds(1) is the minimum bound and
%          bounds(2) is the maximum bound.
%
%OUTPUTS: val The scalar value(s) of the uniform CDF with the given
%             bounds at point(s) in z.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(z-bounds(1))/(bounds(2)-bounds(1));
    val(z<=bounds(1))=0;
    val(z>=bounds(2))=1;
end

function x=rand(N,bounds)
%%RAND Generate multivariate uniform random variables with a given set of
%      bounds
%
%INPUTS: N The scalar number of random variables to generate.
%  bounds A 2XnumDim matrix of the minimum and maximum integer bounds of
%         each of the numDim dimensions of the uniform distribution.
%         bounds(1,:) is the set of minimum bounds and bounds(2,:) is
%         the set of maximum bounds.
%
%OUTPUT: x A sample of an numDimXN set of random uniform variables with the
%          given bounds.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDim=size(bounds,2);

    x=(bounds(2,:)-bounds(1,:))'.*rand(xDim,N)+bounds(1,:)';
end

function entropyVal=entropy(bounds)
%%ENTROPY Obtain the differential entropy of the multivariate uniform
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: bounds A 2XnumDim matrix of the minimum and maximum integer bounds 
%               of each of the numDim dimensions of the uniform
%               distribution. bounds(1,:) is the set of minimum bounds and
%               bounds(2,:) is the set of maximum bounds.
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

    entropyVal=log(prod(bounds(2,:)-bounds(1,:)));
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
