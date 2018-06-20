classdef UniformDiscreteD
%%UNIFORMDISCRETED Functions to handle the scalar and multivariate (hyper-
%     rectangularly bounded) discrete uniform distributions.
%Implemented methods are: mean, cov, PMF, CDF (scalar only), rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static) 
    
function val=mean(bounds)
%%MEAN Obtain the mean of the discrete uniform distribution
%
%INPUTS: bounds A 2XnumDim matrix of the integer minimum and maximum
%               bounds of each of the numDim dimensions of the discrete
%               uniform distribution. bounds(1,:) is the set of minimum
%               bounds and bounds(2,:) is the set of maximum bounds.
%
%OUTPUTS: val  The numDimX1 mean of the discrete uniform distribution.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(bounds(2,:)+bounds(1,:))'/2;
end
    
function val=cov(bounds)
%%COV  Obtain the covariance matrix of the discrete uniform distribution
%      (the variance if scalar).
%
%INPUTS: bounds A 2XnumDim matrix of the integer minimum and maximum
%               bounds of each of the numDim dimensions of the discrete
%               uniform distribution. bounds(1,:) is the set of minimum
%               bounds and bounds(2,:) is the set of maximum bounds.
%
%OUTPUTS: val  The covariance matrix of the discrete uniform distribution.
%              This is a diagonal matrix.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=diag((1/12)*((bounds(2,:)-bounds(1,:)+1).^2-1));
end

function vals=PMF(z,bounds)
%%PMF Evaluate a scalar or multivariate discrete uniform probability mass
%     function (PMF) at a certain point given the bounds of the PMF.
%
%INPUTS: z The integer points at which the PMF should be evaluated. If the
%          PMF is multivariate, then this is a column vector. If evaluation
%          at multiple points are desired, then this is a matrix with each
%          column being the a point (a vector).
%   bounds A 2XnumDim matrix of the integer minimum and maximum bounds of
%          each of the numDim dimensions of the discrete uniform
%          distribution. bounds(1,:) is the set of minimum bounds and
%          bounds(2,:) is the set of maximum bounds.
%
%OUTPUTS: vals The scalar value of the discrete uniform PMF with the given
%              bounds. If multiple points are passed (z is a matrix), then
%              val is a column vector.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numDim=size(z,1);
    numPoints=size(z,2);
    val=1/prod(bounds(2,:)-bounds(1,:)+1);
    vals=repmat(val,[numPoints,1]);
    
    sel=false(1,numPoints);
    for curDim=1:numDim
        sel=sel|z(curDim,:)<bounds(1,curDim)|z(curDim,:)>bounds(2,curDim);
    end
    vals(sel)=0;
end

function val=CDF(z,bounds)
%%CDF Evaluate a scalar discrete uniform cumulative distribution function
%     (CDF) at a certain point given the bounds of the distribution.
%
%INPUTS: z The scalar point(s) at which the CDF should be evaluated.
%   bounds A 2X1 vector of the integer minimum and maximum bounds of the
%          discrete uniform distribution. bounds(1) is the minimum bound
%          and bounds(2) is the maximum bound.
%
%OUTPUTS: val The scalar value(s) of the discrete uniform CDF with the
%             given bounds at point(s) in z.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(fix(z)-bounds(1)+1)/(bounds(2)-bounds(1)+1);
    val(z<bounds(1))=0;
    val(z>bounds(2))=1;
end

function x=rand(N,bounds)
%%RAND Generate multivariate discrete uniform random variables with a given
%      set of integer bounds
%
%INPUTS: N The scalar number of random variables to generate
%   bounds A 2XnumDim matrix of the integer minimum and maximum bounds of
%          each of the numDim dimensions of the discrete uniform
%          distribution. bounds(1,:) is the set of minimum bounds and
%          bounds(2,:) is the set of maximum bounds.
%
%OUTPUT: xA sample of an numDimXN set of random discrete uniform variables
%           with the given bounds.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    xDim=size(bounds,2);
    %The max of eps(0) makes sure that if rand returns 0, we don't get a
    %number below the valid range. Thus, the minimum of the ceil function
    %is 1.
    x=ceil((bounds(2,:)-bounds(1,:))'.*max(rand(xDim,N),eps(0)))+bounds(1,:)'-1;
end

function entropyVal=entropy(bounds)
%%ENTROPY Obtain the Shannon entropy of the discrete uniform distribution
%         given in nats. The Shannon entropy of a discrete distribution is
%         entropy=-sum_x Pr(x)*log(Pr(x)) where the sum is over all
%         discrete values of x. Units of nats mean that the natural
%         logarithm is used in the definition.
%
%INPUTS: bounds A 2XnumDim matrix of the integer minimum and maximum bounds
%               of each of the numDim dimensions of the discrete uniform
%               distribution. bounds(1,:) is the set of minimum bounds and
%               bounds(2,:) is the set of maximum bounds.
%
%OUTPUTS: entropyVal The value of the Shannon entropy in nats.
%
%Shannon entropy is defined in Chapter 2 of [1].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    entropyVal=log(prod(bounds(2,:)-bounds(1,:)+1));
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
