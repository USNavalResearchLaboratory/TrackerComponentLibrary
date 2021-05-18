classdef LaplaceD
%%LAPLACED Functions to handle the multivariate Laplace distribution. In
%       one dimension, it reduces to the traditional Laplace distribution.
%Implemented methods are: mean, cov, PDF, CDF (for scalar distributions),
%                         momentGenFun, cumGenFun, rand, entropy (for
%                         scalar distributions)
%
%The multivariate (and univariate) Laplce distributions are described in
%detail in [1]. They arise from certain transformations of Gaussian random
%variables. Specifically, they arise from
%y=mu+sqrt(z)*chol(Gamma,'lower')*x
%x is a normal 0-I random  vector, z is an exponential random variable with
%rate parameter lambda, mu is a vector, and Gamma is a positive definite
%matrix with det(Gamma)=1.
%
%Simple multivariate laplace distributions often arise when considering
%compressive sensing problems as in [2].
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%[2] S. D. Babacan, R. Molina, and A. K. Katsaggelos, "Bayesian compressive
%    sensing using Laplace priors," IEEE Transactions on Image Processing,
%    vol. 19, no. 1, pp. 53-63, Jan. 2010.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(mu)
%%MEAN  Obtain the mean of the multivariate Laplace distribution.
%
%INPUTS: mu The location parameter of the PDF. If the PDF is multivariate,
%            then this is a column vector.
%
%OUTPUTS: val The mean of the multivariate Laplace distribution.
%  
%The mean is given by Equation 10 of [1]. The distribution is aprameterized
%by its mean, so this function just returns the mean it is given.
%
%REFERENCES:  
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp. 300-
%    303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end

function val=cov(lambda,Gamma)
%%COV Obtain the covariance matrix of the multivariate Laplace distribution
%     (the variance if scalar).
%
%INPUTS: lambda A scale parameter for the distribution. Note that with
%               respect to the traditional parameterization of a scalar
%               distribution, sqrt(2/lambda)=1/b, where b is the
%               traditional scale parameter. Here, lambda stems from the
%               rate parameter of an exponential distribution when deriving
%               the multivariate Laplace.
%         Gamma The scale matrix of the PDF. This must be symmetric and
%               positive definite with determinant det(Gamma)=1. When
%               dealing with scalar distributions, this is just 1.
%
%OUTPUTS: val The covariance matrix of the multivariate Laplace
%             distribution.
%
%The covariance is given by Equation 11 of [1].
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=lambda*Gamma;
end


function vals=PDF(x,lambda,mu,Gamma)
%%PDF Evaluate a scalar or multivariate Laplace PDF.
%
%INPUTS: x The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluation at
%          multiple points are desired, then this is a matrix with each
%          column being the a point (a vector).
%   lambda A scale parameter for the distribution. Note that with respect
%          to the traditional parameterization of a scalar distribution,
%          sqrt(2/lambda)=1/b, where b is the traditional scale parameter.
%          Here, lambda stems from the rate parameter of an exponential
%          distribution when deriving the multivariate Laplace.
%       mu The location parameter of the PDF. If the PDF is multivariate,
%          then this is a column vector.
%    Gamma The scale matrix of the PDF. This must be symmetric and positive
%          definite with determinant det(Gamma)=1. When dealing with scalar
%          distributions, this is just 1.
%
%OUTPUTS: vals The scalar value of the multivariate Laplace PDF. If
%              multiple points are passed (z is a matrix), then val is a
%              row vector.
%
%The multivariate Laplace PDF is given in Equation 9 of [1].
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(x,1);
    numX=size(x,2);
    
    vals=zeros(1,numX);
    
    for curX=1:numX
        y=x(:,curX);
        
        diff=y-mu;
        qy=invSymQuadForm(diff,Gamma);
        
        vals(curX)=((2*pi)^(-d/2))*(2/lambda)*besselk(d/2-1,sqrt((2/lambda)*qy))/(sqrt((lambda/2)*qy)^((d/2)-1));
    end

    %If any of the values are not finite, then the value is probably very
    %close to the mean, which means that the ratio involving the Bessel
    %function is Inf/Inf. For d=1, this ratio has a finite asymptotic
    %value and the asymptotic result can be substituted. For higher
    %dimensions, this ratio is just Inf, so no substitution will be made.
    if(d==1)
        sel=~isfinite(vals);
        vals(sel)=1/sqrt(2*lambda);
    end
end

function vals=CDF(x,lambda,mu,Gamma)
%%CDF Evaluate a scalar Laplace CDF at given points.
%
%INPUTS: z The points at which the PDF should be evaluated. This can be a
%          column vector, a row vector, or a matrix.
%   lambda A scale parameter for the distribution. Note that with respect
%          to the traditional parameterization of a scalar distribution,
%          sqrt(2/lambda)=1/b, where b is the traditional scale parameter.
%          Here, lambda stems from the rate parameter of an exponential
%          distribution when deriving the multivariate Laplace.
%       mu The location parameter of the PDF. If the PDF is multivariate,
%          then this is a column vector.
%    Gamma The scale matrix of the PDF. This must be symmetric and positive
%          definite with determinant det(Gamma)=1. When dealing with scalar
%          distributions, this is just 1.
%
%OUTPUTS: vals The CDF of the laplace distribution evaluated at the
%              specified points.
%
%The CDF of the Laplace distribution can be found explicitly by
%integrating the scalar PDF in Equation 5 of [1].
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin>2&&Gamma~=1)
   error('The scale matrix of a scalar Laplace distribution must be one.') 
end

b=1/sqrt(2/lambda);

sel=x<mu;
vals=zeros(size(x));

%For x<mu
vals(sel)=(1/2)*exp((x-mu)/b);

%For x>mu
vals(~sel)=1-(1/2)*exp((mu-x)/b);
end

function momentVal=momentGenFun(lambda,mu,Gamma,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the multivariate Laplace distribution. 
%              Taking the ith, jth, kth... derivative of the moment
%              generating function with respect to the first, second,
%              components of the argument and evaluating it at t=0 provides
%              third... the noncentral moment of the multivariate Laplace
%              distribution involving the ith, jth, kth power of the first
%              second, third... components of the random vector.
%
%INPUTS: lambda A scale parameter for the distribution. Note that with
%               respect to the traditional parameterization of a scalar
%               distribution, sqrt(2/lambda)=1/b, where b is the
%               traditional scale parameter. Here, lambda stems from the
%               rate parameter of an exponential distribution when deriving
%               the multivariate Laplace.
%            mu The location parameter of the PDF. If the PDF is
%               multivariate, then this is a column vector.
%         Gamma The scale matrix of the PDF. This must be symmetric and
%               positive definite with determinant det(Gamma)=1. When
%               dealing with scalar distributions, this is just 1.
%     numDerivs A numDimX1 or 1XnumDim vector indicating the number of
%               derivatives to take with respect to each of the components
%               of the argument of the moment generating function.
%               numDerivs>=0.
%             t The numDimXnumPoints argument of the moment generating
%               function at which the derivatives of the moment generating
%               function should be evaluated. If this parameter is omitted
%               or an empty matrix is passed, then a default of
%               t=zeros(numDim,1) is used.
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random vector is defined to be
%E(exp(t'*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter having the same dimensionality as the
%random variable. In Equation 13 of [1], the moment generating function of
%the multivariate Laplace distribution is shown to be
%exp(mu'*t)/(1-(lambda/2)*t'*Gambda*t)
%
%Derivatives can be evaluated systematically. First, it is clear that all
%terms of all derivatives are multiplied by exp(mu'*t). Also, all
%derivatives involve a term of the form (1-(lambda/2)*t'*Lambda*t)^(-n) for
%some positive integer n. All terms are also multiplied by multivariate
%polynomials in t. Thus, thus function evaluates the derivatives by keeping
%track of the multivariate polynomials in t and the power of the
%denominators. It thus recursively finds the derivatives one after the
%other.
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

numDim=size(mu,1);

if(nargin<5||isempty(t))
    t=zeros(numDim,1);
end

%The initial function (no derivatives) just has the exp(mu'*t) term divided
%by the denominator (so first order denominator term).
numeratPolys{1}=1;
denomOrder=1;

%Start taking derivatives with respect to each component.
for curDim=1:numDim
    %When taking the derivative of the denominator, this is the derivative
    %polynomial of lambda/2*t'*Gamma*t with respect to dimension curDim.
    %The derivative is (lambda/2)*(t'*Gamma(column curDim)+Gamma(row
    %curDim)*t). Thus, the result is linear in all of the t variables with
    %no zero-order terms.
    idxDims=2*ones(numDim,1);
    denomDerivTerm=zeros(idxDims(:)');%Hold the multidimensional polynomial.
    idxVals=ones(numDim,1);
    for curIdx=1:numDim
        idxVals(curIdx)=2;
        curEl=nDim2Index(idxDims,idxVals);
        
        denomDerivTerm(curEl)=(lambda/2)*(Gamma(curIdx,curDim)+Gamma(curDim,curIdx));
        idxVals(curIdx)=1;
    end

    %Take the required number of derivatives with respect to the current
    %component.
    for derivNum=1:numDerivs(curDim)
    %The maximum number of terms in the sum after differentiating with
    %respect to the current component is 3*numTerms, so preallocate that
    %amount of space.
        numTerms=length(numeratPolys);
    
        newNumeratPolys=cell(3*numTerms,1);
        newDenomOrders=zeros(3*numTerms,1);
        numNewAdded=0;
        
        %First, we will take the derivative with respect to the implicit
        %exp(mu'*t) out front. This multiplied all of the polynomials by
        %mu(curDim) and it does not change the orders of the denominators.
        %if mu(curDim) is zero, then no terms are added.
        if(mu(curDim)~=0)
            newNumeratPolys(1:numTerms)=numeratPolys;
            newDenomOrders(1:numTerms)=denomOrder;
            
            for curPoly=1:numTerms
                newNumeratPolys{curPoly}=newNumeratPolys{curPoly}*mu(curDim);
            end
            numNewAdded=numTerms;
        end
        
        %In the second half of the product rule, we must taken the
        %derivative of all of the polynomials. By the chain rule, this
        %means two derivatives: one with respect to the numerator, and one
        %with respect to the denominator.
        for curPoly=1:numTerms
            %First add the term that is with respect to the numerator. This
            %is just the derivative of the multivariate polynomial in the
            %numerator; the denominator is unchanged.
            theDer=polyDerMultiDim(numeratPolys{curPoly},curDim);
            
            %If the derivative did not reduce the polynomial to zero, then
            %keep it.
            if(~isscalar(theDer)||theDer~=0)
                numNewAdded=numNewAdded+1;
                newNumeratPolys{numNewAdded}=theDer;
                newDenomOrders(numNewAdded)=denomOrder(curPoly);
            end
            
            %Next, take the derivative of the denominator term. This
            %multiplies the existing polynomial in the numerator by a new
            %term and it increases the order of the denominator.
            numNewAdded=numNewAdded+1;
            newNumeratPolys{numNewAdded}=denomOrder(curPoly)*convn(denomDerivTerm,numeratPolys{curPoly});
            newDenomOrders(numNewAdded)=denomOrder(curPoly)+1;
        end
        
        %We want to combine the terms with common denominator orders. The
        %idea is that reducing the number of terms should reduce finite
        %precision errors if many derivatives are taken and there are lots
        %of terms. We sort by order, then we add all of the terms of
        %the same order. First, we reduce the size of the arrays to the
        %number actually added so as to simplify things.
        newDenomOrders=newDenomOrders(1:numNewAdded);
        newNumeratPolys=newNumeratPolys(1:numNewAdded);
        [newDenomOrders,idxSort]=sort(newDenomOrders,'ascend');
        newNumeratPolys=newNumeratPolys(idxSort);
        
        numTerms=length(newDenomOrders);
        numeratPolys=cell(numTerms,1);
        denomOrder=zeros(numTerms,1);
        numAdded=0;
        
        curOrder=newDenomOrders(1);
        curOrderIdx=1;
        for curTerm=2:numTerms
            %If the end of a stretch of equal orders has been reached.
            if(newDenomOrders(curTerm)~=curOrder)
                combinedPoly=newNumeratPolys{curOrderIdx};
                for curIdx=(curOrderIdx+1):(curTerm-1)
                    combinedPoly=polySumMultiDim(combinedPoly,newNumeratPolys{curIdx});
                end
                
                numAdded=numAdded+1;
                numeratPolys{numAdded}=combinedPoly;
                denomOrder(numAdded)=curOrder;

                curOrder=newDenomOrders(curTerm);
                curOrderIdx=curTerm;
            end
        end

        %Add the final term in the sequence
        combinedPoly=newNumeratPolys{curOrderIdx};
        for curIdx=(curOrderIdx+1):numTerms
            combinedPoly=polySumMultiDim(combinedPoly,newNumeratPolys{curIdx});
        end
        numAdded=numAdded+1;
        numeratPolys{numAdded}=combinedPoly;
        denomOrder(numAdded)=curOrder;
        
        %Shrink to fit.
        numeratPolys=numeratPolys(1:numAdded);
        denomOrder=denomOrder(1:numAdded);
    end
end

%Now that we have constructed the terms of the multivariate series,
%evaluate it at the given value of t to get the moments
numTerms=length(denomOrder);
numPoints=size(t,2);

momentVal=zeros(numPoints,1);

for curPoint=1:numPoints
    %The value of the denominator that is raised to various powers in the
    %different terms.
    denomVal=1-(lambda/2)*(t(:,curPoint)'*Gamma*t(:,curPoint));
    for curTerm=1:numTerms
        numerVal=polyValMultiDim(numeratPolys{curTerm},t(:,curPoint));
        momentVal(curPoint)=momentVal(curPoint)+numerVal/denomVal^(denomOrder(curTerm));
    end
    %The term that multiplies all of the fractions in the series.
    momentVal(curPoint)=exp(mu'*t(:,curPoint))*momentVal(curPoint);
end
end

function cumVal=cumGenFun(lambda,mu,Gamma,numDerivs,t)
%%CUMGENFUN Evaluate the cumulant generating function (or one of its
%           derivatives) of the multivariate Laplace distribution. Taking
%           the ith, jth, kth... derivative of the cumulant generating
%           function with respect to the first, second, third...
%           components of the argument and evaluating it at t=0 provides
%           the cumulant of the multivariate Laplace distribution involving
%           the ith, jth, kth power of the first second, third...
%           components of the random vector. The cumulant generating
%           function is the natural logarithm of the moment generating
%           function.
%
%INPUTS: lambda A scale parameter for the distribution. Note that with
%               respect to the traditional parameterization of a scalar
%               distribution, sqrt(2/lambda)=1/b, where b is the
%               traditional scale parameter. Here, lambda stems from the
%               rate parameter of an exponential distribution when deriving
%               the multivariate Laplace.
%            mu The location parameter of the PDF. If the PDF is
%               multivariate, then this is a column vector.
%         Gamma The scale matrix of the PDF. This must be symmetric and
%               positive definite with determinant det(Gamma)=1. When
%               dealing with scalar distributions, this is just 1.
%     numDerivs A numDimX1 or 1XnumDim vector indicating the number of
%               derivatives to take with respect to each of the components
%               of the argument of the moment generating function.
%               numDerivs>=0.
%             t The numDimXnumPoints argument of the moment generating
%               function at which the derivatives of the moment generating
%               function should be evaluated. If this parameter is omitted
%               or an empty matrix is passed, then a default of
%               t=zeros(numDim,1) is used.
%
%OUTPUTS: cumVal A numPointsX1 vector of the values of the derivatives
%                of the cumulant generating function given at the points
%                in t or at t=0 if t is omitted.
%
%Cumulants are useful in interpolating probability distributions through
%the use of, for example, an Edgeworth series. The cumulant generating
%function is defined as the natural logarithm of the moment generating
%function. In Equation 13 of [1], the moment generating function of
%the multivariate Laplace distribution is shown to be
%exp(mu'*t)/(1-(lambda/2)*t'*Gamma*t)
%Thus, the cumulant generating function is just
%mu'*t+log(1-(lambda/2)*t'*Gamma*t)
%
%Derivatives are evaluated systematically in a similar manner to what was
%done in the LambdaD.momentGenFun function . The first derivative with
%respect to component i1 of t is
%mu(i1)-(lambda/2)*D{t'*Lambda*t,i1}/(1-(lambda/2)*(t'*Gamma*t))
%where D{x[i1],i1} is an operator to differentiate the argument x with
%respect to the variable i1. After the first derivative, the mu term goes
%away and the subsequent terms all contain powers of
%1-(lambda/2)*t'*Lambda*t in the denominator and multivariate polynomials
%in the numerator. This function keeps track of the polynomial terms and
%the powers of the denominator.
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300- 303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

numDim=size(mu,1);

if(nargin<5||isempty(t))
    t=zeros(numDim,1);
end

numPoints=size(t,2);
cumVal=zeros(numPoints,1);

%We need to know where to start taking derivatives.
startIdx=find(numDerivs);

%If no derivatives are taken, then evaluate the cumulant generating
%function directly including the logarithmic term.
if(isempty(startIdx))
    for curPoint=1:numPoints
    	cumVal(curPoint)=mu'*t(:,curPoint)-log(1-(lambda/2)*(t(:,curPoint)'*Gamma*t(:,curPoint)));
    end

	return;
end

%If we are here, then at least one derivative is taken. If only a single
%derivative with respect to a single component is taken, then we treat it
%as a special case, because unlike all other derivatives, it includes a mu
%term.
if(numDerivs(startIdx(1))==1&&length(startIdx)==1)
    derivIdx=startIdx(1);
    for curPoint=1:numPoints
    	cumVal(curPoint)=mu(derivIdx)+(lambda/2)*(t(:,curPoint)'*Gamma(:,derivIdx)+Gamma(derivIdx,:)*t(:,curPoint))/(1-(lambda/2)*(t(:,curPoint)'*Gamma*t(:,curPoint)));
    end

    return;
end

%If we are here, then more than one derivative is taken and thus there is
%no extra mu term in front of the series. We will manually insert
%the value of the first derivative (without the mu term). Then, we can
%systematically evaluate the other derivatives in the same manner as in
%LambdaD.momentGenFun if more derivatives are taken.
idxDims=2*ones(numDim,1);
firstDerivTerm=zeros(idxDims(:)');%Hold the multidimensional polynomial.
idxVals=ones(numDim,1);
derivIdx=startIdx(1);
for curIdx=1:numDim
    idxVals(curIdx)=2;
    curEl=nDim2Index(idxDims,idxVals);

    firstDerivTerm(curEl)=(lambda/2)*(Gamma(curIdx,derivIdx)+Gamma(derivIdx,curIdx));
    idxVals(curIdx)=1;
end

numeratPolys{1}=firstDerivTerm;
denomOrder=1;

%Start taking derivatives.
for curDim=startIdx(1):numDim
    %When taking the derivative of the denominator, this is the derivative
    %polynomial of lambda/2*t'*Gamma*t with respect to dimension curDim.
    %The derivative is (lambda/2)*(t'*Gamma(column curDim)+Gamma(row
    %curDim)*t). Thus, the result is linear in all of the t variables with
    %no zero-order terms.
    idxDims=2*ones(numDim,1);
    denomDerivTerm=zeros(idxDims(:)');%Hold the multidimensional polynomial.
    idxVals=ones(numDim,1);
    for curIdx=1:numDim
        idxVals(curIdx)=2;
        curEl=nDim2Index(idxDims,idxVals);
        
        denomDerivTerm(curEl)=(lambda/2)*(Gamma(curIdx,curDim)+Gamma(curDim,curIdx));
        idxVals(curIdx)=1;
    end

    %Take the required number of derivatives with respect to the current
    %component.
    if(curDim==startIdx(1))
        derivStart=2;
    else
        derivStart=1;
    end

    for derivNum=derivStart:numDerivs(curDim)
    %The maximum number of terms in the sum after differentiating with
    %respect to the current component is 2*numTerms, so preallocate that
    %amount of space.
        numTerms=length(numeratPolys);
    
        newNumeratPolys=cell(2*numTerms,1);
        newDenomOrders=zeros(2*numTerms,1);
        numNewAdded=0;

        %We must taken the derivative of both of the polynomials in the
        %ratio. By the chain rule, this means two derivatives: one with
        %respect to the numerator, and one with respect to the denominator.
        for curPoly=1:numTerms
            %First add the term that is with respect to the numerator. This
            %is just the derivative of the multivariate polynomial in the
            %numerator; the denominator is unchanged.
            theDer=polyDerMultiDim(numeratPolys{curPoly},curDim);
            
            %If the derivative did not reduce the polynomial to zero, then
            %keep it.
            if(~isscalar(theDer)||theDer~=0)
                numNewAdded=numNewAdded+1;
                newNumeratPolys{numNewAdded}=theDer;
                newDenomOrders(numNewAdded)=denomOrder(curPoly);
            end
            
            %Next, take the derivative of the denominator term. This
            %multiplies the existing polynomial in the numerator by a new
            %term and it increases the order of the denominator.
            numNewAdded=numNewAdded+1;
            newNumeratPolys{numNewAdded}=denomOrder(curPoly)*convn(denomDerivTerm,numeratPolys{curPoly});
            newDenomOrders(numNewAdded)=denomOrder(curPoly)+1;
        end
        
        %We want to combine the terms with common denominator orders. The
        %idea is that reducing the number of terms should reduce finite
        %precision errors if many derivatives are taken and there are lots
        %of terms. We sort by order, then we add all of the terms of
        %the same order. First, we reduce the size of the arrays to the
        %number actually added so as to simplify things.
        newDenomOrders=newDenomOrders(1:numNewAdded);
        newNumeratPolys=newNumeratPolys(1:numNewAdded);
        [newDenomOrders,idxSort]=sort(newDenomOrders,'ascend');
        newNumeratPolys=newNumeratPolys(idxSort);
        
        numTerms=length(newDenomOrders);
        numeratPolys=cell(numTerms,1);
        denomOrder=zeros(numTerms,1);
        numAdded=0;
        
        curOrder=newDenomOrders(1);
        curOrderIdx=1;
        for curTerm=2:numTerms
            %If the end of a stretch of equal orders has been reached.
            if(newDenomOrders(curTerm)~=curOrder)
                combinedPoly=newNumeratPolys{curOrderIdx};
                for curIdx=(curOrderIdx+1):(curTerm-1)
                    combinedPoly=polySumMultiDim(combinedPoly,newNumeratPolys{curIdx});
                end
                
                numAdded=numAdded+1;
                numeratPolys{numAdded}=combinedPoly;
                denomOrder(numAdded)=curOrder;

                curOrder=newDenomOrders(curTerm);
                curOrderIdx=curTerm;
            end
        end

        %Add the final term in the sequence
        combinedPoly=newNumeratPolys{curOrderIdx};
        for curIdx=(curOrderIdx+1):numTerms
            combinedPoly=polySumMultiDim(combinedPoly,newNumeratPolys{curIdx});
        end
        numAdded=numAdded+1;
        numeratPolys{numAdded}=combinedPoly;
        denomOrder(numAdded)=curOrder;
        
        %Shrink to fit.
        numeratPolys=numeratPolys(1:numAdded);
        denomOrder=denomOrder(1:numAdded);
    end
end

%Now that we have constructed the terms of the multivariate series,
%evaluate it at the given value of t to get the cumulants.
numTerms=length(denomOrder);
numPoints=size(t,2);

cumVal=zeros(numPoints,1);

for curPoint=1:numPoints
    %The value of the denominator that is raised to various powers in the
    %different terms.
    denomVal=1-(lambda/2)*(t(:,curPoint)'*Gamma*t(:,curPoint));
    for curTerm=1:numTerms
        numerVal=polyValMultiDim(numeratPolys{curTerm},t(:,curPoint));
        cumVal(curPoint)=cumVal(curPoint)+numerVal/denomVal^(denomOrder(curTerm));
    end
end

end

function vals=rand(N,lambda,mu,Gamma)
%%RAND Generate multivariate Laplace random variables.
%
%INPUTS: N The number of random variables to generate.   
%   lambda A scale parameter for the distribution. Note that with respect
%          to the traditional parameterization of a scalar distribution,
%          sqrt(2/lambda)=1/b, where b is the traditional scale parameter.
%          Here, lambda stems from the rate parameter of an exponential
%          distribution when deriving the multivariate Laplace.
%      mu  The location parameter of the PDF. If the PDF is multivariate,
%          then this is a column vector.
%    Gamma The scale matrix of the PDF. This must be symmetric and positive
%          definite with determinant det(Gamma)=1. When dealing with scalar
%          distributions, this is just 1.
%
%OUTPUTS: vals An numDimXN matrix of random instances of the multivariate
%              Laplace distribution. The dimensionality (numDim) is
%              inferred from the dimensionality fo the location parameter
%              mu.
%
%The generation of random variables follows from the definition of the
%multivariate Laplace distribution in Section II of [1].
%
%REFERENCES:
%[1] T. Elftoft, T. Kim, and T.-W. Lee, "On the multivariate Laplace
%    distribution," IEEE Signal Processing Letters, vol. 13, no. 5, pp.
%    300-303, May 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    numDim=size(mu,1);
    Z=ExponentialD.rand([1,N],1/lambda);
    SGamma=chol(Gamma,'lower');
    
    vals=bsxfun(@plus,mu,bsxfun(@times,sqrt(Z),SGamma*randn(numDim,N)));
end

function entropyVal=entropy(lambda)
%%ENTROPY Obtain the differential entropy of the scalar Laplace
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: lambda A scale parameter for the distribution. Note that with
%               respect to the traditional parameterization of a scalar
%               distribution, sqrt(2/lambda)=1/b, where b is the
%               traditional scale parameter.
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
    
    b=1/sqrt(2/lambda);

    entropyVal=log(2*b*exp(1));
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
