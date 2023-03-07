classdef GaussianD
%%GAUSSIAND Functions to handle the scalar and multivariate Gaussian
%           distribution.
%Implemented methods are: mean,cov, PDF, PDFI, PDFIDerivs (for multivariate
%                         derivatives of the PDF), PDFS, logPDFS,
%                         PDFSGradHessVechS (for the gradient and Hessian
%                         of the elements of a lower-triangular square-root
%                         of the covariance matrix), CDF (for scalar
%                         distributions), invCDF (for scalar
%                         distributions), normProdDist, normConvDist (for
%                         scalar distributions), momentGenFun
%                         (multivariate, including derivatives), cumGenFun
%                         (multivariate, including derivatives), rand,
%                         randS, integralOverRegion, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(mu)
%%MEAN Obtain the mean of the Gaussian distribution.
%
%INPUTS: mu The mean of the PDF. If the PDF is multivariate, then this is
%           a numDimX1 column vector.
%
%OUTPUTS: val The numDimX1 mean of the Gaussian distribution.
%
%The Gaussian distribution is parameterized by its mean and covariance
%matrix. Thus, this function just returns the mean it is given.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end

function val=cov(Sigma)
%%COV Obtain the covariance matrix of the Gaussian distribution (the
%     variance if scalar).
%
%INPUTS: Sigma The variance (if scalar) or covariance matrix (if
%               multidimensional) of the PDF. The variance cannot be zero
%               and the covariance matrix cannot be singular.
%
%OUTPUTS: val The covariance matrix of the Gaussian distribution.
%
%The Gaussian distribution is parameterized by its mean and covariance
%matrix. Thus, this function just returns the covariance matrix it is
%given.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=Sigma;
end

function vals=PDF(z,mu,Sigma)
%%PDF Evaluate a scalar or multivariate Gaussian (normal) PDF at specified
%     points given the mean and the covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluation at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDimX1 column vector. If omitted or an empty matrix is passed,
%          a zero mean is used.
%    Sigma The variance (if scalar) or numDimXnumDim covariance matrix
%          (if multidimensional) of the PDF. The variance cannot be zero
%          and the covariance matrix cannot be singular. If omitted or an
%          empty matrix is passed, the identity matrix is used as the
%          covariance matrix.
%
%OUTPUTS: vals The scalar values of the normal PDF with mean mu and
%              covariance matrix Sigma evaluated at the points in z. If
%              multiple points are passed (z is a matrix), then val is a
%              row vector.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    numDim=size(z,1);
    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(Sigma))
       Sigma=eye(numDim,numDim); 
    end
    
    diff=bsxfun(@minus,z,mu);
    vals=(2*pi)^(-numDim/2)*(det(Sigma))^(-1/2)*exp(-0.5*invSymQuadForm(diff,Sigma));
end

function vals=logPDF(z,mu,Sigma)
%%LOGPDF Evaluate the natural logarithm of a scalar or multivariate
%        Gaussian (normal) PDF at a certain points given the mean and the
%        covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluation at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDimX1 column vector. If omitted or an empty matrix is passed,
%          a zero mean is used.
%    Sigma The variance (if scalar) or numDimXnumDim covariance matrix
%          (if multidimensional) of the PDF. The variance cannot be zero
%          and the covariance matrix cannot be singular. If omitted or an
%          empty matrix is passed, the identity matrix is used as the
%          covariance matrix.
%
%OUTPUTS: vals The scalar values of the natural logarithm of the normal PDF
%              with mean mu and covariance matrix Sigma evaluated at the
%              points in z. If multiple points are passed (z is a matrix),
%              then val is a row vector.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    numDim=size(z,1);
    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(Sigma))
       Sigma=eye(numDim,numDim); 
    end
    
    diff=bsxfun(@minus,z,mu);
    vals=-(1/2)*log(det(2*pi*Sigma))-(1/2)*invSymQuadForm(diff,Sigma);
end

function vals=PDFI(z,mu,SigmaInv,SigmaInvDet)
%%PDFI Evaluate a scalar or multivariate Gaussian (normal) PDF at specified
%      points given the mean and the inverse of the covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluations at
%          multiple points are desired, then this is a numDimXN matrix with
%          each column being the a point (a vector).
%       mu The mean of the PDF. If the PDF is multivariate, then this is a
%          numDimX1 column vector. If omitted or an empty matrix is passed,
%          a zero mean is used.
% SigmaInv The inverse variance (if scalar) or numDimXnumDim inverse
%          covariance matrix (if multidimensional) of the PDF. SigmaInv can
%          be singular. If omitted or an empty matrix is passed, the
%          identity matrix is used as the covariance matrix.
% SigmaInvDet Optionally, a length-N set of determinants of the matrices
%          in SigmaInv can be passed so as to speed up the computation. If
%          omitted or an empty matrix is passed, determinants will be taken
%          as needed.
%
%OUTPUTS: val The scalar value of the normal PDF with mean mu and inverse
%             covariance matrix SigmaInv evaluated at the point z. If
%             multiple points are passed (z is a matrix), then val is a row
%             vector.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numPoints=size(z,2);
    vals=zeros(1,numPoints);
    n=size(z,1);
    if(nargin<2||isempty(mu))
        mu=zeros(numDim,1);
    end
    
    if(nargin<3||isempty(SigmaInv))
       SigmaInv=eye(numDim,numDim); 
    end
    
    if(nargin<4||isempty(SigmaInvDet))    
        constVal=sqrt(det(SigmaInv)/(2*pi)^(n));
    else
        constVal=sqrt(SigmaInvDet/(2*pi)^(n));
    end
    for curPoint=1:numPoints
        diff=z(:,curPoint)-mu;
        %Note that det(A^(-1))=1/det(A) and that det(a*A)=a^n*det(A), where
        %a is a scalar and A is an nXn matrix.

        vals(curPoint)=constVal*exp(-0.5*(diff'*SigmaInv*diff));
    end
end

function [PDFDerivVal,coeffPolyPart]=PDFIDerivs(mu,SigmaInv,numDerivs,x)
%%PDFIDERIVS Compute derivatives of the multivariate Gaussian normal PDF at
%            given the mean and the inverse of the covariance matrix.
%            Derivatives are taken with respect to components of the
%            argument of the PDF, x.
%
%INPUTS: mu The mean of the PDF. If the PDF is multivariate, then this is a
%           numDimX1 column vector.
%  SigmaInv The inverse variance (if scalar) or numDimXnumDim inverse
%           covariance matrix (if multidimensional) of the PDF.
%           SigmaInv can be singular.
% numDerivs A numDimX1 or 1XnumDim vector indicating the number of
%           derivatives to take with respect to each of the dimensions of
%           the state.numDerivs>=0.
%         x The numDimXnumPoints argument of the PDF at which the
%           derivatives of the PDF should be evaluated. If this parameter
%           is omitted or an empty matrix is passed, then a default of
%           x=zeros(numDim,1) is used.
%
%OUTPUTS: PDFDerivVal A numPointsX1 vector of the values of the derivatives
%                   of the PDF function given at the points in x or at x=0
%                   if x is omitted.
%     coeffPolyPart A hypermatrix taking numDim indices that can be
%                   evaluated using the polyValMultiDim at different values
%                   of x to get the polynomial coefficient that multiplies
%                   the exponential term in the given set of derivatives
%                   of the PDF.
%
%All derivatives of the multivariate normal distribution have the form
%sqrt((2*pi)^(-2)*det(SigmaInv))*exp(-1/2*((x-mu)'*SigmaInv*(x-mu))*...
%(polynomial)
%The polynomial can be found using the chain rule each time a derivative is
%taken. The derivative of the exponential term with respect to x(i), the
%ith component of the state is
%sqrt((2*pi)^(-2)*det(SigmaInv))*exp(-1/2*((x-mu)'*SigmaInv*(x-mu)) times
%SigmaInv(i,:)*mu-SigmaInv(i,:)*x.
%This function returns coeffPolyPart, the coefficients of the multivariate
%polynomial that multiplies
%sqrt((2*pi)^(-2)*det(SigmaInv))*exp(-1/2*((x-mu)'*SigmaInv*(x-mu)) when
%computing the derivatives in addition to returning the value of the PDF at
%any desired points.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C. 

numIdx=size(mu,1);

coeffPolyPart=1;

%i is the current dimension we are differentiating.
for i=1:numIdx
    %basePoly shall be the multivariate polynomial that is added every time
    %the exponential term in the moment generating function is
    %differentiated with respect to the ith index. This makes basePoly the
    %derivative with respect to the ith element of t of the argument of the
    %exponent in the multivariate normal PDF.
    
    %Allocate space for the polynomial and make it have the correct shape.
    basePoly=reshape(zeros(2^numIdx,1),2*ones(1,numIdx));
    
    %The additive term
    basePoly(1)=SigmaInv(i,:)*mu;
    
    %The term multiplied by x consists of elements from the inverse
    %covariance matrix. This just consists of terms in the ith row of the
    %inverse covariance matrix (the matrix is symmetric).
    idxVec=ones(numIdx,1);
    %The dimensionalities of all of the variables for the nDim2Index
    %function.
    dims=2*ones(numIdx,1);
    for curDim=1:numIdx
        idxVec(curDim)=2;
        basePoly(nDim2Index(dims,idxVec))=-SigmaInv(i,curDim);
        idxVec(curDim)=1;
    end
    
    %Now, we enter into a loop to evaluate derivatives of the PDF with
    %respect to the current dimension.
    for derivsLeft=numDerivs(i):-1:1
        %The chain rule means that there are two terms to consider (both of
        %which are multiplied by the same exponential term).
        %The first term is the current coeffPolyPart times the derivative 
        %of the exponential term with respect to the ith dimensions. This
        %is the product of coeffPolyPart and basePoly. The second term is
        %the exponential term times the derivative of coeffPolyPart. The
        %two terms then must be added.
        term1=convn(basePoly,coeffPolyPart);%Multiply the polynomials.
        term2=polyDerMultiDim(coeffPolyPart,i);
        
        %Add the multivariate polynomials.
        coeffPolyPart=polySumMultiDim(term1,term2);
    end
end

%Get rid of redundant parts that arose due to the multiplication.
coeffPolyPart=shrinkMultiDimPoly2Fit(coeffPolyPart);

%coeffPolyPart now contains the multivariate polynomial that is multiplied
%by the exponential term. If no value of t is given, then just evaluate it
%at x=0. In this instance, the exponential term is zero and only the
%constant term from the polynomial appears.
if(nargin<4||isempty(x))
    PDFDerivVal=GaussianD.PDFI([0;0],mu,SigmaInv)*coeffPolyPart(1);
else
    numPoints=size(x,2);
    PDFDerivVal=zeros(numPoints,1);
    for curPoint=1:numPoints
        xCur=x(:,curPoint);
        PDFDerivVal(curPoint)=GaussianD.PDFI(xCur,mu,SigmaInv)*polyValMultiDim(coeffPolyPart,xCur);
    end
end
end

function val=PDFS(z,mu,S)
%%PDFS Evaluate a scalar or multivariate Gaussian (normal) PDF at specifed
%      points given the mean and the lower-triangular square root of the
%      covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluations at
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
%             root covariance matrix S evaluated at the points in z. If
%             multiple points are passed (z is a matrix), then val is a row
%             vector.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

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
    %those terms is not unique due to the squaring.
    val = (1/((2*pi)^(n/2)*abs(det(S))))*exp(-0.5*sum(diff.*diff,1)); 
end

function val=logPDFS(z,mu,S)
%%LOGPDFS Evaluate the natural logarithm of a scalar or multivariate
%         Gaussian (normal) PDF at specified points given the mean and the
%         lower-triangular square root of the covariance matrix.
%
%INPUTS: z The points at which the PDF should be evaluated. If the PDF is
%          multivariate, then this is a column vector. If evaluations at
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
%OUTPUTS: val The scalar value(s) of the natural logarithm of the normal
%             PDF with mean mu and square root covariance matrix S
%             evaluated at the points in z. If multiple points are passed
%             (z is a matrix), then val is a row vector.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

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
    %those terms is not unique due to the squaring.
    % val=-(n/2)*log((2*pi))-log(abs(det(S)))-0.5*sum(diff.*diff,1); 
    val=log((1/((2*pi)^(n/2)*abs(det(S))))*exp(-0.5*sum(diff.*diff,1)));
end

function [grad,CDetGrad,Hess,CDetHess]=PDFSGradHessVechS(x,mu,C)
%%PDFGRADHESSVECHS Find the gradient and (if requested Hessian) of the
%            normal (Gaussian) probability density function (PDF), taken
%            with respect to the vech(C), where C is the lower-triangular
%            square root of the covariance matrix of the distribution.
%
%INPUTS: x The dXnumPoints points at which the normal PDF is considered.
%       mu The dX1 mean of the normal PDF.
%        C The dXd lower-triangular square root covariance matrix of the
%          normal PDF. This cannot be singular.
%
%OUTPUTS: grad The gradient of the normal PDF with respect to vech(C)
%              (vector of first partial derivatives  with respect to the
%              elements of vech(C)). This is an (n*(n+1)/2)XnumPoints set
%              of vectors.
%     CDetGrad The (n*(n+1)/2)X1 gradient of 1/sqrt(C*C') with respect to
%              vech(C). This term is needed to compute grad and Hess and is
%              often needed in algorithms using grad and Hess.
%         Hess The Hessian of the normal PDF with respect to vech(C)
%              (vector of second partial derivatives  with respect to the
%              elements of vech(C)). This is an
%              (n*(n+1)/2)X(n*(n+1)/2)XnumPoints set of symmetric matrices.
%              Element i,j in a matrix is the second derivative with
%              respect to elements i and j of vech(C).
%     CDetHess The (n*(n+1)/2)X1 Hessian of 1/sqrt(C*C') with respect to
%              vech(C). This term is needed to compute grad and Hess and is
%              often needed in algorithms using grad and Hess.
%
%The gradient and Hessian provided by this function play a role in
%multivariate kernel bandwidth estimation algorithms, such as the
%cross-validation bandwidth estimation algorithms of [1], which optimize
%over the components of a lower-triangular square root matrix to avoid
%obtaining invalid covariance matrix estimates.
%
%This function is implemented based on the rules of matrix calculus.
%
%A 5-dimensional example.
% c=[17;23;4;10;11;5;6;12;18;13;19;25;21;2;9];
% C=vech2Mat(c,0);
% x=[5;-10;33;24;-18];
% mu=[0;0;0;0;0];
% [grad,Hess]=GaussianD.PDFSGradHessVechS(x,mu,C)
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Cross-validation bandwidth matrices for
%    multivariate kernel density estimation," Scandinavian Journal of
%    Statistics, vol. 32, no. 3, pp. 485-506, Sep. 2005.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Since the derivatives are not with respect to x, it does not matter if the
%distribution is not zero-mean; just shift the value.
x=bsxfun(@minus,x,mu);

%Problem dimensionality
d=size(x,1);
numPoints=size(x,2);

vechC=vech(C);
numVechEls=length(vechC);

CDetRoot=1/sqrt(det(C*C'));

%Find the gradient of 1/sqrt(C*C') with respect to vech(C). Only the
%derivatives corresponding to the diagonals are non-zero.
CDetGrad=vech(diag(-CDetRoot./diag(C)));

CInv=inv(C);

%The exponential term
expTerm=exp(-0.5*invSymQuadForm(x,C,1));

%Find the gradient of C^(-1) with respect to vech(C).
CInvGrad=zeros(d,d,numVechEls);
CDeriv=zeros(d,d);
curEl=1;
for j=1:d
    for i=j:d
        CDeriv(i,j)=1;
        
        CInvGrad(:,:,curEl)=CInv*CDeriv*CInv;
        curEl=curEl+1;
        CDeriv(i,j)=0;
    end
end

%Find the gradient of (C*C')^(-1) with respect to vech(C).
CCpInvGrad=zeros(d,d,numVechEls);
for i=1:numVechEls
    temp=CInv'*CInvGrad(:,:,i);
    CCpInvGrad(:,:,i)=temp+temp';
end

%Find the gradient of the exponential term with respect to vech(C).
expTermGrad=zeros(numVechEls,numPoints);
for i=1:numVechEls
    expTermGrad(i,:)=(1/2)*expTerm.*sum(bsxfun(@times,x,CCpInvGrad(:,:,i)*x),1);
end

%Now, find the derivatives of the PDF with respect to every element of
%of vech(C).
grad=zeros(numVechEls,numPoints);
for i=1:numVechEls
    grad(i,:)=CDetGrad(i)*expTerm+CDetRoot*expTermGrad(i,:);
end
grad=grad/(2*pi)^(d/2);

if(nargout>2)%If the Hessian is desired.
    Hess=zeros(numVechEls,numVechEls,numPoints);
    CDetHess=zeros(numVechEls,numVechEls);
    for m=1:numVechEls
        CDerivM=zeros(d,d);
        [i,j]=vechInd2Sub(d,m);%--derivative indices
        CDerivM(i,j)=1;
        
        %This term is used in the Hessian of exp((1/2)*x'*inv(C*C')*x)
        term1H=sum(bsxfun(@times,x, CCpInvGrad(:,:,m)*x),1);
        for n=m:numVechEls
            CDerivN=zeros(d,d);
            [k,l]=vechInd2Sub(d,n);%--derivative indices
            CDerivN(k,l)=1;

            %The Hessian of inv(C) with respect to indices i,j and k,l.
            CInvHess=CInv*(CDerivN*CInv*CDerivM+CDerivM*CInv*CDerivN)*CInv;

            term1=CInv'*CInvHess;
            term2=CInvGrad(:,:,m)'*CInvGrad(:,:,n);
            %The Hessian of inv(C*C') with respect to indices i,j and k,l.
            CCpInvHess=term1+term1'+term2+term2';
            %The Hessian of exp((1/2)*x'*inv(C*C')*x)
            term2=sum(bsxfun(@times,x, CCpInvGrad(:,:,n)*x),1);
            term3=sum(bsxfun(@times,x,CCpInvHess*x),1);
            expTermHess=(expTerm/2).*((1/2)*term1H.*term2-term3);
            
            Hess(m,n,:)=expTermGrad(m,:)*CDetGrad(n)+expTermGrad(n,:)*CDetGrad(m)+CDetRoot*expTermHess;
            if(i==j&&l==k)%The second derivative term of the Hessian
                CDetHess(m,n)=CDetRoot/(C(i,i)*C(k,k));
                if(i==k)
                    CDetHess(m,n)=CDetHess(m,n)*2;
                end
                %Due to the symmetry of the Hessian.
                CDetHess(n,m)=CDetHess(m,n);
                
                Hess(m,n,:)=Hess(m,n,:)+reshape(expTerm*CDetHess(m,n),1,1,numPoints);
            end

            Hess(m,n,:)=Hess(m,n,:)*(2*pi)^(-d/2);
            %Due to the symmetry of the Hessian, the upper triangular part
            %is also known.
            Hess(n,m,:)=Hess(m,n,:);
        end
    end
end

end

function val=CDF(z,mu,varVal)
%%CDF Evaluate cumulative distribution function (CDF) of a a scalar
%     Gaussian (normal) distribution at a specified points given the mean
%     and the variance, or for a normal(0,1) distribution if the mean and
%     variance are omitted.
%
%INPUTS: z A matrix of the point(s) at which the CDF should be evaluated.
%       mu The mean of the distribution. If omitted or an empty matrix is
%          passed, a mean of 0 is used.
%   varVal The variance of the distribution. If omitted or an empty matrix
%          is passed, a variance of 1 is used.
%
%OUTPUTS: val The scalar value(s) of the normal CDF with mean mu and
%             variance varVal evaluated at the point(s) z.
%
%This just uses the relation between the normal CDF and the error function
%along with the erf function in Matlab.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

	if(nargin<2||isempty(mu))
        mu=0; 
	end

	if(nargin<3||isempty(varVal))
        varVal=1; 
	end

    x=(z-mu)/sqrt(varVal);
    val=(1+erf(x/sqrt(2)))/2;
end

function val=invCDF(prob,mu,varVal)
%%CDF Evaluate the inverse CDF of a scalar Gaussian (normal) distribution
%     at a given point given the mean and the variance, or for a
%     normal(0,1) distribution if the mean and variance are omitted. When
%     considering a normal(0,1) distribution, this is also known as the
%     probit function.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          mu The mean of the distribution. If omitted or an empty matrix
%             is passed, a mean of 0 is used.
%      varVal The variance of the distribution. If omitted or an empty
%             matrix is passed, a variance of 1 is used.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%This just uses the relation between the normal CDF and the error function
%along with the erfinv (inverse error function) command in Matlab.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    if(nargin<2||isempty(mu))
       mu=0; 
    end
    
    if(nargin<3||isempty(varVal))
       varVal=1; 
    end

    val=sqrt(2*varVal)*erfinv(2*prob-1)+mu;
end

function [mu,SigmaInv,multiConst]=normProdDist(mu1,SigmaInv1,mu2,SigmaInv2)
%%NORMPRODDIST The product of two multivariate normal distributions is an
%              unnormalized Gaussian distribution. This finds the
%              parameters of the product distribution.
%
%INPUTS: mu1,SigmaInv1 The mean and inverse of the covariance matrix
%                      (inverse of the variance for a scalar distribution)
%                      of the first normal distribution.
%        mu2,SigmaInv2 The mean and inverse of the covariance matrix
%                      (inverse of the variance for a scalar distribution)
%                      of the second normal distribution.
%
%OUTPUTS: mu, SigmaInv The mean and inverse covariance matrix of the
%                      product distribution.
%           multiConst The product distribution is not normalized. This is
%                      the multiplicative constant that is multiplied by a
%                      normalized distribution.
%
%The derivation of the product of two normal PDFs is a standard exercise in
%many statistics classes. The product distribution is also given
%explicitly in [1].
%
%Note that while SigmaInv1 and SigmaInv2 can each be singular, the sum must
%be non-singular. Note that for products involving distributions with
%distant means and small variances, multiConst might be numerically zero.
%
%REFERENCES:
%[1] K. B. Petersen and M. S. Pedersen, "The matrix cookbook," Technical
%    University of Denmark, Tech. Rep., 15 Nov. 2012. [Online]. Available:
%    http://www2.imm.dtu.dk/pubdb/views/publication details.php?id=3274
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C. 

SigmaInv=SigmaInv1+SigmaInv2;
mu=SigmaInv\(SigmaInv1*mu1+SigmaInv2*mu2);

multiConst=GaussianD.PDFI(mu1,mu2,SigmaInv);
end

function [mu,Sigma]=normConvDist(mu1,Sigma1,mu2,Sigma2)
%%NORMCONVDIST The convolution of two scalar normal distributions is also a
%              normal distribution. This provides the parameters for the
%              convoluted distribution.
%
%INPUTS: mu1,Sigma1 The mean and variance of the first scalar normal
%                   distribution.
%        mu2,Sigma2 The mean and variance of the second scalar normal
%                   distribution.
%
%OUTPUTS: mu, SigmaInv The mean and variance of the product distribution.
%
%Note that the product distribution is normalized. The derivation of the
%product distribution is straightforward nothing that the Fourier transform
%of a normal distribution is also a fully-normalized (integrates to one) 
%normal distribution.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C. 

mu=mu1+mu2;
Sigma=Sigma1+Sigma2;

end

function [momentVal,coeffPolyPart]=momentGenFun(mu,Sigma,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the multivariate normal distribution. Taking
%              the ith, jth, kth... derivative of the moment generating
%              function with respect to the first, second, third...
%              components of the argument and evaluating it at t=0 provides
%              the noncentral moment of the multivariate normal
%              distribution involving the ith, jth, kth power of the first
%              second, third... components of the random vector.
%
%INPUTS: mu The mean of the PDF. If the PDF is multivariate, then this is a
%           numDimX1 column vector.
%     Sigma The variance (if scalar) or numDimXnumDim covariance matrix 
%           (if multidimensional) of the PDF.
% numDerivs A numDimX1 or 1XnumDim vector indicating the number of
%           derivatives to take with respect to each of the components of
%           the argument of the moment generating function. numDerivs>=0.
%         t The numDimXnumPoints argument of the moment generating
%           function at which the derivatives of the moment generating
%           function should be evaluated. If this parameter is omitted or
%           an empty matrix is passed, then a default of
%           t=zeros(numDim,1) is used.
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%     coeffPolyPart A hypermatrix taking numDim indices that can be
%                   evaluated using the polyValMultiDim at different values
%                   of t to get the polynomial coefficient that multiplies
%                   the exponential term in the given set of derivatives
%                   of the moment generating function.
%
%The moment generating function of a random vector is defined to be
%E(exp(t'*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter having the same dimensionality as the
%random variable. It can be shown that the moment generating function of a
%multivariate normal distribution is 
%E(exp(t'*x))=exp(t'*mu+(1/2)*t'*Sigma*t)
%Derivatives of this can be evaluated systematically. First,
%differentiating the exponential with respect to the ith component of t
%leads to the original exponential term times
%mu(i)+sum(Sigma(:,i).*t(:))
%Thus, all derivatives include the original exponential term times a
%multivariate polynomial. This function keeps track of the polynomial
%through all of the derivatives, using the chain rule, the fact that
%multivariate polynomial multiplication can be performed using the convn
%function, and using polyDerMultiDim and polySumMultiDim for
%multidimensional polynomial differentiation and addition.
%
%As an example, consider finding the noncentral third moment E(x1*x2*x3) of
%a 3D multivariate normal distribution. This can be done using
%If one solves for it by hand, one gets
%mu(1)*Sigma(2,3)+mu(2)*Sigma(1,3)+mu(3)*Sigma(1,2)+prod(mu)
%However, the moment can also be found by evaluating the derivatives of the
%moment generating function with respect to the first, second and third
%variables at t=[0;0;0]. For example, consider
% mu=[1;2;3];
% Sigma=[8, 3, 2;
%        3,11, 9;
%        2, 9,18];
% numDerivs=[1;1;1];
% momentVal=GaussianD.momentGenFun(mu,Sigma,numDerivs)
%One will find the value of the noncentral moment to be 28, which is
%correct.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

%The number of dimensions.
numIdx=length(mu);

%Derivative of the moment generating function produce a polynomial times
%the original exponential term. We shall keep track of that polynomial
%while differentiating. The loops below take successive derivatives with
%respect to the different dimensions of the argument vector of the moment
%generating function. Each component is differentiated for the required
%number of derivatives. coeffPolyPart holds the accumulated polynomial that
%multiplies the exponential part of the differentiated moment generating
%function. Initially, coeffPolyPart=1 to signify that no derivatives have
%been taken. It is initially implemented with more elements than necessary
%to simplify the addition of the terms after the use of the polyDerMultiDim
%function below.
coeffPolyPart=1;

%i is the current dimension we are differentiating.
for i=1:numIdx
    %basePoly shall be the multivariate polynomial that is added every time
    %the exponential term in the moment generating function is
    %differentiated with respect to the ith index. This makes basePoly the
    %derivative with respect to the ith element of t of the argument of the
    %exponent in the multivariate normal moment generating function
    
    %Allocate space for the polynomial and make it have the correct shape.
    if(numIdx>1)
        basePoly=reshape(zeros(2^numIdx,1),2*ones(1,numIdx));
    else
        basePoly=zeros(2^numIdx,1);
    end
    
    %The additive term is the component of the mean that was multiplied by 
    basePoly(1)=mu(i);
    
    %The term multiplied by x consists of elements from the covariance
    %matrix. This just consists of terms in the ith row of the
    %covariance matrix (the matrix is symmetric).
    idxVec=ones(numIdx,1);
    %The dimensionalities of all of the variables for the nDim2Index
    %function.
    dims=2*ones(numIdx,1);
    for curDim=1:numIdx
        idxVec(curDim)=2;
        basePoly(nDim2Index(dims,idxVec))=Sigma(i,curDim);
        idxVec(curDim)=1;
    end
    
    %Now, we enter into a loop to evaluate derivatives of the moment
    %generating function with respect to the current dimension.
    for derivsLeft=numDerivs(i):-1:1
        %The chain rule means that there are two terms to consider (both of
        %which are multiplied by the same exponential term).
        %The first term is the current coeffPolyPart times the derivative 
        %of the exponential term with respect to the ith dimensions. This
        %is the product of coeffPolyPart and basePoly. The second term is
        %the exponential term times the derivative of coeffPolyPart. The
        %two terms then must be added.
        term1=convn(basePoly,coeffPolyPart);%Multiply the polynomials.
        term2=polyDerMultiDim(coeffPolyPart,i);
        
        %Add the multivariate polynomials.
        coeffPolyPart=polySumMultiDim(term1,term2);
    end
end

%coeffPolyPart now contains the multivariate polynomial that is multiplied
%by the exponential term. If no value of t is given, then just evaluate it
%at t=0. In this instance, the exponential term is zero and only the
%constant term from the polynomial appears.
if(nargin<4||isempty(t))
    momentVal=coeffPolyPart(1);
else
    numPoints=size(t,2);
    momentVal=zeros(numPoints,1);
    for curPoint=1:numPoints
        tCur=t(:,curPoint);
        momentVal(curPoint)=exp(tCur'*mu+tCur'*Sigma*tCur)*polyValMultiDim(coeffPolyPart,tCur);
    end
end

end

function cumVal=cumGenFun(mu,Sigma,numDerivs,t)
%%CUMGENFUN Evaluate the cumulant generating function (or one of its
%           derivatives) of the multivariate normal distribution. Taking
%           the ith, jth, kth... derivative of the cumulant generating
%           function with respect to the first, second, third...
%           components of the argument and evaluating it at t=0 provides
%           the cumulant of the multivariate normal distribution involving
%           the ith, jth, kth power of the first second, third...
%           components of the random vector. The cumulant generating
%           function is the natural logarithm of the moment generating
%           function.
%
%INPUTS: mu The mean of the PDF. If the PDF is multivariate, then this
%           is a numDimX1 column vector.
%     Sigma The variance (if scalar) or numDimXnumDim covariance matrix 
%           (if multidimensional) of the PDF.
% numDerivs A numDimX1 or 1XnumDim vector indicating the number of
%           derivatives to take with respect to each of the components of
%           the argument of the cumulant generating function.
%           numDerivs>=0.
%         t The numDimXnumPoints argument of the cumulant generating
%           function at which the derivatives of the cumulant generating
%           function should be evaluated. If this parameter is omitted or
%           an empty matrix is passed, then a default of
%           t=zeros(numDim,1) is used.
%
%OUTPUTS: cumVal A numPointsX1 vector of the values of the derivatives
%                of the cumulant generating function given at the points
%                in t or at t=0 if t is omitted.
%
%Cumulants are useful in interpolating probability distributions through
%the use of, for example, an Edgeworth series. The cumulant generating
%function is defined as the natural logarithm of the moment generating
%function. It can be shown that the moment generating function of the
%multivariate normal distribution is 
%E(exp(t'*x))=exp(t'*mu+(1/2)*t'*Sigma*t)
%Thus, the cumulant generating function is just
%t'*mu+(1/2)*t'*Sigma*t
%Consequently, for derivatives higher than two, the cumulant generating
%function (and hence the cumulants) of the multivariate normal distribution
%are all zero.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(t))
    numDim=length(mu);
    t=zeros(numDim,1);
end

numPoints=size(t,2);
cumVal=zeros(numPoints,1);

sumVal=sum(numDerivs);
switch(sumVal)
    case 0%No derivatives.
        for curPoint=1:numPoints
            tCur=t(:,curPoint);
            cumVal(curPoint)=tCur'*mu+(1/2)*tCur'*Sigma*tCur;
        end
    case 1
        %Find the nonzero term.
        derivIdx=find(numDerivs);
        for curPoint=1:numPoints
            tCur=t(:,curPoint);
            cumVal(curPoint)=mu(derivIdx)+sum(Sigma(:,derivIdx).*tCur(:));
        end
    case 2
        derivIdx=find(numDerivs);
        cumVal(:)=Sigma(derivIdx(1),derivIdx(2));%The same for all t.
    otherwise%No third or higher order cumulants.
       cumVal(:)=0; 
end
end

function x=rand(N,mu,P)
%%RAND Generate multivariate Gaussian random variables with a given mean
%      vector and covariance matrix.
%
%INPUTS: N The number of random variables to generate.
%       mu The xDim X1 mean of the multivariate Gaussian to generate.
%        P The xDim X xDim positive definite covariance matrix of the
%          multivariate Gaussian to generate. If this parameter is omitted
%          or an empty matrix is passed, then the identity matrix will be
%          used.
%
%OUTPUTS: x An xDimXN matrix of random instances of the multivariate
%           Gaussian distribution.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDim=size(mu,1);
    if(nargin<3||isempty(P))
        P=zeros(xDim,xDim);
    end

    xDim=size(mu,1);
    x=bsxfun(@plus,mu,chol(P,'lower')*randn(xDim,N));
end

function x=randS(N,mu,S)
%%RANDS Generate multivariate Gaussian random variables with a given mean
%       vector and lower-triangular square root covariance matrix.
%
%INPUTS: N The number of random variables to generate.
%       mu The xDim X1 mean of the multivariate Gaussian to generate.
%        S The xDim X xDim lower triangular square root covariance matrix
%          of the multivariate Gaussian to generate.
%
%OUTPUTS: x An xDimXN matrix of random instances of the multivariate
%           Gaussian distribution.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDim=size(mu,1);
    x=bsxfun(@plus,mu,S*randn(xDim,N));
end

function [P,error,curIter]=integralOverRegion(mu, Sigma,minVals,maxVals,epsVal,alpha,maxIter)
%%INTEGRALOPVERREGION Compute the probability of a Gaussian probability
%                  density function (PDF) within a (hyper-)rectangular
%                  region. The probability is computed using a transformed
%                  Monte Carlo method designed for this specific integral
%                  that converges significantly faster than generic,
%                  textbook Monte Carlo integration techniques.
%
%INPUTS:  mu The NX1 mean vector of the multivariate normal distribution.
%      Sigma The NXN positive definite covariance matrix of the
%            distribution.
%    minVals An NX1 or 1XN vector of the lower integration bounds for all
%            of the dimensions.
%    maxVals An NX1 or 1XN vector of the upper integration bounds. Note
%            that minVals(i)<maxVals(i) for all elements.
%     epsVal The desired error tolerance for the probability computation.
%            If omitted or an empty matrix is passed, the defaut of 1e-6
%            (1e-4%) is used.
%      alpha The confidence factor for the standard error tolerance. Thus,
%            for the estimate of epsVal used to determine termination to be
%            correct 99% of the time, one should use 
%            alpha=GaussianD.invCDF(0.99); If this parameter is omitted or
%            an empty matrix is passed, the default value of 3.5 is used.
%    maxIter The maximum number of iterations to use. If this parameter is
%            omitted or an empty matrix is passed, the default value of
%            1000 is used.
%
%OUTPUTS: P The approximate probability within the region (0-1).
%     error The estimated error in P, with a confidence interval determined
%           by alpha.
%   curIter The number of the last iteration performed before termination.
%
%The algorithm is that of [1]. Various transformations are applied to the
%function to achieve better convergence than when performing typical Monte
%Carlo integration.
%
%EXAMPLE:
%This is the 3D example in the paper. The probability should be about
%0.8279.
% minBounds=[-Inf;-Inf;-Inf];
% maxBounds=[1;4;2];
% R=[1,3/5,1/3;
%    3/5,1,11/15;
%    1/3,11/15,1];
% mu=[0;0;0];
% PApprox=GaussianD.integralOverRegion(mu,R,minBounds,maxBounds)
%
%REFERENCES:
%[1] A.Genz, "Numerical computation of multivariate normal probabilities,"
%    Journal of Computational and Graphical Statistics, vol. 1, no. 2, pp.
%    141-149, June 1992.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    if(nargin<7||isempty(maxIter))
       maxIter=1000; 
    end
    
    if(nargin<6||isempty(alpha))
       alpha=3.5; 
    end
    
    if(nargin<5||isempty(epsVal))
       epsVal=1e-6; 
    end
    
    %Recenter things for a zero-mean distribution.
    minVals=minVals-mu;
    maxVals=maxVals-mu;
    
    %First, we reorder the variables so that the largest integration
    %regions come last. This is suggested in [1] to make the algorithm
    %faster.
    [~,idxSort]=sort(maxVals-minVals,'ascend');
    a=minVals(idxSort);
    b=maxVals(idxSort);
    C2=Sigma(idxSort,idxSort);

    C=chol(C2,'lower');

    numDim=size(C,1);
    
    intSum=0;
    varSum=0;
    
    %Allocate space
    d=zeros(numDim,1);
    e=zeros(numDim,1);
    f=zeros(numDim,1);
    
    d(1)=GaussianD.CDF(a(1)/C(1,1));
    e(1)=GaussianD.CDF(b(1)/C(1,1));
    f(1)=e(1)-d(1);
    
    y=zeros(numDim-1,1);%Allocate space.
    for curIter=0:(maxIter-1)
        
        w=rand(numDim-1,1);
        
        for i=2:numDim
            y(i-1)=GaussianD.invCDF(d(i-1)+w(i-1)*(e(i-1)-d(i-1)));
            deltaVal=C(i,1:(i-1))*y(1:(i-1));
            
            if(a(i)==-Inf)
                d(i)=0;
            else
                d(i)=GaussianD.CDF((a(i)-deltaVal)/C(i,i));
            end
            
            if(b(i)==Inf)
                e(i)=1;
            else
                e(i)=GaussianD.CDF((b(i)-deltaVal)/C(i,i));
            end
            f(i)=(e(i)-d(i))*f(i-1);
        end
        
        N=curIter+1;

        delta=(f(numDim)-intSum)/N;
        intSum=intSum+delta;
        varSum=(N-2)*varSum/N+delta^2;
        error=alpha*sqrt(varSum);
        
        if(error<epsVal)
            break;
        end
    end
    
    %The probability estimate.
    P=min(1,max(0,intSum));
end

function entropyVal=entropy(Sigma)
%%ENTROPY Obtain the differential entropy of the multivariate Gaussian
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: Sigma The NXN positive definite covariance matrix of the
%              distribution.
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

    entropyVal=(1/2)*log(det(2*pi*exp(1)*Sigma));
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
