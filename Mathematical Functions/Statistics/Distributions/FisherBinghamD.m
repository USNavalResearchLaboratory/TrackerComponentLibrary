classdef FisherBinghamD
%%FISHERBINGHAMD Functions to handle the Fisher-Bingham distribution which
%                is a common, general directional distribution on the
%                n-sphere. The 4D distribution is often used for
%                orientation estimation using quaternions.
%Implemented methods are: PDF, (only for the Bingham distribution)
%                         covApprox, normProdDist, Gauss2FisherBingham,
%                         normConstApprox, (only for the Bingham
%                         distribution) normConstDerivApprox
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    function val=PDF(x,A,g,normConst)
    %%PDF Evaluate the probability density function (PDF) of the
    %     Fisher-Bingham distribution at a given set of points. The
    %     Fisher-Bingham distribution comes from constainting an N.2
    %     dimensional Gaussian distribution to the unit sphere. However, it
    %     is normally parameterized slightly differently than a typical
    %     multivariate Gaussian and its parameters can take more general
    %     values. The PDF is f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x).
    %     One can find the A and g for a particular Gaussian
    %     distribution using the function Gauss2FisherBingham.
    %
    %INPUTS: x An xDimXN set of N xDim-dimensional points on the unit
    %          sphere (magnitude must be one) at which the PDF should be
    %          evaluated. xDim>=2.
    %        A A type of inverse covariance matrix for the distribution. It
    %          must be symmetric. It does not have to be positive definite.
    %        g A real mean parameter, which need not have unit magnitude.
    %          If this parameter is omitted or an empty matrix is passed,
    %          then zeros will be used.
    %normConst Optionally the normalizing constant of the distribution can 
    %          be provided. If omitted or an empty matrix is passed,
    %          the function FisherBinghamD.normConstApprox is
    %          used to approximate it. This is offered as an option,
    %          because computing the normalizing constant can be difficult
    %          in general scenarios.
    %
    %OUTPUTS: val The value(s) of the Fisher-Bingham distribution with the
    %             specified parameters evaluated at x.
    %
    %A common formulation of the Fisher-Bingham distribution is given in
    %[1] and is the form used here. Note that definitions of A and g can
    %vary in the literature.
    %
    %For very peaky distributions, e.g. A=[1e20,0;0,1e20], the
    %normalization constant will evaluate to be numerically zero and the
    %PDF will return NaNs, except when at the mean, when Inf will be
    %returned.
    %
    %REFERENCES:
    %[1] A. Kume and A. T. A. Wood, "Saddlepoint approximations for the
    %    Bingham Fisher-Bingham normalizing constants," Biometrika, vol.
    %    92, no. 2, pp. 465-476, Jun. 2005.
    %   
    %August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
        
        xDim=size(x,1);
        numX=size(x,2);
        
        if(nargin<3||isempty(g))
           g=zeros(xDim,1); 
        end
        
        if(nargin<4||isempty(normConst))
            normConst=FisherBinghamD.normConstApprox(A,g);
        end

        val = zeros(numX,1);
        for curX=1:numX
            xCur = x(:,curX);
            val(curX)=(1/normConst)*exp(-xCur'*A*xCur+g'*xCur);
        end
    end
    
    function P=covApprox(A)
    %%COVAPPROX  Obtain the approximate linear covariance matrix of the
    %            Bingham distribution (the Bingham distribution is a
    %            special case of the Fisher-Bingham distribution with g=0).
    %            Note that this is the traditional definition of the
    %            covarinace, not one that is specially tuned for parameters
    %            constrained to the unit (hyper-)sphere. The solution is
    %            approximate as it is based on approximations for the
    %            normalizing constant of the distribution and its
    %            derivatives.
    %
    %INPUTS: A A type of inverse covariance matrix for the distribution. It
    %          must be symmetric. It does not have to be positive definite.
    %
    %OUTPUTS: P The covariance matrix for the distribution.
    %
    %The expression for the covariance matrix is taken from Section II of
    %[1].
    %
    %The PDF for the Bingham distribution is
    %f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x).
    %with g=0.
    %
    %REFERENCES:
    %[1] I. Gilitschenski, G. Kurz, S. J. Julier, and U. D. Hanebeck,
    %    "Unscented orientation estimation based on the Bingham
    %    distribution," IEEE Transactions on Automatic Control, accepted 2015.
    %   
    %September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

        numDim=size(A,1);
        
        [M,Z]=eig(A);
        
        normConst=FisherBinghamD.normConstApprox(Z);
        
        numDerivVec=zeros(numDim,1);
        derivVec=zeros(numDim,1);
        
        for curDim=1:numDim
            numDerivVec(curDim)=1;
            derivVec(curDim)=FisherBinghamD.normConstDerivApprox(Z,numDerivVec);
            numDerivVec(curDim)=0;
        end
        
        d=derivVec./normConst;
        %The diagonal is always 1, which, due to the use of approximations,
        %might not be the case here. Thus, we force it to be 1.
        d=d/sum(d);
        
        %The -1, which is not in [1], deals with the definition of the
        %Bingham distribution lacking the leading - term, which is assumed
        %here.
        P=M*diag(d)*M';
        %Force it to be symmetric.
        P=(P+P')/2;
    end
   
    function [A,g,normConst]=Gauss2FisherBingham(mu,SigmaInv)
    %%GAUSS2FISHERBINGHAM The Fisher-Bingham distribution comes from
    %           conditioning a multivariate normal distribution to the unit
    %           hypersphere. This function converts the mean and inverse
    %           covariance matrix of a multivariate normal distribution
    %           into the A and g parameters of the Fisher-Bingham
    %           distribution. Note that this is not the same as the
    %           distribution obtained by taking a normal random vector and
    %           dividing by its magnitude (projecting it onto the unit
    %           sphere).
    %
    %INPUTS: mu The xDimX1 mean of the normal distribution with xDim>=2.
    %  SigmaInv The xDimXxDim inverse of the covariance matrix of the
    %           normal distribution.
    %
    %OUTPUTS: A, g The covariance and mean parameters used in the
    %                  Fisher-Bingham distribution. For example, in
    %                  FisherBinghamD.PDF.
    %        normConst An approximation of the normalizing constant of the
    %                  distribution.
    %
    %The conversion just comes from looking for equivalent parameters
    %between a multivariate Gaussian PDF and the Fisher-Bingham
    %distribution.
    %
    %August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

        A=(1/2)*SigmaInv;
        g=SigmaInv*mu;
        if(nargout>2)
            normConst=FisherBinghamD.normConstApprox(A,g);
        end
    end
    
    
    function [A,g,normConst]=normProdDist(A1,g1,A2,g2)
    %%NORMPRODDIST The product of two Fisher-Bingham distributions is an
    %              unnormalized Fisher-Bingham distribution. This finds the
    %              parameters of the normalzied product distribution.
    %
    %INPUTS: A1, g1 The parameters of the first Fisher-Bingham
    %                   distribution (A is symmetric).
    %        A2, g2 The parameters of the first Fisher-Bingham
    %                   distribution (A is symmetric).
    %
    %OUTPUTS: A, g The parameters of the normalized product
    %                  distribution.
    %        normConst The normalizing constant of the normalized product
    %                  distribution.
    %
    %The product distribution is very simple to find by multiplying both of
    %the PDFs of the form f(x|Ai,gi)=(1/c)*exp(-x'*Ai*x+gi'*x)
    %
    %August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

        A=A1+A2;
        g=g1+g2;
        if(nargout>2)
            normConst=FisherBinghamD.normConstApprox(A,g);
        end
    end
    
    function c=normConstApprox(A,g)
    %%NORMCONSTAPPROX Approximate the normalization constant
    %               of the Fisher-Bingham distribution (the Bingham
    %               distribution is a special case of the Fisher-Bingham
    %               distribution with g=0). The Fisher-Bingham
    %               distribution is what one obtains when a multivariate
    %               normal random vector is conditioned to have unit
    %               magnitude. f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x)
    %               where c is the normalization constant found by this
    %               function, x is the random variable with unit magnitude
    %               and A and g are parameters.
    %
    %INPUTS: A A type of inverse covariance matrix for the distribution. It
    %          must be symmetric. It does not have to be positive definite.
    %          For the normalization constant, only the eigenvalues of A
    %          matter.
    %        g A real mean parameter, which need not have unit magnitude.
    %          If this parameter is omitted or an empty matrix is passed,
    %          a zero vector shall be used for g.
    %
    %OUTPUTS: c An approximation to the normalizing constant of the PDF
    %           f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x)
    %
    %The second-order saddle point approximation of [1] is used. Numerically
    %computing the exact value for arbitrary A and g is actually quite
    %difficult and requires the evaluation of matrix hypergeometric
    %functions. When g=0, Newton's method is proven to be optimal in [2],
    %as the convexivity of the cost function is proven. Here, we assume
    %that the cost function is convex for a general g and use Newton's
    %method as well.
    %
    %Note that the distribution can be rewritten as
    %f(x|A,g)=(1/cp)*exp(-(1/2)*(x-gp)'*Ap*(x-gp))
    %after adjusting A and g and changing the normalziing constant. This
    %formulation more closely reflects the origin of the Fisher-Bingham
    %distribution as being a normal distribution consitioned on x having
    %unit length.
    %
    %REFERENCES:
    %[1] A. Kume and A. T. A. Wood, "Saddlepoint approximations for the
    %    Bingham and Fisher-Bingham normalizing constants," Biometrika,
    %    vol. 92, no. 2, pp. 465-476, Jun. 2005.
    %[2] I. Gilitschenski, G. Kurz, S. Julier, and U. D. Hanebeck,
    %    "Efficient Bingham filtering based on saddlepoint approximations,"
    %    in Proceedings of the International Conference on Multisensor
    %    Fusion and Information Integration for Intelligent Systems,
    %    Beijing, China, 28-29 Sep. 2014.
    %
    %August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
        if(nargin<2||isempty(g))
           numDim=size(A,1);
           g=zeros(numDim,1);
        end
    
        maxIter=1000;
        cumFuncPrec=1e-14;

        lambda=eig(A);
        p=length(lambda);

        %The saddlepoint approximation is for all positive values of the
        %eigenvalues of A. However, as described after Equation 16 in [1],
        %we can offset the values and then scale the distribution
        %accordingly.
        deltaLambda=min(lambda)-0.1;
        lambda=lambda-deltaLambda;%min(lambda) is now 0.1
        %The final result is scaled by this factor to undo the change in
        %lambda.
        scalFact=exp(-deltaLambda);

        %The 0.1 is the minimum value of lambda.
        upperBound=0.1-(1/4)-(1/2)*sqrt(1/4+min(g.^2));

        %If the noncentrality parameter is zero, use Newton's method to
        %find t as it was shown in [2] that the function is optimal as the
        %cost function was proven to be convex. here, we assume that the
        %cost function is also convex with g~=0.
        t=upperBound-0.5;%Initial estimate for Newton's method.
        for curIter=1:maxIter
            %One step of Newton's method
            cumFunc=KFunc(1,t,lambda,g);

            t=t-(cumFunc-1)/KFunc(2,t,lambda,g);
            if(abs((cumFunc-1))<=cumFuncPrec)
                break;%If it has converged.
            end
        end

        k2=KFunc(2,t,lambda,g);
        p3=KFunc(3,t,lambda,g)/k2^(3/2);
        p4=KFunc(4,t,lambda,g)/k2^(4/2);

        T=(1/8)*p4-(5/24)*p3^2;%Equation 13 in [1].

        %Equation 15 for the first-order saddlepoint approximation.
        c1=sqrt(2)*pi^((p-1)/2)*(1/sqrt(k2))*prod(1./sqrt(lambda-t))*exp(-t+(1/4)*sum(g.^2./(lambda-t)));
        %The modified second order saddlepoint approxiamtion from Equation
        %16
        c3=c1*exp(T);

        %The scaling factor undoes the initial adjustments to lambda.
        c=scalFact*c3;

        function val=KFunc(j,t,lambda,g)
            %The jth derivative of the cumulant generating function of a
            %linear combination of independent noncentral chi square
            %distributions (j>=0).
            if(j==0)
                %The cumulant generating function (before Equation 9).
                val=sum(-(1/2)*log(1-t./lambda)+(1/4)*g.^2./(lambda-t)-g.^2./(4*lambda));
            elseif(j==1)
                %The first derivative of the cumulant generating function.
                %Equation 9 in [1].
                val=sum((1/2)*(1./(lambda-t))+(1/4)*(g.^2./(lambda-t).^2));
            else
                %The jth derivative (j>=2) of the cumulant generating
                %function. This is Equation 10 in [1]
                val=sum(factorial(j-1)./(2*(lambda-t).^j)+(factorial(j)/4)*(g.^2./(lambda-t).^(j+1)));
            end
        end
    end

    function cDeriv=normConstDerivApprox(Z,numDeriv)
    %NORMCONSTDERIVAPPROX Approximate an arbitrary derivative of the
    %               normalization constant of the Bingham distribution (the
    %               Bingham distribution is a special case of the Fisher-
    %               Bingham distribution with g=0). The derivative is taken
    %               with respect to an eigenvalue of the inverse covariance
    %               matrix parameterizing the distribution. This function
    %               takes the eigenvalues rather than the inverse
    %               covariance matrix itself so as to eliminate
    %               uncertainty in which derivative is with repect to which
    %               eigenvalue. The derivatives of the normalization
    %               constant play a role in parameter estimation.
    %
    %INPUTS: Z The eigenvalues of the inverse covariance matrix in the
    %          Bingham distribution. This is the eigenvalues of A in
    %          the PDF f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x), where g=0, because
    %          this is the Bingham distribution and not the Fisher-Bingham
    %          distribution. This can be a vector or a diagonal matrix.
    %          Given A, one can get Z using [M,Z]=eig(A);
    % numDeriv If Z contains numDim eigenvalues, numDeriva is a numDimX1 or
    %          1XnumDim vector where numDeriv(k) is the non-negative
    %          integer number of derivatives that should be taken with
    %          respect to the kth eigenvalue in Z.
    %
    %The formula for derivatives of the normalization constant with respect
    %with respect to the eigenvalues of A expressed as a scaling of the
    %normalization constant itself, is given in [1] and is implemented
    %here. The implementation is an approxiamtion, because the function
    %FisherBinghamD.normConstApprox is used as opposed to using the exact,
    %true normaliztion constant. An application of how these derivatives
    %can be used with parameter estimation methods is given in [2]. Note
    %that the definition of the distribution in [2] has absorbed the minus
    %sign in the exponent into the Z matrix.
    %
    %REFERENCES:
    %[1] A. Kume and A. T. A. Wood, "On the derivatives of the normalizing
    %    constant of the Bingham distribution," Statistics and Probability
    %    Letters, vol. 77, no. 8, pp. 832-837, 15 Apr. 2007.
    %[2] I. Gilitschenski, G. Kurz, S. Julier, and U. D. Hanebeck,
    %    "Efficient Bingham filtering based on saddlepoint approximations,"
    %    in Proceedings of the International Conference on Multisensor
    %    Fusion and Information Integration for Intelligent Systems,
    %    Beijing, China, 28-29 Sep. 2014.
    %
    %August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
        
        %Equation 6 in [1]
        sumVal=sum(numDeriv);
        a=(-1)^sumVal*pi^(-sumVal)*prod(gamma((1+2*numDeriv)/2)/gamma(1/2));
    
        %If the input is a diagonal matrix.
        if(size(Z,1)==size(Z,2))
            lambda=diag(Z);
        else
            lambda=Z(:);
        end
        
        ZTilde=diag(runLenDecode(lambda,2*numDeriv(:)+1));
        cDeriv=a*FisherBinghamD.normConstApprox(ZTilde);
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
