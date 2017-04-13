classdef StudentTD
%%STUDENTD Functions to handle the multivariate Student-t distribution.
%    Note that the multivariate Student-t distribution with one degree of
%    freedom is the same as the multivariate Cauchy distribution.
%    Also the scalar student-t distribution is the distribution of
%    x/sqrt(y/n) where x is a normal random variable, y is a central
%    chi squared random variable and n is the number of degrees of freedom
%    of y. In terms of tracking, the distribution can be of interest,
%    because, as noted in [1], angular measurements in the presence of
%    glint have been fit to Cauchy distributions in some instances.
%Implemented methods are: mean, cov, PDF, rand (only for scalar, central,
%                                               non-scaled distributions).
%
%REFERENCES:
%[1] U. Nickel, "Angular superresolution with phased array radar: A review
%    of algorithms and operational constraints," IEEE Proceedings, vol.
%    134, Part F, no. 1, pp. 53-59, Feb. 1987.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(mu,nu)
%%MEAN  Obtain the mean of the multivariate Student-t distribution for 
%       given location and scale parameters and degrees of freedom.
%
%INPUTS: mu The DX1 location vector of the student's-t distribution.
%        nu The scalar number of degrees of freedom of the Student-t
%           distribution. nu>=0.
%
%OUTPUTS: val The mean of the multivariate Student-t distribution.
%
%The mean is undefined if nu<=1.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
    if(nu>1)
        val=mu;
    else
        val=NaN;
    end
end 

function val=cov(Sigma,nu)
%%MEAN  Obtain the covariance matrix of the multivariate Student-t 
%       distribution for given location and scale parameters and degrees 
%       of freedom.
%
%INPUTS: Sigma The DXD symmetric, positive definite scale matrix.
%           nu The scalar number of degrees of freedom of the Student-t
%              distribution. nu>=0.
%
%OUTPUTS: val The covariance matrix of the multivariate Student-t
%             distribution.
%
%The covariance matrix is undefined if nu<=2.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
    if(nu>2)
        val=nu/(nu-2)*Sigma;
    else
        val=NaN;
    end
end 
    
function val=PDF(x,mu,Sigma,nu)
%%PDF Evaluate the monovariate or multivariate Student-t distribution at
%     the desired point.
%
%INPUTS: x The DX1 vector at which the (possibly multivariate) Student's-t
%          distribution is to be evaluated.
%       mu The DX1 location vector of the student's-t distribution.
%    Sigma The DXD symmetric, positive definite scale matrix.
%       nu The scalar number of degrees of freedom of the Student-t
%          distribution. nu>=0.
%
%OUTPUTS: val The PDF of the Student's-t distribution at x with the given
%             parameters.
%
%The vector version of the Student-t Distribution is given in Appendix B of
%[1]. As nu->Inf, the distribution reduces to a multivariate Gaussian
%distribution with mean mu and covariance matrix Sigma.
%
%Logarithms are used in the implementation to reduce the effect of
%precision problems that can arise if nu is large. The problems arise due 
%to the ratio of gamma functions.
%
%REFERENCES:
%[1] C. M. Bishop, Pattern Recognition and Machine Learning. Cambridge,
%    United Kingdom: Springer, 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    D=size(x,1);

    diff=x-mu;
    Delta2=diff'*pinv(Sigma)*diff;

    num=gammaln((nu+D)/2);
    denom=gammaln(nu/2)+log(sqrt(det(nu*pi*Sigma)))+((D+nu)/2)*log(1+Delta2/nu);

    val=num-denom;
    val=exp(val);
end

function vals=rand(N,nu)
%%RAND Generate scalar Student-t random variables with unit scale factor
%      and a given number of degrees of freedom.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       nu The scalar number of degrees of freedom of the Student-t
%          distribution. nu>=0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated scalar Student-t random variables.
%
%The algorithm implemented is the TIR algorithm in [1].
%
%REFERENCES:
%[1] A. J. Kinderman, J. F. Monahan, and J. G. Ramage, "Computer methods
%    for sampling from student?s t distribution," Mathematics of
%    Computation, vol. 31, no. 140, pp. 1009-1018, Oct. 1977.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(isscalar(N))
    dims=[N, N];
else
    dims=N;
end

vals=zeros(dims);
numVals=numel(vals);

b=sqrt(2*exp(-1/2)-1);
alpha=nu;
for curVal=1:numVals
    while(1)        
        %Step 1
        u=rand(1);

        if(u<b/2)
            x=4*u-b;
            %Step 2
            v=rand(1);
            if(v<=1-abs(x)/2)
               vals(curVal)=x;
               break;
            end

            uAlpha=(1+x^2/alpha)^(-(alpha+1)/2);

            if(v<=uAlpha)
                vals(curVal)=x;
                break;
            end
            continue;
        end
    
        if(u<0.5)
            %Step 3
            temp=4*u-1-b;
            x=(abs(temp)+b)*sign(temp);
            v=rand(1);

            %Step 4
            if(v<=1-abs(x)/2)
                vals(curVal)=x;
                break;
            end

            if(v>=(1+b^2)/(1+x^2))
                continue;
            end

            uAlpha=(1+x^2/alpha)^(-(alpha+1)/2);
            if(v<=uAlpha)
                vals(curVal)=x;
                break;
            end
            continue;
        end
    
        if(u<0.75)
            %Step 5
            temp=8*u-5;
            x=2/((abs(temp)+1)*sign(temp));
            u1=rand(1);
            v=x^(-2)*u1;

            %Step 4 again
            if(v<=1-abs(x)/2)
                vals(curVal)=x;
                break;
            end

            if(v>=(1+b^2)/(1+x^2))
                continue;
            end

            uAlpha=(1+x^2/alpha)^(-(alpha+1)/2);
            if(v<=uAlpha)
                vals(curVal)=x;
                break;
            end
            continue;
        end

        %Step 6
        x=2/(8*u-7);
        v=rand(1);

        uAlpha=(1+x^2/alpha)^(-(alpha+1)/2);
        if(v<x^2*uAlpha)
            vals(curVal)=x;
            break;
        end
    end
end
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
