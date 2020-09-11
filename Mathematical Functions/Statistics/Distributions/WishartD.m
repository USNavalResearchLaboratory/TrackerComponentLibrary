classdef WishartD
%%WISHARTD Functions to handle the Wishart distribution.
%Implemented methods are: mean, var, cov, covTerm, secondMatMoment,
%                         secondMoments, secondMomentTerm, PDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(A,nu)
%%MEAN Find the mean of a particular Wishart distribution
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%
%OUTPUTS: val The dXd mean of the Wishart distribution.
%
%The Wishart distribution and its generation are discussed in chapters 5
%and 8 of [1].
%
%EXAMPLE:
%Here, we demonstrate that the sample mean is consistent with the mean
%returned by this function.
% A=[55,  7, 12, 17;
%     7, 45, 17, 22;
%    12, 17, 35, 27;
%    17, 22, 27, 25];
% nu=6;
% numMC=1e5;
% XMean=zeros(4,4);
% for k=1:numMC
%     XMean=XMean+WishartD.rand(A,nu);
% end
% XMean=XMean/numMC;
% XMeanTrue=WishartD.mean(A,nu);
% RelErr=max(abs((XMean(:)-XMeanTrue(:))./XMeanTrue(:)))
%The relative error will be on the order of 0.005, which is good.
%Increasing the number of Monte Carlo runs will decrease the relative
%error.
%
%REFERENCES:
%[1] M. L. Eaton, Multivariate Statistics: A Vector Space Approach, ser.
%    Institute of Mathematical Statistics Lecture Notes-Monograph Series.
%    Beachwood, Ohio: Institute of Mathematical Statistics, 2007, vol. 53.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=A*nu;
end

function val=var(A,nu)
%%VAR Obtain the variance of the Wishart distribution. Though the random
%     variable is a matrix, the variance is E{(X-xMean)*(X-XMean)}, where
%     xMean is the mean of the distribution. Note that this is different
%     than the definition used for the method cov.
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%
%OUTPUTS: val The dXd variance of the Wishart distribution.
%
%%An expression for the second moment is given in Equation 3.3.15 on page 98
%of [1]. Using the mean and the second moment, the variance is
%E{X^2}-E{X}^2.
%
%EXAMPLE:
%Here, we demonstrate that the sample variance is consistent with the
%variance returned by this function.
% A=[55,  7, 12, 17;
%     7, 45, 17, 22;
%    12, 17, 35, 27;
%    17, 22, 27, 25];
% nu=6;
% numMC=1e5;
% 
% XMean=zeros(4,4);
% XVals=zeros(4,4,numMC);
% for k=1:numMC
%     X=WishartD.rand(A,nu);
%     XMean=XMean+X;
%     XVals(:,:,k)=X;
% end
% XMean=XMean/numMC;
% 
% XVar=zeros(4,4);
% for k=1:numMC 
%     diff=XVals(:,:,k)-XMean;
% 
%     XVar=XVar+diff*diff;
% end
% XVar=XVar/numMC;%Sample variance
% XVarTrue=WishartD.var(A,nu);
% RelErr=max(abs((XVar(:)-XVarTrue(:))./XVarTrue(:)))
%The relative error will typically be from 0.01 to 0.03, indicating good
%agreement.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
Ex=nu*A;%The mean
Ex2=(nu+nu^2)*A*A+nu*trace(A)*A;

val=Ex2-Ex*Ex;
end

function matVal=cov(A,nu)
%%COV Calculate the covariance of all of the elements of the dXd Wishart
%     matrix with each other. Since the Wishart distribution is a matrix
%     distribution, the covariance is stored in a tensor (4D matrix).
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%
%OUTPUTS: matVal The dXdXdXd matrix of covariance values
%                matVal(i1,i2,i3,i4) is the covariance of X(i1,i2) with
%                X(i3,i4).
%
%An expression for the covariance is given on page 99 of [1].
%
%EXAMPLE:
%Here, we demonstrate that the covariance values in this function are
%consistent with those from the first and second moments, using the
%identity that Cov(X,Y)=E(X*Y)+E(X)*E(Y).
% A=[55,  7, 12, 17;
%     7, 45, 17, 22;
%    12, 17, 35, 27;
%    17, 22, 27, 25];
% nu=6;
% 
% d=size(A,1);
% meanVal=WishartD.mean(A,nu);
% prodVals=zeros(d,d,d,d);
% for i=1:d
%     for j=1:d
%         for k=1:d
%             for l=1:d
%                 prodVals(i,j,k,l)=meanVal(i,j)*meanVal(k,l);
%             end
%         end
%     end
% end
% covMat=WishartD.secondMoments(A,nu)-prodVals;
% max(abs(vec(covMat-WishartD.cov(A,nu))))
%One will see that the maximum difference between the covariance valuefrom
%this function and that computed from the first and second noncentral
%moments is zero.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(A,1);
    matVal=zeros(d,d,d,d);
    
    for i=1:d
        for j=i:d 
            for k=1:d
                for l=k:d
                    val=nu*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
                    matVal(i,j,k,l)=val;
                    matVal(j,i,k,l)=val;
                    matVal(i,j,l,k)=val;
                    matVal(j,i,l,k)=val;
                end
            end
        end
    end
end

function val=covTerm(ijkl,A,nu)
%%COVTERM Compute the covariance of X(i,j) and X(k,l), which are two
%         elements in the Wishard random matrix.
%
%INPUTS: ijkl A 4X1 of 1X4 vector of the indices.
%           A The dXd positive-definite, symmetric scale matrix of the
%             distribution.
%          nu The number of degrees of freedom of the distribution. Note
%             that nu>=d.
%
%OUTPUTS: val The desired covariance value.
%
%An expression for the covariance is given on page 99 of [1].
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    i=ijkl(1);
    j=ijkl(2);
    k=ijkl(3);
    l=ijkl(4);

    val=nu*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
end

function matVal=secondMatMoment(A,nu)
%%SECONDMATMOMENT Compute the second matrix moment of the Wishart
%                 distribution. That is E{X*X}.
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.  
%
%OUTPUTS: matVal The dXd second moment matrix.
%
%An expression for the second moment is given on page 98 of [1].
%
%EXAMPLE:
%Here, we compare the second-order moments obtained from Monte Carlo runs
%to those obtained from this function.
% A=[55,  7, 12, 17;
%     7, 45, 17, 22;
%    12, 17, 35, 27;
%    17, 22, 27, 25];
% nu=6;
% numMC=1e5;
% EX2=zeros(4,4);
% for k=1:numMC
%     X=WishartD.rand(A,nu);
%     EX2=EX2+X*X;
% end
% EX2=EX2/numMC;
% EX2Exact=WishartD.secondMatMoment(A,nu);
% RelErr=max(abs((EX2(:)-EX2Exact(:))./EX2Exact(:)))
%One will see that the Monte Carlo solution tends to have around 0.01
%maximum relative error and thus agrees with the exact solution.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    matVal=(nu+nu^2)*A*A+nu*trace(A)*A;
end


function matVal=secondMoments(A,nu)
%%SECONDMOMENTS Calculate the second moment (and all cross moments) of all
%     of the elements of the dXd Wishart distribution. Since the Wishart
%     distribution is a matrix distribution, the second moments are stored
%     in a tensor (4D matrix).
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%
%OUTPUTS: matVal The dXdXdXd matrix of second moments matVal(i1,i2,i3,i4)
%                is the expected value E(X(i1,i2)*X(i3,i4)).
%
%An expression for the second moment is given in Equation 3.3.15 on page 98
%of [1].
%
%EXAMPLE:
%Here, we compare the second-order moments obtained from Monte Carlo runs
%to those obtained from this function.
% A=[55,  7, 12, 17;
%     7, 45, 17, 22;
%    12, 17, 35, 27;
%    17, 22, 27, 25];
% nu=6;
% d=size(A,1);
% 
% numMonteCarlo=1e6;
% m2Mat=zeros(4,4,4,4);
% for curRun=1:numMonteCarlo
%     X=WishartD.rand(A,nu);
%     
%     %Find all second moments.
%     for i=1:d
%         for j=1:d 
%             for k=1:d
%                 for l=1:d
%                     val=m2Mat(i,j,k,l)+X(i,j)*X(k,l);
%                     m2Mat(i,j,k,l)=val;
%                 end
%             end
%         end
%     end
% end
% m2Mat=m2Mat/numMonteCarlo;
% m2Exact=WishartD.secondMoments(A,nu);
% %Compare to the analytic value.
% max(abs(m2Mat(:)-m2Exact(:))./abs(m2Exact(:)))
%One will see that the Monte Carlo solution tends to have around 0.3%
%maximum relative error and thus agrees with the exact solution.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(A,1);
    matVal=zeros(d,d,d,d);
    
    for i=1:d
        for j=i:d 
            for k=1:d
                for l=k:d
                    val=nu^2*A(i,j)*A(k,l)+nu*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
                    matVal(i,j,k,l)=val;
                    matVal(j,i,k,l)=val;
                    matVal(i,j,l,k)=val;
                    matVal(j,i,l,k)=val;
                end
            end
        end
    end
end

function val=secondMomentTerm(ijkl,A,nu)
%%SECONDMOMENTTERM Compute the expected value E(X(i,j)*X(k,l)), where the
%         variables are two elements in the Wishard random matrix.
%
%INPUTS: ijkl A 4X1 of 1X4 vector of the indices.
%           A The dXd positive-definite, symmetric scale matrix of the
%             distribution.
%          nu The number of degrees of freedom of the distribution. Note
%             that nu>=d.
%
%OUTPUTS: val The desired second moment value.
%
%An expression for the second moment is given in Equation 3.3.15 on page 98
%of [1].
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions. Boca Raton:
%    Chapman & Hall/CRC, 2000.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    i=ijkl(1);
    j=ijkl(2);
    k=ijkl(3);
    l=ijkl(4);

    val=nu*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
end

function val=PDF(X,A,nu)
%%PDF Evaluate the Wishart distribution for a particular matrix.
%
%INPUTS: X The dXd symmetric matrix at which the PDF of the distribution
%          is to be evaluated.
%        A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%
%OUTPUTS: val The value of the PDF of the Wishart distribution at X with
%             the given parameters.
%
%An expression for this Wishart distribution with the normalization
%constant is given in Chapter 2.3.6 of [1].
%
%REFERENCES:
%[1] C. M. Bishop, Pattern Recognition and Machine Learning. Cambridge,
%    United Kingdom: Springer, 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(X,1);
    
    if(rank(X)<d)
        val=0;
        return;
    end
    
    %The normalization constant
    B=1/(det(A)^(nu/2)*2^(nu*d/2)*gammaMultiDim(nu/2,d));
    
    val=B*det(X)^((nu-d-1)/2)*exp(-(1/2)*trace(A\X));
end
    
function X=rand(A,nu,algorithm)
%%RAND Generate a Wishart random matrix with the given number of degrees of
%      freedom and scale matrix.
%
%INPUTS: A The dXd positive-definite, symmetric scale matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d. This needn't be an integer if algorithm=1.
% algorithm An optional parameter specifying the algorithm to use to
%          generate the random variables. If nu is an integer and
%          nu<=d*(d-1)/2 (the number of values on the lower half of a dXd
%          matrix, including the diagonal), then the default if omitted or
%          an empty matrix is passed is 0; otherwise it is 1. Possible
%          values are:
%          0 Generate the random Wishart matrix from the outer product of
%            matrices of normal random variables as described in Chapters 5
%            and 8 [1] (and in numerous other textbooks). This can only be
%            done if nu is an integer.
%          1 Generate the random Wishart matrix using Bartlett's
%            decomposition as in [2].
%
%OUTPUTS: X A dXd Wishart random matrix.
%
%The Wishart distribution and its generation are discussed in Chapters 5
%and 8 of [1]. The Wishart distribution is just the outer product of dXnu
%normally generated random matrices, where each row is generated having
%zero-mean and covariance matrix A.
%
%However, using such a method is inefficient if nu is large. The technique
%based on Bartlett's decomposition in [2], which is algorithm 1, requires
%fewer random variables. That method generates random samples with the
%correct number of degrees of freedom, but having A=eye(d,d). Thus, a
%transformation is needed to obtain the correct distribution. Since a
%Wishart random variable with scale matrix Sigma arises from the non-
%normalized sum of the outer product of zero-mean normal random variables,
%if Xorig is the sample with A=eye(d,d), then to get a desired A, we use
%SA*X*SA', where SA is the lower-triangular Cholesky decomposition of A.
%
%A different derivation of Bartlett's decomposition is given in [3], where
%it becomes clear that the formulation  in [2], in terms of normal
%distributions and chi-squared distributions, is also valid if nu is not an
%integer.
%
%REFERENCES:
%[1] M. L. Eaton, Multivariate Statistics: A Vector Space Approach, ser.
%    Institute of Mathematical Statistics Lecture Notes-Monograph Series.
%    Beachwood, Ohio: Institute of Mathematical Statistics, 2007, vol. 53.
%[2] A. M. Kshirsagar, "Bartlett decomposition and Wishart distribution,"
%    The Annals of Statistics, vol. 30, no. 1, pp. 239-241, Mar. 1959.
%[3] D. G. Kabe, "A note on the Bartlett decomposition of a Wishart
%    matrix," Journal of the Royal Statistical Society Series B
%    (Methodological), vol. 26, no. 2, pp. 270-273, 1963.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

d=size(A,1);

if(nargin<3||isempty(algorithm))
    if(fix(nu)~=nu||nu>d*(d-1)/2)
        algorithm=1; 
    else
        algorithm=0;
    end
end

switch(algorithm)
    case 0
        S=chol(A,'lower')*randn(d,nu);
        X=S*S';
    case 1%Use the algorithm of [2].
        B=zeros(d,d);
        B(1,1)=sqrt(ChiSquareD.rand(1,nu));
        for i=2:d
            B(i,i)=sqrt(ChiSquareD.rand(1,nu-(i-1)));
            B(i,1:(i-1))=randn(1,i-1);
        end
        SA=chol(A,'lower');

        X=SA*(B*B')*SA';
    otherwise
        error('Unknown algorithm chosen.')
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
