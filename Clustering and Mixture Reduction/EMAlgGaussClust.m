function [w,mu,P]=EMAlgGaussClust(z,w,mu,P,numIter)
%%EMALGGAUSSCLUST Use the expectation maximization (EM) algorithm to
%                 refine estimates of the components of a Gaussian mixture
%                 with a known number of terms given a set of samples.
%
%INPUTS: z A zDim X numPoints set of samples of the Gaussian mixture.
%        w A KX1 set of initial weight estimates of the K Gaussians in the
%          mixture.
%       mu A zDim X K set of initial mean estimates of the component
%          Gaussians in the mixture.
%        P A zDim X zDim X K hypermatrix of initial covariance matrix
%          estimates of the components of the Gaussians in the mixture.
%  numIter The number of iterations of the EM algorithm to perform.
%
%OUTPUTS: w The refined weights.
%        mu The refined means.
%         P The refined covariance matrix estimates.
%
%The EM algorithm for Gaussian mixtures is an implementation of the
%algorithm described in Chapter 9.2.2 of [1].
%
%REFERENCES:
%[1] C. M. Bishop, Pattern Recognition and Machine Learning. Cambridge,
%    United Kingdom: Springer, 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %zDim=size(z,1);
    numPoints=size(z,2);
    K=size(mu,2);

    gamma=zeros(numPoints,K);
    for curIter=1:numIter
        %Calculate the posterior weights.
        for k=1:K
            gamma(:,k)=GaussianPDF(z,mu(:,k),P(:,:,k))*w(k);
        end
        %Normalize the weights for each measurement.
        gamma=bsxfun(@rdivide,gamma,sum(gamma,2));

        %Update the means, covariances and weights using the posterior
        %weights.
        for k=1:K
            Nk=sum(gamma(:,k));
            w(k)=Nk/numPoints;
            
            [mu(:,k), P(:,:,k)]=calcMixtureMoments(z,gamma(:,k)/Nk);
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
