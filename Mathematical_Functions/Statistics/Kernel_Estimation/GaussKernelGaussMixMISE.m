function MISE=GaussKernelGaussMixMISE(n,H,w,mu,P)
%%GAUSSKERNELGAUSSMIXMISE Given a Gaussian mixture probability density
%            function (PDF) that one is attempting to approximate using a
%            kernel density estimator with Gaussian kernels all having
%            square root covariance matrix H, determine the mean integrated
%            squared error (MISE) of the estimate. The MISE is the expected
%            value of int_{R}(\hat{f}(x)-f(x))^2dx where \hat{f}(x) is the
%            kernel estimator based on a set of random samples of the true
%            density f(x) (since it is based on random samples, the mean
%            part is the expected value being taken over all possible
%            samples). This function is useful for assessing different
%            kernel bandwidth estimators (estimators of H) when using
%            Gaussian mixture sample problems.
%
%INPUTS: n The number of samples that the kernel-density estimator will
%          use.
%        H The dXd bandwidth matrix used in the kernel density estimator
%          (as one would use in the function kernelApprox).
%        w A numCompX1 or 1XnumComp vector of the weights associated with
%          each of the numComp components of the Gaussian mixture.
%       mu The dXnumComp matrix of the mean of each component of the
%          Gaussian mixture.
%        P The dXdXnumComp set of covariance matrices of the components of
%          the Gaussian mixture.
%
%OUTPUTS: MISE The MISE of thekernel estimator with the given number of
%              samples and bandwidth matrix.
%
%This function implements the formula for the MISE in an arbitrary number
%of dimensions in Section 4 of [1].
%
%REFERENCES:
%[1] M. P. Wand and M. C. Jones, "Comparison of smoothing parameterizations
%    in bivariate kernel density estimation," Journal of the American
%    Statistical Association, vol. 88, no. 422, pp. 520-528, Jun. 1993.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of components in the 
k=size(mu,2);
d=size(mu,1);

%The definition of the bandwidth used in the paper is the square of what is
%used here.
H=H*H';

Omega0=zeros(k,k);
Omega1=zeros(k,k);
Omega2=zeros(k,k);

for l=1:k
    for lp=l:k
       Sigma=P(:,:,l)+P(:,:,lp);
       
       Omega0(l,lp)=GaussianD.PDF(mu(:,l),mu(:,lp),Sigma);
       Omega1(l,lp)=GaussianD.PDF(mu(:,l),mu(:,lp),Sigma+H);
       Omega2(l,lp)=GaussianD.PDF(mu(:,l),mu(:,lp),Sigma+2*H);
    end
end

MISE=(1/n)*(4*pi)^(-d/2)*1/sqrt(det(H))+w(:)'*((1-1/n)*Omega2-2*Omega1+Omega0)*w(:);

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
