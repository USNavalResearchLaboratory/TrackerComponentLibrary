function [H,mu,SR]=standardEpanechnikovKernelBW(param1,N)
%%STANDARDEPANECHNIKOVKERNELBW Estimate the (univariate or multivariate)
%           kernel bandwidth when approximating a set of samples whose
%           underlying density is (or is approximated as) Gaussian with an
%           Epanechnikov kernel. The Epanechnikov kernel has certain
%           desirable optimality properties.
%
%INPUTS: param1 If inputType=0, then this is:
%              xi A numDimXN set of samples of the distribution from which
%                 a Gaussian kernel estimator should interpolate a PDF. For
%                 multidimensional estimation, there must be >=numDim
%                 samples. Specifically, the sample covariance matrix must
%                 be positive definite.
%              If inputType=1, then this is:
%              SR A root covariance matrix such that SR*SR' equals the 
%                 (exact or approximated) covariance matrix of the PDF
%                 being approximated.
%            N If this parameter is provided, then param1 is assumed to be
%              SR and N is the number of samples that are being smoothed.
%              Otherwise, if this is omitted or an empty matrix is passed,
%              param1 is xi.
%
%OUTPUTS: H The numDimXnumDim kernel bandwidth.
%         mu If xi is given, then this is the sample mean. Otherwise, this
%            is an empty matrix.
%         SR The root covariance matrix from the input or computed from xi.
%
%The formulae for the optimal bandwidth of an Epanechnikov kernel given
%samples from an underlying Gaussian distribution is given in Equation 2.3
%of [1]. The function EpanechnikovKernel implements the Epanechnikov
%kernel.
%
%REFERENCES:
%[1] C. Musso, N. Oudjane, and F. LeGland, "Improving regularised particle
%    filters," in Sequential Monte Carlo Methods in Practice, A. Doucet, J.
%    F. G. de Freitas, and N. J. Gordon, Eds., Ch. 12, New York:
%    Springer-Verlag, 2001.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    N=[];
end

nx=size(param1,1);%Dimensionality.

if(isempty(N))%xi is given
    xi=param1;
    N=size(xi,2);
    [mu,S]=calcMixtureMoments(xi);
    SR=cholSemiDef(S,'lower',1);%Not fully triangular. but S=SR*SR'
else%SR is given and N must be given.
    mu=[];
    SR=param1;
end

cn=hypersphereVolume(1,nx);
A=((8/cn)*(nx+4)*2*sqrt(pi)^nx)^(1/(nx+4));
hOpt=A*N^(-1/(nx+4));

H=hOpt*SR*eye(nx);

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
