function [xUpdate,PUpdate,xiPart,wPart]=GaussPartFilterUpdate(xPred,PPred,zLike,numParticles,importPDF,importSampFun)
%%GAUSSPARTFILTERUPDATE Perform the measurement update step in the Gaussian
%              particle filter. In this filter, the state is approximated
%              as a Gaussian variable, but the measurement has an arbitrary
%              continuous conditional density function that is provided.
%              The Gaussian particle filter is a Monte Carlo method and
%              the measurement update step tends to be worse than that of
%              the progressive Gaussian filter with a comparable number of
%              points. On the other hand, the Gaussian particle filter with
%              a very large number of points can be viewed as the best that
%              one can get with a Gaussian filter.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim predicted state covariance matrix.
%        zLike A function handle such that zLike(x) returns the conditional
%              likelikoods of the measurement given a matrix of of N target
%              states. That is, an xDimXN matrix. The likelihoods can be
%              returned as a 1XN vector or an NX1 vector. It does not
%              matter if they are all off by a constant multiplicative
%              factor from their true values.
% numParticles The (optional) number of particles to use in the Gaussian
%              particle filter. If this parameter is omitted or an empty
%              matrix is passed, the default of 1e4 particles is used.
%    importPDF The (optional) function handle for the PDF used for
%              importance sampling.  Specifically, importPDF(x,xPred,PPred)
%              is the importance  PDF evaluated at N target states in x (x
%              is an xDimXN  matrix). The parameters xPred,PPred are passed
%              in case the importance sampling function uses them. The
%              importance sampling PDF is ideally as close to the true
%              posterior PDF as possible. However, in practice, one will
%              typically just  use the prior (Gaussian) PDF for  importance
%              sampling, which is what is done if this parameter is omitted
%              or an empty matrix is passed.
% importSampFun If importPDF is provided (and is not an empty matrix), then
%              importSampFun is required. This is a function handle such
%              that importSampFun(numParticles,xPred,PPred) draws
%              numParticles samples of the importance sampling PDF.             
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector. If zLike returned all
%                 when evaluating the particles then an empty matrix will
%                 be returned to indicate a failure of the algorithm.
%         PUpdate The updated xDim X xDim state covariance matrix, unless
%                 an algorithmic failure occurs, in which case an empty
%                 matrix is returned.
%    xiPart,wPart The particles and their associated weights that were
%                 used. These are returned in case one wishes to use the
%                 alternate time update algorithm in the function
%                 GaussPartFilterPred.
%
%The algorithm is implemented as given in Table I of [1]. The algorithm
%take a Gaussian approxiamtion as an input and produces a Gaussian
%approximation as an output, using particles to directly evaluate Bayes'
%theorem.
%
%REFERENCES:
%[1] J. H. Kotecha and P. M. Djuric, "Gaussian particle filtering,"
%    IEEE Transactions on Signal Processing, vol. 51, no. 10, pp.
%    2592-2601, Oct. 2003.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(numParticles))
   numParticles=1e4; 
end

if(nargin<5||isempty(importPDF))
    xDim=size(xPred,1);
    
    %Use the predicted PDF as the proposal density if no proposal density
    %is given.
    importPDF=@(x,xPred,PPred)GaussianD.PDF(x,xPred,PPred);
    importSampFun=@(numPart,xPred,PPred)bsxfun(@plus,xPred,cholSemiDef(PPred,'lower')*randn(xDim,numPart));
end

%1) Draw samples from the importance function.
xiPart=importSampFun(numParticles,xPred,PPred);

%2) Obtain the weights
wPart=zLike(xiPart).*(GaussianD.PDF(xiPart,xPred,PPred)./importPDF(xiPart,xPred,PPred))';

%3) Normalize the weights
wPart=wPart/sum(wPart);

%If zLike returned all zeros, then the update failed and we will return
%empty matrices to indicate that.
if(any(~isfinite(wPart)))
    xUpdate=[];
    PUpdate=[];
    return;
end

%4) Estimate the mean and covariance
[xUpdate, PUpdate]=calcMixtureMoments(xiPart,wPart);

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
