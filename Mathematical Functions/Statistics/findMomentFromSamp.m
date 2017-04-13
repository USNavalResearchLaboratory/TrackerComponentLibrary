function theMoment=findMomentFromSamp(alpha,xSamp,wSamp,normalize)
%%FINDMOMENTFROMSAMP Given a set of samples of a univariate or multivariate
%               random variable, determine a particular (non-central)
%               sample moment. To find all sample moments up to a
%               particular degree, use the function findMomentsFromSamp.
%
%INPUTS: alpha The nX1 or 1Xn vector of exponents for the monomial term
%              for which the moment is desired. The moment will be the
%              weighted average of prod(xSamp.^alpha)
%        xSamp An nXnumSamp matrix of numSamp samples from which moments
%              are to be computed.
%        wSamp If the samples are weighted, then wSamp is a 1XnumSamp or
%              numSampX1 vector of the weights. If the samples are not
%              weighted, then an empty matrix can be passed.
%    normalize An optional boolean value indicating whether the weights in
%              wSamp should be normalized prior to use. The default if
%              this parameter is omitted is false. The weights cannot be
%              normalized if they sum to zero.
%
%OUTPUTS: theMoment The scalar value of the computed moment.
%
%Noncentral moments are unbiased if one just sums the values of the powers
%of the components of each sample and divides by the number of samples. For
%weighted samples, we multiply each product in the sum by the sample
%weight, and then divide by the sum of the weights.
%
%Note, however, that central samples computed in such a manner ARE biased.
%Hence the reason introductory textbooks use a 1/(numSamp-1) in front of
%variance computations. The bias in central moments (when just dividing by
%numSamp) comes from using the sample mean in place of the true mean to
%centralize the moment.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(xSamp,2);

if(nargin<4||isempty(normalize))
    normalize=false; 
end

if(nargin<3||isempty(wSamp))
   wSamp=ones(N,1)/N;
end

if(normalize)
    normConst=sum(wSamp);
else
    normConst=1;
end

theMoment=sum(bsxfun(@times,wSamp(:)',prod(bsxfun(@power,xSamp,alpha(:)),1)))/normConst;

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
