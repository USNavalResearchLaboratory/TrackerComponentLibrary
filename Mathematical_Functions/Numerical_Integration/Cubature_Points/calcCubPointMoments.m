function [mu,P]=calcCubPointMoments(z,S,h,xi,w,innovTrans,meanFun)
%%CALCCUBPOINTMOMENTS Using cubature integration with specified cubature
%                     points and weights, determine the first two moments
%                     of a transformation of a value x that is corrupted
%                     with additive Gaussian. x is a Gaussian distributed
%                     random variable with mean z and covariance matrix
%                     S*S'. The mean vector mu and covariance matrice P
%                     returned are those of h(x), where h is a function
%                     provided. Cubature integration using the provided
%                     cubature points and weights.
%
%INPUTS: z The numDimX1 value added to the zero-mean Gaussian noise with
%          covariance R that is driving the system.
%        S The numDimXnumDim lower-triangular square root covariance matrix
%          of the zero-mean Gaussian noise driving the system. If R is the
%          covariance matrix of the noise, then S=cholSemiDef(R,'lower');
%        h A function handle for how the vector z is transformed. This
%          function must be able to take a numDimXnumPoints matrix of
%          transformed z values and return a numDimOutXnumPoints matrix of
%          the transformed vectors. The return values of the
%          calcCubPointMoments function are the mean and covariance matrix
%          of h(x), where x is Gaussian noise with covariance matrix S*S'
%          and mean z.
%       xi An optional numDimXnumPoints matrix of cubature points for
%          normal 0-I distribution. If this and w are omitted, then a
%          default (probably inefficient) set of cubature points and
%          weights is generated using fifthOrderCubPoints.
%        w A numPointsX1 vector of the weights associated with each of
%          the cubature points in xi. Note that all w>=0 and normally
%          sum(w)=1. If this and xi are omitted, then a default (probably
%          inefficient) set of cubature points and weights is generated
%          using fifthOrderCubPoints.
% innovTrans An optional function handle that computes differences
%          between converted values. This can be useful for example, when
%           taking the difference of circular data, since the difference
%          between, for example, pi and -pi should be 0, not 2*pi. This is
%          called as innovTrans(a,b) and the default if omitted or an empty
%          matrix is passed is @(a,b)bsxfun(@minus,a,b); This function
%          must be able to handle multiple values of a at once.
%  meanFun An optional function handle that, when given N values with
%          weights, produces the weighted average. this is called as
%          meanFun(x,w), where w are the weights. This function only has to
%          be provided if the domain of the converted values is not linear.
%          For example, when averaging angular values, then the function
%          meanAng should be used.
%
%OUTPUTS: mu The mean of h(x) found using cubature integration.
%          P The covariance matrix of h(x) found using cubature
%            integration.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(z,1);

if(nargin<4||isempty(xi))
    [xi,w]=fifthOrderCubPoints(numDim);
end

if(nargin<6||isempty(innovTrans))
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

xi=h(transformCubPoints(xi,z,S));
if(nargin<7||isempty(meanFun))
    mu=sum(bsxfun(@times,xi,w(:).'),2);%The mean
else
    mu=meanFun(xi,w(:).');
end

diffVal=innovTrans(xi,mu);

xDim=size(xi,1);
numCubPts=size(xi,2);
P=zeros(xDim,xDim);
for curPoint=1:numCubPts
    P=P+w(curPoint)*(diffVal(:,curPoint)*diffVal(:,curPoint)');
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
