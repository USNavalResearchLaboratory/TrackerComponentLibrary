function [xi,xPred,PPred]=purePropDiscPred(xi,w,f,Q,stateDiffTrans,stateTrans)
%%PUREPROPDISCPRED Perform the state prediction step of a discrete-time
%                  pure propagation (cubature) filter. Such a filter is
%                  similar to a cubature (unscented) Kalman filter in that
%                  cubature points are propagated, but differs in that the
%                  points are never resampled from the first two moments.
%                  The input and primary output of this function is a set
%                  of cubature points. Thus, information on higher order
%                  moments is, to an extent, preserved. The filter 
%                  propagates the points, which can then be provided to the
%                  measurement update function. The first two moments (mean
%                  and covariance) are also available for display, but do
%                  not directly play a role in the filter.
%
%INPUTS: xi An xDim X numCubPoints matrix of sample (cubature) points that
%           are to be propagated. These are akin to target state
%           hypotheses.
%         w A numCubPoints X 1 vector of the weights associated with the
%           sample points. These are not changed during the propagation of
%           the points.
%         f A function handle for the state transition function that takes
%           the state as its parameter.
%         Q The xDimX xDim process noise covariance matrix for an additive
%           Gaussian process noise model.
% stateDiffTrans An optional function handle that takes an xDimXN matrix
%                of N differences between states estimates and
%                transforms them however might be necessary. For
%                example, a state continaing angular components will
%                generally need differences between angular components
%                wrapped to the range +/-pi.
% stateTrans An optional function that takes a state estimate and
%            transforms it. This is useful if one wishes the elements of
%            the state to be bound to a certain domain. For example, if an
%            element of the state is an angle, one should generally want to
%            bind it to the region +/-pi. This is not applied to the output
%            of f.
%
%OUTPUTS: xi The xDim X numCubPoints matrix of propagated sample points. w
%            is still the associated set of weights.
%      xPred The weighed mean of the predicted sample points.
%      PPred The covariance matrix associated with the predicted sample
%            points. If the propagated sample points before adding in the
%            effects of the additive noise were not full rank, then this
%            might not equal the true predicted sample covariance matrix.
%
%A discrete dynamic model of
%x_{k+1}=f(x_{k})+w_k
%where x_k is the state at discrete time k, and x_{k+1} is the state at
%discrete time k+1. w_k is zero-mean additive Gaussian noise is assumed.
%
%The algorithm is based on that of [1], but in [1], the sigma-point
%covariance increase method is specifically tailored for the points used in
%a standard unscented Kalman filter. We would prefer to use a completely
%arbitrary set of cubature points. We could create a similar contiuous-time
%Riccatti equation to transform
%Pxx=sum_i w(i)*xiBar(i)*xiBar(i)'
%into
%Pxx+Q=sum_i w(i)*(xiBar(i)+U(i))*(xiBar(i)+U(i))'
%for an arbitrary set of points, where xiBar(i)=xi(i)-xiBar, where xiBar is
%the weighted average of the cubature points. However, how one chooses the
%points U(i) will warp the higher moments in different ways. However, one
%would expect that all points should be equally affected. Thus, we choose
%to scale the points about the mean instead. If Pxx is positive definite,
%then there should be a solution.
%
%Basically, we want to find a lower-triangular matrix C such that
%C*P*C'=P+Q
%Say P=SP*SP' and P+Q=SPQ*SPQ', where SP and SPQ can be found via Cholesky
%decompositions. Then, C*SP*SP'*C'=SPQ*SPQ'. Thus, if we say that C*SP=SPQ,
%then C=SPQ*inv(SP). We then replace xi(i) with C*(xi(i)-xiBar)+xiBar. The
%pinv function and the cholSemiDef function are used to help lessen
%problems with nearly singular matrices.
%
%REFERENCES:
%[1] O. Straka, J. Duník, and M. Simandl, "Design of pure propagation
%    unscented Kalman filter," in Proceedings of the 19th World Congress of
%    The International Federation of Automatic Control, Cape Town, South
%    Africa, 24-29 Aug. 2014, pp. 5399-5938.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(stateTrans))
   stateTrans=@(x)x; 
end

if(nargin<5||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

xDim=size(xi,1);
numCubPoints=length(w);

%Propagate the state samples
for curP=1:numCubPoints
    xi(:,curP)=f(xi(:,curP));
end

%Instead of the sigma-point covariance method to increase the covariance,
%we use the scaling described above, which works with an arbitrary sample
%set.
xPred=sum(bsxfun(@times,xi,w'),2);%The sample mean
xi=stateDiffTrans(bsxfun(@minus,xi,xPred));%Center the samples.

%Compute the covariance matrix
P=zeros(xDim,xDim);
for curP=1:numCubPoints
    P=P+w(curP)*(xi(:,curP)*xi(:,curP)');
end

SP=cholSemiDef(P,'lower');
SPQ=cholSemiDef(P+Q,'lower');
C=SPQ*pinv(SP);

%The updated points.
xi=stateTrans(bsxfun(@plus,C*xi,xPred));

if(nargout>2)
    %This assumes that P was positive definite.
    PPred=P+Q;
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
