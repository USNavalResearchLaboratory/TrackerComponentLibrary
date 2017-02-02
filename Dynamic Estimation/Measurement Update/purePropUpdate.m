function [xi,xUpdate,PUpdate]=purePropUpdate(xi,w,z,R,h,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%%PUREPROPUPDATE Perform the measurement update step of a pure propagation
%                (cubature) filter. Such a filter is similar to a cubature
%                (unscented) Kalman filter in that cubature points are
%                updated, but differs in that differs in that the
%                points are never resampled from the first two moments.
%                The input and primary output of this function is a set
%                of cubature points. Thus, information on higher order
%                moments is, to an extent, preserved. The filter updates
%                the points, which can then be given to the propagation
%                function. The first two moments (mean and covariance) are
%                also available for display, but do not directly play a
%                role in the filter.
%
%INPUTS: xi An xDim X numCubPoints matrix of sample (cubature) points that
%           are to be propagated. These are akin to target state
%           hypotheses.
%         w A numCubPoints X 1 vector of the weights associated with the
%           sample points. These are not changed during the propagation of
%           the points.
%         z The zDim X 1 vector measurement.
%         R The zDim X zDim measurement covariance matrix.
%         h A function handle for the measurement function that takes the
%           state as its argument.
%innovTrans An optional function handle that transforms the value of the
%           difference between the observation and any predicted points.
%           This must be able to handle sets of differences. For a zDim
%           measurement, this must be able to handle a zDimXN matrix of N
%           differences. This only needs to be supplied when a measurement
%           difference must to be restricted to a certain range. For
%           example, the innovation between two angles will be 2*pi if one
%           angle is zero and the other 2*pi, even though they are the same
%           direction. In such an instance, a function handle to the
%           wrapRange function with the appropriate parameters should be
%           passed for innovTrans.
% measAvgFun An optional function handle that, when given N measurement
%           values with weights, produces the weighted average. This
%           function only has to be provided if the domain of the
%           measurement is not linear. For example, when averaging angular
%           values, then the function meanAng should be used.
%stateDiffTrans An optional function handle that, like innovTrans does for
%           the measurements, takes an xDimXN matrix of N differences
%           between states and transforms them however might be necessary.
%           For example, a state containing angular components will
%           generally need to be transformed so that the difference
%           between the angles is wrapped to -pi/pi.
%stateTrans An optional function that takes a matrix of N state estimates
%           and transforms them. This is useful if one wishes the elements
%           of the state to be bound to a certain domain. For example, if
%           an element of the state is an angle, one might generally want
%           to bind it to the region +/-pi.
%
%OUTPUTS: xi The xDim X numCubPoints matrix of updated sample points. w
%            is still the associated set of weights.
%    xUpdate The weighed mean of the updated sample points.
%    PUpdate The covariance matrix associated with the updated sample
%            points.
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
%The optional parameters innovTrans and measAvgFun are not described in the
%above reference, but allow for possible modifications to the filter as
%described in [2]. The parameters have been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(innov)[innov(1,:);
%                   wrapRange(innov(2,:),-pi,pi)];
%measAvgFun=@(z,w)[calcMixtureMoments(z(1,:),w);
%                  meanAng(z(2,:),w')];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] O. Straka, J. Duník, and M. Simandl, "Design of pure propagation
%    unscented Kalman filter," in Proceedings of the 19th World Congress of
%    the The International Federation of Automatic Control, Cape Town,
%    South Africa, 24-29 Aug. 2014, pp. 5399-5938.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(x)x;
end

if(nargin<7||isempty(measAvgFun))
   measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
end

if(nargin<8||isempty(stateDiffTrans))
   stateDiffTrans=@(x)x; 
end

if(nargin<9||isempty(stateTrans))
   stateTrans=@(x)x; 
end

zDim=size(z,1);
xDim=size(xi,1);
numCubPoints=size(xi,2);

%Predicted cubature measurement points
zPredPoints=zeros(zDim,numCubPoints);
for curP=1:numCubPoints
    zPredPoints(:,curP)=h(xi(:,curP));
end

%Measurement prediction.
zPred=measAvgFun(zPredPoints,w);

%Centered, predicted cubature measurement points, transformed as
%necessary to keep the values within a desired range.
zPredCenPoints=innovTrans(bsxfun(@minus,zPredPoints,zPred));

xPred=sum(bsxfun(@times,xi,w'),2);%The sample mean
xPredCenPoints=stateDiffTrans(bsxfun(@minus,xi,xPred));%Center the samples.

%The innovation covariance.
Pzz=R;
%The cross covariance.
Pxz=zeros(xDim,zDim);
for curP=1:numCubPoints
    Pzz=Pzz+w(curP)*(zPredCenPoints(:,curP)*zPredCenPoints(:,curP)');
    Pxz=Pxz+w(curP)*(xPredCenPoints(:,curP)*zPredCenPoints(:,curP)');
end

innovPoints=innovTrans(bsxfun(@minus,z,zPredPoints));

%The filter gain.
K=Pxz*pinv(Pzz);

%Propagate the sigma points through the update relation.
xi=stateTrans(xi+K*innovPoints);

%Instead of the sigma-point covariance method to increase the covariance,
%we use the scaling described above, which works with an arbitrary sample
%set. Note that the paper hasa type with respect to Q. This is the correct
%formulation of Q.
Q=K*R*K';

xUpdate=sum(bsxfun(@times,xi,w'),2);%The sample mean
%Center the propagated samples.
xi=stateDiffTrans(bsxfun(@minus,xi,xUpdate));

%Compute the covariance matrix
P=zeros(xDim,xDim);
for curP=1:numCubPoints
    P=P+w(curP)*(xi(:,curP)*xi(:,curP)');
end

SP=cholSemiDef(P,'lower');
SPQ=cholSemiDef(P+Q,'lower');
C=SPQ*pinv(SP);

%The updated points.
xi=stateTrans(bsxfun(@plus,C*xi,xUpdate));

if(nargout>1)
    xUpdate=sum(bsxfun(@times,xi,w'),2);%The sample mean
    %This assumes that stateTrans did not change the covariance.
    PUpdate=P+Q;
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
