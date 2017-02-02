function [xSmooth,SSmooth]=sqrtCubKalSmoother(xInit,SInit,z,h,f,SR,SQ,xi,w,kD,innovTrans,measAvgFun,stateDiffTrans,stateTrans,stateAvgFun)
%%SQRTCUBKALSMOOTHER  Run the forward-backward square root cubature Kalman
%                 smoother on a batch of measurements. The smoothed result
%                 at one time step or along the entire batch are available.
%                 The initial predicted states can not be uninformative.
%INPUTS: xInit  The predicted state at the time of the initial measurement
%               in z.
%        SInit  The square root covariance matrix associated with the
%               predicted state at the time of the initial measurement in
%               z.
%           z   The zDim X N matrix of measurements for the whole batch.
%           h   A NX1 cell array of function handles for the measurement
%               function that transform the state into the measurement
%               domain at each step. If the same measurement function is
%               used for all steps in the batch, then h can just be the
%               single function handle used.
%           f   A NX1 cell array of function handles for the state
%               transition function that transform the state into the
%               measurement domain at each step. If the same measurement
%               function is used for all steps in the batch, then h can
%               just be the single function handle used.
%          SR   The zDim X zDim X N hypermatrix of square root measurement
%               covariance matrices. Alternatively, if all of the
%               measurement covariance matrices are the same, one can just
%               pass a single zDim X zDim matrix.
%          SQ   The xDim X xDim X (N-1) hypermatrix of square root process
%               noise covariance matrices. Alternatively, if all of the
%               process noise covariance matrices are the same, one can
%               just pass a single xDim X xDim matrix.
%          xi   An xDim X numCubPoints matrix of cubature points.        
%           w   A numCubPoints X 1 vector of the weights associated with
%               the cubature points.
%          kD   The discrete time-step at which the smoothed state estimate
%               is desired, where z(:,1) is at discrete time-step 1 (not
%               0). If kD is omitted or an empty matrix is passed, then
%               results along the entire batch are obtained.
%  innovTrans   An optional function handle that transforms the value of
%               the difference between the observation and any predicted
%               points. This must be able to handle sets of differences.
%               For a zDim measurement, this must be able to handle a
%               zDimXN matrix of N differences. This only needs to be
%               supplied when a measurement difference must to be
%               restricted to a certain range. For example, the innovation
%               between two angles will be 2*pi if one angle is zero and
%               the other 2*pi, even though they are the same direction. In
%               such an instance, a function handle to the wrapRange
%               function with the appropriate parameters should be passed
%               for innovTrans.
%  measAvgFun   An optional function handle that, when given N measurement
%               values with weights, produces the weighted average. This
%               function only has to be provided if the domain of the
%               measurement is not linear. For example, when averaging
%               angular values, then the function meanAng should be used.
%stateDiffTrans  An optional function handle that, like innovTrans does for
%               the measurements, takes an xDimXN matrix of N differences
%               between states and transforms them however might be
%               necessary. For example, a state continaing angular
%               components will generally need to be transformed so that
%               the difference between the angles is wrapped to -pi/pi.
%  stateTrans   An optional function that takes a state estimate and
%               transforms it. This is useful if one wishes the elements of
%               the state to be bound to a certain domain. For example, if
%               an element of the state is an angle, one should generally
%               want to bind it to the region +/-pi.
% stateAvgFun   An optional function that given an xDimXN matrix of N state
%               estimates and an NX1 vector of weights, provides the
%               weighted average of the state estimates. This is necessary
%               if, for example, states with angular components are
%               averaged.
%
%OUTPUTS:xSmooth The xDimXN smoothed state estimates at all steps if kD is
%                not provided or the xDimX1 smoothed information state
%                estimate at step kD if kD is provided.
%        SSmooth The lower-triangular square root covariance matrices
%                associated with the smoothed state estimates. This is
%                xDimXxDimXN for the whole batch if kD is not provided and
%                is xDimXxDim if kD is provided.
%
%The mathematics behind the state prediction and update functions are
%described in more detail in Section IX of "Basic Tracking Using Nonlinear
%3D Monostatic and Bistatic Measurements" by David F. Crouse.
%
%The smoothing algorithm is taken from [1], but can use any type of
%cubature points, not just the third-order points used in [1].
%
%REFERENCES:
%[1] I. Arasaratnam and S. Haykin, Cubature Kalman smoothers, Automatica,
%    vol. 47, no. 10, pp. 2245-2250, Oct. 2011.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10)
    kD=[];
end
if(nargin<11 || isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(x)x;
end
if(nargin<12 || isempty(measAvgFun))
    measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
end

if(nargin<13 || isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

if(nargin<14 || isempty(stateTrans))
    stateTrans=@(x)x;
end
if(nargin<15 || isempty(stateAvgFun))
    stateAvgFun=@(x,w)calcMixtureMoments(x,w);
end

xDim=size(xInit,1);
nCub=size(xi,2);
N=size(z,2);

if(isa(h,'function_handle'))
    h=repmat({h},N,1);
end
if(isa(f,'function_handle'))
    f=repmat({f},N,1);
end
if(size(SR,3)==1)
    SR=repmat(SR,1,1,N);
end
if(size(SQ,3)==1)
    SQ=repmat(SQ,1,1,N);
end

%Run the sqrt cubature Kalman filter forward
xPred=zeros(xDim,N);
SPred=zeros(xDim,xDim,N);
xPropCenPoints=zeros(xDim,nCub,N);
xUpd=zeros(xDim,N);
SUpd=zeros(xDim,xDim,N);

%The first step uses the priors
xPred(:,1)=xInit;
SPred(:,:,1)=SInit;

%The rest of the steps
for curStep=1:(N-1)
    [xUpd(:,curStep),SUpd(:,:,curStep)]=sqrtCubKalUpdate(xPred(:,curStep),SPred(:,:,curStep),z(:,curStep),SR(:,:,curStep),h{curStep},xi,w,innovTrans,measAvgFun,stateDiffTrans,stateTrans);
    [xPred(:,curStep+1),SPred(:,:,curStep+1),xPropCenPoints(:,:,curStep+1)]=sqrtDiscCubKalPred(xUpd(:,curStep),SUpd(:,:,curStep),f{curStep},SQ(:,:,curStep),xi,w,stateDiffTrans,stateAvgFun,stateTrans);
end
[xUpd(:,N),SUpd(:,:,N)]=sqrtCubKalUpdate(xPred(:,N),SPred(:,:,N),z(:,N),SR(:,:,N),h{N},xi,w,innovTrans,measAvgFun,stateDiffTrans,stateTrans);

%Run the backwards smoother
xSmooth=zeros(xDim,N);
SSmooth=zeros(xDim,xDim,N);

xSmooth(:,N)=xUpd(:,N);
SSmooth(:,:,N)=SUpd(:,:,N);
for curStep=(N-1):-1:1
    %Forward cubature state points
    xFwdPoints=stateTrans(transformCubPoints(xi,xUpd(:,curStep),SUpd(:,:,curStep)));
    
    %Forward, centered cubature state points
    xFwdCenPoints=bsxfun(@times,stateDiffTrans(bsxfun(@minus,xFwdPoints,xUpd(:,curStep))),sqrt(w)');
    
    %Root adjustment covariance
    Spp=tria([xPropCenPoints(:,:,curStep+1),SQ(:,:,curStep+1)]);
    
    %The cross covariance.
    Pfp=xFwdCenPoints*xPropCenPoints(:,:,curStep+1)';
    
    %The smoother gain
    G=(Pfp/Spp')/Spp;
    
    %Smooth state estimate
    xDiff=stateDiffTrans(xSmooth(:,curStep+1)-xPred(:,curStep+1));
    xSmooth(:,curStep)=stateTrans(xUpd(:,curStep)+G*xDiff);
    
    %Smooth state root covariance
    SSmooth(:,:,curStep)=tria([stateDiffTrans(xFwdCenPoints-G*xPropCenPoints(:,:,curStep+1)),G*SQ(:,:,curStep+1),G*SSmooth(:,:,curStep+1)]);
end

if ~isempty(kD)
    xSmooth=xSmooth(:,kD);
    SSmooth=SSmooth(:,:,kD);
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
