function [xUpdate,PUpdate,innov,Pzz,W]=EKFUpdateWithPred(z,R,zPred,PzPred,otherInfo,numIter,innovTrans,stateDiffTrans,stateTrans)
%%EKFUPDATEWITHPRED Given the output of the measurement prediction step
%           from EKFMeasPred and a measurement, complete the measurement
%           update step of the first or second order Extended Kalman Filter
%           (EKF) or iterated Extended Kalman Filter (IEKF) Kalman filter.
%           Separating the measurement prediction step from the rest of the
%           update step can make the creation of multiple measurement
%           association hypotheses from a single target prediction more
%           efficient. The full measurement update function is EKFUpdate.
%
%INPUTS: z The zDim X 1 vector measurement.
%        R The zDim X zDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimX1 measurement prediction from the filter.
%   PzPred The zDimXzDim covariance matrix associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the EKFMeasPred function.
%  numIter The number of iterations to perform if an iterated EKF is
%          desired. If this parameter is omitted or an empty matrix is
%          passed, the default value of zero is used. That is, just use
%          the standard update without any additional iterations.
% innovTrans An optional function handle that transforms the value of the
%          difference between the observation and any predicted points.
%          This only needs to be supplied when a measurement difference
%          must be restricted to a certain range. For example, the
%          innovation between two angles will be 2*pi if one angle is zero
%          and the other 2*pi, even though they are the same direction. In
%          such an instance, a function handle to the wrapRange function
%          with the appropriate parameters should be passed for innovTrans.
% stateDiffTrans An optional function handle that, like innovTrans does for
%          the measurements, a difference between states and transforms it
%          however might be necessary. For example, a state containing
%          angular components will generally need to be transformed so
%          that the difference between the angles is wrapped to -pi/pi.
% stateTrans An optional function that takes a state estimate and
%          transforms it. This is useful if one wishes the elements of the
%          state to be bound to a certain domain. For example, if an
%          element of the state is an angle, one should generally want to
%          bind it to the region +/-pi.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update.
%
%See the comments to the function EKFMeasPred for an example of usage of
%this function. See the comments to EKFUpdate for more information on
%the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zDim=size(z,1);

if(nargin<6||isempty(numIter))
    numIter=0;
end

if(nargin<7||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(x)x;
end

if(nargin<8||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

if(nargin<9||isempty(stateTrans))
    stateTrans=@(x)x;
end

HJacob=otherInfo.HJacob;
HHessian=otherInfo.HHesssian;
measPredTrans=otherInfo.measPredTrans;
Pxz=otherInfo.Pxz;
h=otherInfo.h;
H=otherInfo.H;
xPred=otherInfo.xPred;
PPred=otherInfo.PPred;

Pzz=PzPred+R;

innov=innovTrans(z-zPred);
W=Pxz/Pzz;

xUpdate=stateTrans(xPred+W*innov);

temp=W*H;
temp=eye(size(temp))-temp;
PUpdate=temp*PPred*temp'+W*R*W';
%Ensure symmetry
PUpdate=(PUpdate+PUpdate')/2;

for curIter=1:numIter
    if(isa(HJacob,'function_handle'))
        H=HJacob(xUpdate);
    else
        H=HJacob;%If the Jacobian matrix was directly given.
    end

    if(~isempty(HHessian))
        if(isa(HHessian,'function_handle'))
            HHess=HHessian(xUpdate);
        else
            HHess=HHessian;%If the hessian matrix was directly given.
        end

        [zPredHessTerm,PzzHessTerm]=getHessianTerms(HHess,PPred,zDim);
    else
        zPredHessTerm=0;
        PzzHessTerm=0;
    end
    
    %Update the PEstimate
    Pzz=R+H*PPred*H'+PzzHessTerm;
    W=PPred*H'/Pzz;
    temp=W*H;
    temp=eye(size(temp))-temp;
    PUpdate=temp*PPred*temp'+W*R*W';
    %Ensure symmetry
    PUpdate=(PUpdate+PUpdate')/2;
    
    zPred=measPredTrans(h(xUpdate)+zPredHessTerm);
    
    %Update the x estimate
    xUpdate=stateTrans(xUpdate+PUpdate*H'*lsqminnorm(R,innovTrans(z-zPred))-PUpdate*lsqminnorm(PPred,stateDiffTrans(xUpdate-xPred)));
end

end

function [zPredHessTerm,PzzHessTerm]=getHessianTerms(HHess,PPred,zDim)
%The function returns the terms associated with the second derivative in
%the second-order EKF.

PzzHessTerm=zeros(zDim,zDim);
zPredHessTerm=zeros(zDim,1);
for n=1:zDim
    en=zeros(zDim,1);
    en(n)=1;
    
    HPProdn=HHess(:,:,n)*PPred;
    zPredHessTerm=zPredHessTerm+en*trace(HPProdn);
    
    for m=1:zDim
        em=zeros(zDim,1);
        em(m)=1;
        PzzHessTerm=PzzHessTerm+en*em'*trace(HPProdn*HHess(:,:,m)'*PPred);
    end
end
PzzHessTerm=PzzHessTerm/2;
zPredHessTerm=zPredHessTerm/2;
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
