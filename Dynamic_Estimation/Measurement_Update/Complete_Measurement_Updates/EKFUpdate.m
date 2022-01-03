function [xUpdate,PUpdate,innov,Pzz,W]=EKFUpdate(xPred,PPred,z,R,h,HJacob,numIter,HHessian,innovTrans,measPredTrans,stateDiffTrans,stateTrans)
%%EKFUPDATE Perform the measurement update step in the first- or second-
%           order Extended Kalman Filter (EKF) or iterated Extended Kalman
%           filter (IEKF), which is what one gets with numIter>0.
%
%INPUTS: xPred The xDimX1 predicted target state.
%        PPred The xDimXxDim predicted state covariance matrix.
%            z The zDimX1 vector measurement.
%            R The zDimXzDim measurement covariance matrix.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%       HJacob A function handle for the measurement Jacobian matrix that
%              takes the target state as a parameter or the zDimXxDim
%              measurement Jacobian matrix itself. Each column is the
%              derivative vector of h with respect to the corresponding
%              element of x. If not supplied or an empty matrix is passed,
%              then HJacob will be found using numerical differentiation
%              via the numDiff function with default parameters.
%      numIter The number of iterations to perform if an iterated EKF is
%              desired. If this parameter is omitted or an empty matrix is
%              passed, the default value of zero is used. That is, just
%              use the standard update without any additional iterations.
%     HHessian This parameter is only provided if a second-order EKF is
%              desired. This is either a function handle for the
%              measurement Hessian hypermatrix, or it is the measurement
%              Hessian hypermatrix itself. The matrix is
%              xDimXxDimXzDim. The matrix HHess=HHessian(x) is such that
%              HHess(i,j,k) is the second derivative of the kth element
%              of the vector returned by h with respect to the ith and jth
%              components of x. For each k, the Hessian matrix is
%              symmetric. If this parameter is omitted or an empty matrix
%              is passed, a first-order EKF update is used.
%   innovTrans An optional function handle that computes and optionally
%              transforms the value of the difference between the
%              observation and any predicted points. This is called as
%              innovTrans(a,b) and the default if omitted or an empty
%              matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%              able to handle sets of values. For a zDimX1 measurement,
%              either of the inputs could be zDimXN in size while one of
%              the inputs could be zDimX1 in size.  This only needs to be
%              supplied when a measurement difference must be restricted
%              to a certain range. For example, the innovation between two
%              angles will be 2*pi if one angle is zero and the other
%              2*pi, even though they are the same direction. In such an
%              instance, a function handle to the
%              wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%              appropriate parameters should be passed for innovTrans.
% measPredTrans An optional function handle that transforms the predicted
%              measurement into a particular range. The second order EKF
%              has a linear correction to the prediction of the
%              measurement, which might not be appropriate for all types of
%              measurements.
% stateDiffTrans An optional function handle that, like innovTrans does for
%              the measurements, a difference between states and transforms
%              it however might be necessary. For example, a state
%              containing angular components will generally need to be
%              transformed so that the difference between the angles is
%              wrapped to -pi/pi.
%   stateTrans An optional function handle that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one should generally
%              want to bind it to the region +/-pi.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of
%[1]. The Joseph-form covariance update given in Chapter 5 of the same
%book is used for improved numerical stability. The expressions for the
%second-order EKF come from Chapter 10.3.2 of [1].
%
%The iteration of the measurement equation for the iterated EKF is given in
%Chapter 10.5.2 of the same book. In [2], it is noted that the iteration
%need only be done over each measurement update, not over an entire batch
%of measurements.
%
%The optional parameter innovTrans is not described in the above reference,
%but allow for possible modifications to the filter as described in [3].
%The parameters have been added to allow the filter to be used with
%angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] T. H. Kerr, "Streamlining Measurement Iteration for EKF Target
%    Tracking," IEEE Transactions on Aerospace and Electronic Systems,
%    vol. 27, no. 2, pp. 408-421, Mar. 1991.
%[3] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zDim=size(z,1);
if(nargin<6||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

if(nargin<7||isempty(numIter))
    numIter=0;
end

if(nargin<8)
    HHessian=[]; 
end

if(nargin<9||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

if(nargin<10||isempty(measPredTrans))
    measPredTrans=@(x)x;
end

if(nargin<11||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

if(nargin<12||isempty(stateTrans))
    stateTrans=@(x)x;
end

if(isa(HJacob,'function_handle'))
    H=HJacob(xPred);
else
    H=HJacob;%If the Jacobian matrix was directly given.
end

if(~isempty(HHessian))
    if(isa(HHessian,'function_handle'))
        HHess=HHessian(xPred);
    else
        HHess=HHessian;%If the Hessian matrix was directly given.
    end
    
    [zPredHessTerm,PzzHessTerm]=getHessianTerms(HHess,PPred,zDim);
else
    zPredHessTerm=0;
    PzzHessTerm=0;
end

zPred=measPredTrans(h(xPred)+zPredHessTerm);
innov=innovTrans(z,zPred);

Pzz=R+H*PPred*H'+PzzHessTerm;
%Ensure symmetry
Pzz=(Pzz+Pzz')/2;
W=PPred*H'/Pzz;

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
    %Ensure symmetry
    Pzz=(Pzz+Pzz')/2;
    W=PPred*H'/Pzz;
    temp=W*H;
    temp=eye(size(temp))-temp;
    PUpdate=temp*PPred*temp'+W*R*W';
    %Ensure symmetry
    PUpdate=(PUpdate+PUpdate')/2;
    zPred=measPredTrans(h(xUpdate)+zPredHessTerm);
    
    %Update the x estimate
    xUpdate=stateTrans(xUpdate+PUpdate*H'*lsqminnorm(R,innovTrans(z,zPred))-PUpdate*lsqminnorm(PPred,stateDiffTrans(xUpdate-xPred)));
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
