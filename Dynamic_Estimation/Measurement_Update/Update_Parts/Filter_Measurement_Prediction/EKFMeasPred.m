function [zPred,PzPred,otherInfo]=EKFMeasPred(xPred,PPred,zDim,h,HJacob,HHessian,measPredTrans)
%%EKFMEASPRED Perform the measurement prediction part of the measurement
%           update step of the first or second order Extended Kalman Filter
%           (EKF) or iterated Extended Kalman Filter (IEKF) Kalman filter.
%           The function EKFUpdateWithPred can be used to complete the
%           measurement update. Separating the measurement prediction step
%           from the rest of the update step can make the creation of
%           multiple measurement association hypotheses from a single
%           target prediction more efficient. The full measurement update
%           function is EKFUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states.
%        PPred The xDimXxDimXnumComp predicted state covariance matrices
%         zDim The dimensionality of the measurement.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%       HJacob A function handle for the zDimXxDim measurement Jacobian
%              matrix that takes the target state as a parameter or the
%              zDimXxDimXnumComp collection of all numComp measurement
%              Jacobian matrices for all components itself. Each column is
%              the derivative vector of h with respect to the corresponding
%              element of x. If not supplied or an empty matrix is passed,
%              then HJacob will be found using numerical differentiation
%              via the numDiff function with default parameters.
%     HHessian This parameter is only provided if a second-order EKF is
%              desired. This is either a function handle for the
%              xDimXxDimXzDim measurement Hessian hypermatrix, or it is the
%              xDimXxDimXzDimXnumComp collection of all measurement Hessian
%              hypermatrices itself. The matrix HHess=HHessian(x) is such
%              that HHess(i,j,k) is the second derivative of the kth
%              element of the vector returned by h with respect to the ith
%              and jth components of x. For each k, the Hessian matrix is
%              symmetric. If this parameter is omitted or an empty matrix
%              is passed, a first-order EKF update is used.
% measPredTrans An optional function handle that transforms the predicted
%              measurement into a particular range. The second order EKF
%              has a linear correction to the prediction of the
%              measurement, which might not be appropriate for all types of
%              measurements.
%
%OUTPUTS: zPred The zDimXnumComp measurement predictions from the filter.
%        PzPred The zDimXzDimXnumComp covariance matrices associated with
%               zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to EKFUpdateWithPred
%               when updating with a measurement.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of
%[1]. The Joseph-form covariance update given in Chapter 5 of the same
%book is used for improved numerical stability. The expressions for the
%second-order EKF come from Chapter 10.3.2 of [1]. The iteration of the
%measurement equation for the iterated EKF is given in Chapter 10.5.2 of
%the same book. See the comments to EKFUpdate for more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%EKFUpdate in one step as with using EKFMeasPred followed by
%EKFUpdateWithPred.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% HJacob=@(x)[2*x(1), 2*x(2), 2*x(3), 2*x(4);
%                  1,      0,      0, -(3*sqrt(x(4)))/2];
% PPred=[28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13];
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% R=eye(zDim,zDim);%Measurement covariance matrix.
% numIter=1;
% %The update in one step.
% [xUpdate,PUpdate,innov,Pzz,W]=EKFUpdate(xPred,PPred,z,R,h,HJacob,numIter);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=EKFMeasPred(xPred,PPred,zDim,h,HJacob);
% [xUpdate1,PUpdate1,innov1,Pzz1,W1]=EKFUpdateWithPred(z,R,zPred,PzPred,otherInfo,numIter);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;PUpdate1(:)-PUpdate(:);innov1(:)-innov;Pzz1(:)-Pzz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);
numComp=size(xPred,2);

if(nargin<5||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

if(nargin<6||isempty(HHessian))
    HHessian=[]; 
end

if(nargin<7||isempty(measPredTrans))
    measPredTrans=@(x)x;
end

if(isa(HJacob,'function_handle'))
    H=zeros(zDim,xDim,numComp);
    for k=1:numComp
        H(:,:,k)=HJacob(xPred(:,k));
    end
else
    H=HJacob;%If the Jacobian matrix was directly given.
end

if(~isempty(HHessian))
    if(isa(HHessian,'function_handle'))
        HHess=zeros(xDim,xDim,zDim,numComp);
        for k=1:numComp
            HHess(:,:,:,k)=HHessian(xPred(:,k));
        end
    else
        HHess=HHessian;%If the Hessian matrix was directly given.
    end
    
    [zPredHessTerm,PzzHessTerm]=getHessianTerms(HHess,PPred);
end

zPred=zeros(zDim,numComp);
Pxz=zeros(xDim,zDim,numComp);
PzPred=zeros(zDim,zDim,numComp);
if(isempty(HHessian))
    for k=1:numComp
        zPred(:,k)=measPredTrans(h(xPred(:,k)));
        Pxz(:,:,k)=PPred(:,:,k)*H(:,:,k)';
        PzPred(:,:,k)=H(:,:,k)*PPred(:,:,k)*H(:,:,k)';
        %Ensure symmetry
        PzPred(:,:,k)=(PzPred(:,:,k)+PzPred(:,:,k)')/2;
    end
else
    for k=1:numComp
        zPred(:,k)=measPredTrans(h(xPred(:,k))+zPredHessTerm(:,:,k));
        Pxz(:,:,k)=PPred(:,:,k)*H(:,:,k)';
        PzPred(:,:,k)=H(:,:,k)*PPred(:,:,k)*H(:,:,k)'+PzzHessTerm(:,:,k);
        %Ensure symmetry
        PzPred(:,:,k)=(PzPred(:,:,k)+PzPred(:,:,k)')/2;
    end
end

otherInfo.HJacob=HJacob;
otherInfo.h=h;
otherInfo.Pxz=Pxz;
otherInfo.HHesssian=HHessian;
otherInfo.measPredTrans=measPredTrans;
otherInfo.H=H;
otherInfo.xPred=xPred;
otherInfo.PPred=PPred;

end

function [zPredHessTerm,PzzHessTerm]=getHessianTerms(HHess,PPred)
%The function returns the terms associated with the second derivative in
%the second-order EKF.

zDim=size(HHess,3);
numComp=size(HHess,4);

PzzHessTerm=zeros(zDim,zDim,numComp);
zPredHessTerm=zeros(zDim,numComp);
for k=1:numComp
    for n=1:zDim
        en=zeros(zDim,1);
        en(n)=1;

        HPProdn=HHess(:,:,n,k)*PPred(:,:,k);
        zPredHessTerm=zPredHessTerm+en*trace(HPProdn);
        for m=1:zDim
            em=zeros(zDim,1);
            em(m)=1;
            PzzHessTerm(:,:,k)=PzzHessTerm(:,:,k)+en*em'*trace(HPProdn*HHess(:,:,m,k)'*PPred(:,:,k));
        end
    end
    PzzHessTerm(:,:,k)=PzzHessTerm(:,:,k)/2;
    zPredHessTerm(:,:,k)=zPredHessTerm(:,:,k)/2;
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
