function [zPred,PzPred,otherInfo]=sqrtEKFMeasPred(xPred,SPred,zDim,h,HJacob)
%%SQRTEKFMEASPRED Perform the measurement prediction part of the
%           measurement update step of the first-order Extended Kalman
%           Filter (EKF).  The function sqrtEKFUpdateWithPred can be used
%           to complete the measurement update. Separating the measurement
%           prediction step from the rest of the update step can make the
%           creation of multiple measurement association hypotheses from a
%           single target prediction more efficient. The full measurement
%           update function is sqrtEKFUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states.
%        SPred The xDimXxDimXnumComp lower-triangular square root predicted
%              state covariance matrices for the components in xPred.
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
%
%OUTPUTS: zPred The zDimX1 measurement prediction from the filter.
%        PzPred The zDimXzDim covariance matrix associated with zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to KalmanUpdateWithPred
%               when updating with a measurement.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of 
%[1]. See the comments to sqrtEKFUpdate for more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%sqrtEKFUpdate in one step as with using sqrtEKFMeasPred followed by
%sqrtEKFUpdateWithPred.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% HJacob=@(x)[2*x(1), 2*x(2), 2*x(3), 2*x(4);
%                  1,      0,      0, -(3*sqrt(x(4)))/2];
% SPred=chol([28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13],'lower');
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% SR=eye(zDim,zDim);
% %The update in one step.
% [xUpdate,SUpdate,innov,Szz,W]=sqrtEKFUpdate(xPred,SPred,z,SR,h,HJacob);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=sqrtEKFMeasPred(xPred,SPred,zDim,h,HJacob);
% [xUpdate1,SUpdate1,innov1,Szz1,W1]=sqrtEKFUpdateWithPred(z,SR,zPred,otherInfo);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;SUpdate1(:)-SUpdate(:);innov1(:)-innov;Szz1(:)-Szz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

xDim=size(xPred,1);
numComp=size(xPred,2);
if(isa(HJacob,'function_handle'))
    H=zeros(zDim,xDim,numComp);
    for k=1:numComp
        H(:,:,k)=HJacob(xPred(:,k));
    end
else
    H=HJacob;%If the Jacobian matrix was directly given.
end

zPred=zeros(zDim,numComp);
PzPred=zeros(zDim,zDim,numComp);
Pxz=zeros(xDim,zDim,numComp);
for k=1:numComp
    zPred(:,k)=h(xPred(:,k));
    temp=H(:,:,k)*SPred(:,:,k);
    PzPred(:,:,k)=temp*temp';
    Pxz(:,:,k)=SPred(:,:,k)*SPred(:,:,k)'*H(:,:,k)';
    %Pxz is not needed for the measurement prediction, but we compute it
    %here, so that it need not be recomputed again and again if
    %sqrtEKFUpdateWithPred is called for multiple measurements.
end

otherInfo.H=H;
otherInfo.Pxz=Pxz;
otherInfo.xPred=xPred;
otherInfo.SPred=SPred;
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
