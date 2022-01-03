function [zPred,PzPred,otherInfo]=KalmanMeasPred(xPred,PPred,H)
%%KALMANMEASPRED Perform the measurement prediction part of the measurement
%           update step of the Kalman filter. The function
%           KalmanUpdateWithPred can be used to complete the measurement
%           update. Separating the measurement prediction step from the
%           rest of the update step can make the creation of multiple
%           measurement association hypotheses from a single target
%           prediction more efficient. The full measurement update function
%           is KalmanUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states. If one only wants
%              PzPred and otherInfo, an empty matrix can be passed and
%              zPred will be returned as an empty matrix.
%        PPred The xDimXxDimXnumComp predicted state covariance matrices.
%            H The zDimXxDim measurement matrix. The measurement is modeled
%              as z=H*x+noise and the model is the same for all components
%              given in xPred and PPred.
%
%OUTPUTS: zPred The zDimXnumComp measurement prediction from the filter or
%               an empty matrix if xPred was an empty matrix.
%        PzPred The zDimXzDimXnumComp covariance matrix associated with
%               zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to KalmanUpdateWithPred
%               when updating with a measurement.
%
%The Kalman filter is derived in Chapter 5 of [1]. See the comments to
%KalmanUpdate for more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%KalmanUpdate versus calling KalmanMeasPred and then KalmanUpdateWithPred. 
% xPred=[1e3;-2e3;100;200];
% PPred=[28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13];
% z=1e3*[-5.498856156296510;
%        1.199241491470584];
% R=eye(2);
% H=[0, 4, 9, 8;
%    6, 3, 0, 6];
% %The update in one step.
% [xUpdate,PUpdate,innov,Pzz,W]=KalmanUpdate(xPred,PPred,z,R,H);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=KalmanMeasPred(xPred,PPred,H);
% [xUpdate1,PUpdate1,innov1,Pzz1,W1]=KalmanUpdateWithPred(z,R,zPred,PzPred,otherInfo);
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

if(~isempty(xPred))
    zPred=H*xPred;
else
    zPred=[];
end

zDim=size(H,1);
xDim=size(H,2);
numComp=size(PPred,3);
Pxz=zeros(xDim,zDim,numComp);
PzPred=zeros(zDim,zDim,numComp);

for k=1:numComp
    Pxz(:,:,k)=PPred(:,:,k)*H';
    PzPredCur=H*PPred(:,:,k)*H';
    %Ensure symmetry.
    PzPred(:,:,k)=(PzPredCur+PzPredCur')/2;
end

otherInfo.xPred=xPred;
otherInfo.Pxz=Pxz;
otherInfo.PPred=PPred;
otherInfo.H=H;

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
