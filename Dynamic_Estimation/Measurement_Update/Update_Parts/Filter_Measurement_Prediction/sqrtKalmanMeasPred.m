function [zPred,PzPred,otherInfo]=sqrtKalmanMeasPred(xPred,SPred,H)
%%SQRTKALMANMEASPRED Perform the measurement prediction part of the
%           measurement update step of the square-root Kalman filter. The
%           function sqrtKalmanUpdateWithPred can be used to complete the
%           measurement update. Separating the measurement prediction step
%           from the rest of the update step can make the creation of
%           multiple measurement association hypotheses from a single
%           target prediction more efficient. The full measurement update
%           function is sqrtKalmanUpdate.
%
%INPUTS: xPred The xDimXnumComp set of numComp predicted target states.
%        SPred The xDimXxDimXnumComp lower-triangular predicted square-root
%              state covariance matrices.
%            H The zDimXxDim measurement matrix. The measurement is modeled
%              as z=H*x+noise and the model is the same for all components
%              given in xPred and SPred.
%
%OUTPUTS: zPred The zDimXnumComp measurement predictions from the filter.
%        PzPred The zDimXzDimXnumComp covariance matrix associated with
%               zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to sqrtKalmanUpdate
%               when updating with a measurement.
%
%The mathematics behind the specific square root implementation used here
%are described in Appendix G of [1].
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%KalmanUpdate versus calling sqrtKalmanMeasPred and then
%sqrtKalmanUpdateWithPred. 
% xPred=[1e3;-2e3;100;200];
% SPred=chol([28,   3.5,    6,  8.5;
%            3.5,    23,  8.5,   11;
%              6,   8.5,   18, 13.5;
%            8.5,    11, 13.5,   13],'lower');
% z=1e3*[-5.498856156296510;
%        1.199241491470584];
% SR=eye(2);
% H=[0, 4, 9, 8;
%    6, 3, 0, 6];
% 
% %The update in one step.
% [xUpdate,SUpdate,innov,Szz,W]=sqrtKalmanUpdate(xPred,SPred,z,SR,H);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=sqrtKalmanMeasPred(xPred,SPred,H);
% [xUpdate1,SUpdate1,innov1,Szz1,W1]=sqrtKalmanUpdateWithPred(z,SR,zPred,otherInfo);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;SUpdate1(:)-SUpdate(:);innov1(:)-innov;Szz1(:)-Szz(:);W1(:)-W(:)]))
%
%REFERENCES
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zDim=size(H,1);
xDim=size(H,2);
numComp=size(xPred,2);

zPred=H*xPred;

PzPred=zeros(zDim,zDim,numComp);
Pxz=zeros(xDim,zDim,numComp);
for k=1:numComp
    temp=H*SPred(:,:,k);
    PzPred(:,:,k)=temp*temp';

    Pxz(:,:,k)=SPred(:,:,k)*SPred(:,:,k)'*H';
    %Pxz is not needed for the measurement prediction, but we compute it
    %here, so that it need not be recomputed again and again if
    %sqrtKalmanUpdateWithPred is called for multiple measurements.
end
otherInfo.xPred=xPred;
otherInfo.SPred=SPred;
otherInfo.Pxz=Pxz;
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
