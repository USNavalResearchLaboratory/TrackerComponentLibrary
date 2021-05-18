function [zPred,PzPred,otherInfo]=reducedStateMeasPred(xPred,PPred,H)
%%REDUCEDSTATEMEASPRED Perform the measurement prediction part of the
%           measurement update step of the reduced state estimator. The
%           function reducedStateMeasUpdateWithPred can be used to complete
%           the measurement update. Separating the measurement prediction
%           step from the rest of the update step can make the creation of
%           multiple measurement association hypotheses from a single
%           target prediction more efficient. The full measurement update
%           function is reducedStateUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states.
%        PPred The xDimXxDimXnumComp predicted state covariance matrices.
%            H The zDimXxDim measurement matrix. The measurement is modeled
%              as z=H*x+noise.
%
%OUTPUTS: zPred The zDimXnumComp measurement predictions from the filter.
%        PzPred The zDimXzDimXnumComp covariance matrix associated with the
%               values in zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to KalmanUpdateWithPred
%               when updating with a measurement.
%
%The measurement prediction step in the measurement update step of the
%reduced state estimator is the same as that in the Kalman filter. Thus,
%this function just calls KalmanMeasPred. 
%
%The filter is taken from [1]. See the comments to reducedStateUpdate for
%more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%reducedStateUpdate versus calling reducedStateMeasPred and then
%reducedStateMeasUpdateWithPred. 
% xPred=[1e3;-2e3;100;200];
% MPred=[28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13];
% DPred=[8, 0;
%        0, 8;
%        8, 0;
%        0, 8];
% PPred=MPred+(1/2)*(DPred*DPred');
% z=1e3*[-5.498856156296510;
%        1.199241491470584];
% R=eye(2);
% H=[0, 4, 9, 8;
%    6, 3, 0, 6];
% %The update in one step.
% [xUpdate,MUpdate,DUpdate,innov,Pzz,W]=reducedStateUpdate(xPred,PPred,MPred,DPred,z,R,H);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=reducedStateMeasPred(xPred,PPred,H);
% [xUpdate1,MUpdate1,DUpdate1,innov1,Pzz1,W1]=reducedStateMeasUpdateWithPred(z,R,zPred,PzPred,otherInfo,MPred,DPred);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;MUpdate1(:)-MUpdate(:);DUpdate1(:)-DUpdate(:);innov1(:)-innov;Pzz1(:)-Pzz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] P. Mookerjee and F. Reifler, "Reduced state estimator for systems with
%    parametric inputs," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 2, pp. 446-461, Apr. 2004.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[zPred,PzPred,otherInfo]=KalmanMeasPred(xPred,PPred,H);

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
