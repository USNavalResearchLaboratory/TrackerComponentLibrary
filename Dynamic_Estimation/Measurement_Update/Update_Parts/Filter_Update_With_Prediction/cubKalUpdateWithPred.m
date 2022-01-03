function [xUpdate,PUpdate,innov,Pzz,W]=cubKalUpdateWithPred(z,R,zPred,PzPred,otherInfo)
%%CUBKALUPDATEWITHPRED Given the output of the measurement prediction step
%           from cubKalMeasPred and a measurement, complete the measurement
%           update step of the cubature Kalman filter with additive
%           measurement noise. Separating the measurement prediction step
%           from the rest of the update step can make the creation of
%           multiple measurement association hypotheses from a single
%           target prediction more efficient. The full measurement update
%           function is cubKalUpdate.
%
%INPUTS:  z The zDimX1 measurement vector.
%         R The zDimXzDim measurement covariance matrix in the native
%           coordinate system of the measurement.
%     zPred The zDimXnumComp measurement predictions from the filter.
%    PzPred The zDimXzDimXnumComp covariance matrices associated with
%           zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%           the cubKalMeasPred function.
%
%OUTPUTS: xUpdate The xDimXnumComp updated state vectors.
%         PUpdate The updated xDimXxDimXnumComp state covariance matrices.
%      innov, Pzz The zDimXnumComp innovations and the zDimXzDimXnumComp
%                 innovation covariance matrices are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The xDimXzDimXnumComp gains used in the update.
%
%See the comments to the function cubKalMeasPred for an example of usage of
%this function. See the comments to cubKalUpdate for more information on
%the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    innovTrans=otherInfo.innovTrans;
    stateDiffTrans=otherInfo.stateDiffTrans;
    stateTrans=otherInfo.stateTrans;
    xPredCenPoints=otherInfo.xPredCenPoints;
    zPredCenPoints=otherInfo.zPredCenPoints;
    xPred=otherInfo.xPred;
    w=otherInfo.w;
    Pxz=otherInfo.Pxz;
    
    numCubPoints=length(w);
    xDim=size(xPred,1);
    numComp=size(xPred,2);
    zDim=size(z,1);

    xUpdate=zeros(xDim,numComp);
    PUpdate=zeros(xDim,xDim,numComp);
    innov=zeros(zDim,numComp);
    Pzz=zeros(zDim,zDim,numComp);
    W=zeros(xDim,zDim,numComp);

    for k=1:numComp
        Pzz(:,:,k)=PzPred(:,:,k)+R;

        %The innovation, transformed as necessary to keep values in a
        %desired range.
        innov(:,k)=innovTrans(z,zPred(:,k));

        %The filter gain
        W(:,:,k)=Pxz(:,:,k)/Pzz(:,:,k);
    
        %Updated state estimate
        xUpdate(:,k)=stateTrans(xPred(:,k)+W(:,:,k)*innov(:,k));
    
        %Updated state covariance matrix
        %We could just do a simple one-line solution as in [1] and [2].
        %However, that does not guarantee that PUpdate will always be
        %positive (semi)definite. Thus, we use an equivalent but more
        %complicated update formula based on Equation in Appendix C of [1].
        PUpdate(:,:,k)=W(:,:,k)*R*W(:,:,k)';
        for curP=1:numCubPoints
            diff=stateDiffTrans(xPredCenPoints(:,curP,k)-W(:,:,k)*zPredCenPoints(:,curP,k));
            PUpdate(:,:,k)=PUpdate(:,:,k)+w(curP)*(diff*diff');
        end
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
