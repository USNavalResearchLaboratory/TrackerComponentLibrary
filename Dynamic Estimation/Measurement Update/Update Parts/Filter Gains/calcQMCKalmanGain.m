function W=calcQMCKalmanGain(R,PzPred,otherInfo)
%%CALCQMCKALMANGAIN Using the output of the quasi-Monte Carlo Kalman filter
%           measurement prediction function QMCKalMeasPred and the
%           covariance matrix of a measurement R, compute the gain of a
%           quasi-Monte Carlo Kalman filter. This function can be useful
%           for creating a gain based on some type of "maximum" covariance
%           matrix for a scene for purposes of giving it to the
%           calcMissedGateCov function to approximate how the covariance
%           matrix of the target state prediction should be increased for a
%           missed-detection event when gating is performed.
%
%INPUTS: R The zDimXzDim measurement covariance matrix.
%   PzPred The zDimXzDim covariance matrix of the measurement predicted by
%          the filter. 
% otherInfo The structure returned by the QMCKalMeasPred function that
%          includes various terms that can be reused.
%
%OUTPUTS: W The xDimXzDim gain matrix used in the filter.
%
%The function just calls calcCubKalGain to get the gain. This produces the
%gain in the same manner as is done in QMCKalUpdateWithPred and
%QMCKalUpdate. See the comments to those functions for more information.
%
%EXAMPLE:
%Here, we demonstrate that the gain used in a complete update is the same
%as that obtained using just the prediction via QMCKalMeasPred and then
%calcQMCKalmanGain for the gain without using the measurement z itself.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% PPred=[28,   3.5,    6,  8.5;
%       3.5,    23,  8.5,   11;
%         6,   8.5,   18, 13.5;
%       8.5,    11, 13.5,   13];
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% R=eye(zDim,zDim);%Measurement covariance matrix.
% %We use the same random number generator seed for both approaches so that
% %the random part of the results is identical.
% rng('default');
% rng(1);
% [~,~,~,~,W]=QMCKalUpdate(xPred,PPred,z,R,h);
% %Reset the random number generator.
% rng(1);
% [~,PzPred,otherInfo]=QMCKalMeasPred(xPred,PPred,zDim,h);
% W1=calcQMCKalmanGain(R,PzPred,otherInfo);
% %One will see that the result below is true (1).
% all(W1(:)==W(:))
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

W=calcCubKalGain(R,PzPred,otherInfo);

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
