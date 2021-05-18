function [xUpdate,PUpdate,innov,W]=HInfinityUpdate(xPred,PPred,z,RInv,H,gammaVal,QBar)
%%HINFINITYUPDATE Perform the measurement update step of the finite-horizon
%                 H-infinity-filter. This filter is similar to the Kalman
%                 filter, but is more robust to model mismatch errors. Note
%                 that the prediction step of a completely linear,
%                 discrete-time H-infinity filter is the same as that of
%                 the Kalman filter, and thus the function discKalPred can
%                 be used for prediction.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        PPred The xDim X xDim predicted state covariance matrix. This must
%              be positive definite.
%            z The zDim X 1 vector measurement.
%         RInv The zDim X zDim inverse of the measurement covariance
%              matrix. 
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%     gammaVal The inverse of the positive scalar bound on the cost
%              function as given in Equation 4 of [1]. If this is zero
%              (unbounded), then the update is equavalent to the Kalman
%              filter update. This can be seen by comparing the update
%              equation with those of the information filter in Chapter
%              7.2.2 of [2].
%         QBar A symmetric matrix that goes into the cost function. If one
%              is not interested in the state x itself, but is actually
%              interested in L*x, where L is some matrix, then
%              QBar=L'*QT*L, where QT is a matrix that goes into the cost
%              function. If this parameter is omitted or an empty matrix is
%              passed, then the identity matrix will be used.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%           innov The zDimX1 innovation of the filter. This is sometimes
%                 used in gating and likelihood evaluation.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%The H-infinity filter optimizes over a cost function given in Equation 3
%in [1]. The numerator of the cost function includes the QT term and is
%sum((l_k-\hat{l}_{k})*inv(QT)'*(l_k-\hat{l}_{k}))/
%where l_k=L*x_k, where x_k is the state but we are actually interested in
%the linear combination of the state given by l and \hat{l}_{k} is the
%estimate of the state. The QT matrix affects the type of penalty. The
%default using QT as the identity matrix is equivalent to saying that we
%are interested in all components of the state with equal weighting.
%
%In Section 7.2.2 of [2], equations for the measurement update step of the
%information filter are given. It can be seen that the formulation of the
%H-infinity filter in [1] is the same as that of the information filter,
%with the exception that an extra term is present reducting the information
%available (inflating the covariance matrix).
%
%The pair QT and gammaVal are design parameters, but the proper of choice
%of them can be problematic. For a fixed QT, if gammaVal is chosen to be
%too large, then the inverse covariance matrix (which must be inverted to
%get PUpdate) will have negative eigenvalues. In systems where the sampling
%rate is constant, then choosing the correct pair is not particularly
%difficult (choose a QT and then adjust gammaVal for performance
%constraints, such that the filter works). However, given a variable
%sampling rate, the situation is more difficult.
%
%Note that if QBar can be written in the form QBar=H*QT*H', then this
%filter is equivalent to the standard Kalman filter just using a covariance
%matrix of R=inv(inv(R)-gammaVal*QT). In the more general case, this filter
%is equivalent to inflating the Kalman filter predicted covariance value
%by as PPred=inv(inv(PPred)-gammaVal*QBar) prior to performing a
%measurement update. In some instances, this can be mathematically
%equivalent to just inflating the process noise covariance matrix and the
%measurement covariance matrix.
%
%REFERENCES:
%[1] F. Yang, Z. Wang, S. Lauria, and X. Liu, "Mobile robot localization
%    using robust extended H-infinity filtering," Journal of Systems and
%    Control Engineering, vol. 223, no. 8, pp. 1067-1080, 1 Dec. 2008.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(QBar))
%If QBar is not given, then just make the information reduction from the
%gamma term affect all components of the state equally.
    xDim=size(xPred,1);
    QBar=eye(xDim,xDim);
end

%H-infinity covariance update, Equation 7 in [1] (without the prediction
%part), modified slightly to look more like the information filter
%covariance update in Section 7.2.2 of [2]. That is, for gammaVal=0, this
%is the same as the measurement update step of the information filter
%(Kalman filter). If gammaVal*QBar becomes too large, the inverse will not
%exist.
PUpdateInv=inv(PPred)-gammaVal*QBar+H'*RInv*H;
PUpdate=inv(PUpdateInv);
%Ensure symmetry
PUpdate=(PUpdate+PUpdate')/2;

%H-infinity filter gain, Equation 10 in [1] (without the prediction part),
%modified slightly to look more like the information filter gain in Section
%7.2.2 of [2]. That is, for gammaVal=0, this is the same as the gain in the
%information filter (Kalman filter).
W=PUpdateInv\H'*RInv;

%Predicted measurement
zPred=H*xPred;

%Innnovation
innov=z-zPred;

%H-infinity update, Equation 9 in [1] (without the prediction part).
xUpdate=xPred+W*innov;

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
