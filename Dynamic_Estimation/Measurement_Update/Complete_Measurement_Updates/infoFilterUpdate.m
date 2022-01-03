function [yUpdate, PInvUpdate]=infoFilterUpdate(yPred,PInvPred,z,R,H)
%%INFOFILTERUPDATE Perform the measurement update step in the standard 
%                  linear information filter. Using an information filter
%                  means that instead of propagating a state x and its
%                  covariance matrix P, one propagates inv(P) and 
%                  y=inv(P)*x.
%
%INPUTS: yPred The xDimX1 predicted information state. The information
%              state is the inverse covariance matrix times the target
%              state.
%     PInvPred The xDimXxDim inverse of the predicted state covariance
%              matrix.
%            z The zDimX1 vector measurement.
%            R The zDimXzDim measurement covariance matrix. This must be
%              positive definite.
%            H The zDimXxDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%
%OUTPUTS: yUpdate The xDimX1 updated (posterior) information state vector.
%      PInvUpdate The updated xDimXxDim inverse state covariance matrix.
%
%The information filter is algebraically equivalent to the standard linear
%Kalman filter, but allows for track propagation with very uncertain
%states, such as when starting tracks. The implementation of the update
%given here is from the flow chart given in Appendix H of [1].
%
%The standard information filter has the implied linear measurement
%equation
%z=H*x+w
%where z is the measurement, w is the zero-mean Gaussian measurement noise
%with covariance matrix R and x is the true state. However, instead of
%propagating the true state x with its covariance matrix P, the information
%filter propagates
%y=inv(P)*x
%and instead of propagating the true covariance matrix P, it propagates
%PInv=inv(P)
%which means that the filter can be used even when PInv is singular.
%
%More information on information filtering is given in Chapter 7.2 of [2].
%
%EXAMPLE:
%In this example, we feed two measurements to the information filter and
%show that the result has a consistent NEES and the RMSE is the same as the
%second measurement.
% T=1;%Sample period.
% numRuns=1000;
% %Parameters for the noise process dynamics and measurement.
% H=[eye(2,2),zeros(2,2)];
% zDim=size(H,1);
% xDim=size(H,2);
% R=diag([40;40]);
% SR=chol(R,'lower');
% 
% %Statistics for the initial state.
% x0Mean=[1e3;0;75;-50];
% x0Cov=diag([1e3^2;100^2;25^2;25^2]);
% x0S=chol(x0Cov,'lower');
% 
% %Parameters for the state dynamics.
% q=processNoiseSuggest('PolyKal-ROT',9.8,1);
% F=FPolyKal(T,xDim,1);
% Q=QPolyKal(T,xDim,1,q);
% SQ=chol(Q,'lower');
% 
% RMSEMeas=0;
% RMSE=0;
% NEES=0;
% for curRun=1:numRuns
%     %Draw the initial state.
%     x1True=x0Mean+x0S*randn(xDim,1);
%     %Get the first measurement.
%     z1=H*x1True+SR*randn(zDim,1);
% 
%     x2True=F*x1True+SQ*randn(xDim,1);
%     z2=H*x2True+SR*randn(zDim,1);
%     
%     %Two-point initialization.
%     y0=zeros(xDim,1);
%     PInv0=zeros(xDim,xDim);
%     [yUpdate,PInvUpdate]=infoFilterUpdate(y0,PInv0,z1,R,H);
%     [yPred, PInvPred]=infoFilterDiscPred(yUpdate,PInvUpdate,F,Q);
%     [yUpdate,PInvUpdate]=infoFilterUpdate(yPred,PInvPred,z2,R,H);
%     xEst2=PInvUpdate\yUpdate;
%     
%     diff=xEst2-x2True;
%     NEES=NEES+diff'*PInvUpdate*diff;
%     RMSE=RMSE+sum(diff(1:2).^2);
%     diff=z2-x2True(1:2);
%     RMSEMeas=RMSEMeas+sum(diff(1:2).^2);
% end
% RMSEMeas=sqrt(RMSEMeas/numRuns)
% RMSE=sqrt(RMSE/numRuns)
% NEES=NEES/(xDim*numRuns)
%One will see that RMSEMeas equal RMSE, because with just two measurements,
%one cannot smmooth the position estimate syet. Additionally, NEES will be
%close to 1 indicating covariance consistency.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    yUpdate=yPred+H'/R*z;
    PInvUpdate=PInvPred+H'/R*H;
    %Ensure symmetry
    PInvUpdate=(PInvUpdate+PInvUpdate')/2;
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
