function [yPred,PInvPred]=discEIFPred(yPrev,PInvPrev,f,FJacob,Q,FHessian)
%%DISCEIFPRED Perform the discrete-time prediction step that comes with 
%             the first and second-order extended information filter (EIF).
%             As the extended information filter as presented in [1] has no
%             real propagation step, this just extracts the state, calls
%             discEKFPred, and converts the result back to an information
%             state. This function is useful if one is using an information
%             state because of how the measurement update in an information
%             filter is performed.
%
%INPUTS:    yPrev   The xDimX1 information state at the previous time-step.
%                   The information state is the inverse covariance matrix
%                   times the target state.
%        PInvPrev   The xDimXxDim inverse of the state covariance matrix at
%                   the previous time-step.
%               f   A function handle for the state transition function
%                   that takes the state as its parameter.
%           FJacob  A function handle for calculating the xDim X xDim
%                   Jacobian of f, or the xDim X xDim Jacobian matrix
%                   itself. If an empty matrix is passed, then
%                   FJacob will be found using numerical differentiation 
%                   via the numDiff function with default parameters.
%               Q   The xDimX xDim process noise covariance matrix.
%         FHessian  This parameter is only provided if a second-order EKF
%                   is desired. This is either a function handle for the
%                   state transition Hessian hypermatrix, or it is the
%                   state transition Hessian hypermatrix itself. The matrix
%                   is xDim X xDim X xDim. The matrix FH=FHessian(x) is
%                   such that FH(i,j,k) is the second derivative of the
%                   kth element of the vector returned by f with
%                   respect to the ith and jth components of x. The Hessian
%                   matrix is symmetric. If this parameter is omitted, a
%                   first-order filter is used.
%
%OUTPUTS:   yPred    The xDim X 1 predicted information state vector.
%           PInvPred The predicted xDim X xDim inverse state covariance
%                    matrix.
%
%REFERENCES:
%[1] K. P. B. Chandra, D.-W. Gu, and I. Postlethwaite, "Square root
%    cubature information filter," IEEE Sensors Journal, vol. 13, no. 2,
%    pp. 750-758, Feb. 2013.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6)
    FHessian=[];
end

PPrev=pinv(PInvPrev);
xPrev=PPrev*yPrev;

[xPred, PPred]=discEKFPred(xPrev,PPrev,f,FJacob,Q,FHessian);

PInvPred=pinv(PPred);
yPred=PInvPred*xPred;

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
