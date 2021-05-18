function [yPred, PInvPred]=infoFilterDiscPred(yPrev,PInvPrev,F,Q,u)
%%INFOFILTERDISCPRED Perform the discrete-time prediction step that comes 
%                    with the standard linear information filter with
%                    additive process noise.
%
%INPUTS: yPrev The xDimX1 information state at the previous time-step. The
%              information state is the inverse covariance matrix times the
%              target state.
%     PInvPrev The xDimXxDim inverse of the state covariance matrix at the
%              previous time-step.
%            F An xDim X xDim state transition matrix.
%            Q The xDimX xDim process noise covariance matrix. This can be
%              singular.
%            u An optional xDimX1 vector that is the control input. If
%              omitted, no control input is used.
%
%OUTPUTS: yPred The xDim X 1 predicted information state vector.
%      PInvPred The predicted xDim X xDim inverse state covariance matrix.
%
%The implementation of the prediction step given here is from the flow
%chart given in [1]. The matrix G in the paper is omitted, since any
%control input can be pre-multipled by the matrix. The implied discrete-
%time dynamic model is
%x(k)=F(k-1)*x(k-1)+u(k-1)+noise
%where x(k) is the state at time k. However, an information filter
%propagates
%y=inv(P)*x
%and 
%PInv=inv(P)
%instead of propagating P and x, where P is the covariance matrix
%associated with the Gaussian state x. This allows the filter to be used
%even when PInv is singular.
%
%More information on information filtering is given in Chapter 7.2 of [2].
%
%REFERENCES:
%[1] D. F. Crouse, P. Willett, and Y. Bar-Shalom, "A low-complexity
%    sliding-window Kalman FIR smoother for discrete-time models," IEEE
%    Signal Processing Letters, vol. 17, no. 2, pp. 177-180, Feb. 2009.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    DInv=F'+PInvPrev/(F)*Q;
    PInvPred=DInv\PInvPrev/(F);
    
    %Handle possible loss of symmetry due to order of operations and finite
    %precision limitations.
    PInvPred=(PInvPred+PInvPred')/2;

    if(nargin<5||isempty(u))
        yPred=DInv\yPrev;
    else
        yPred=DInv\yPrev+PInvPred*u;
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
