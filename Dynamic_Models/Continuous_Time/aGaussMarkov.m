function [aVal,aJacob,aHess,papt]=aGaussMarkov(x,tau,order)
%AGAUSSMARKOV The drift function for an arbitrary-order Gauss Markov
%             process. Gauss-Markov processes have exponentially-correlated
%             noise in their highest moment. A zeroth-order process is an
%             Ornstein-Uhlenbeck model. A first-order Gauss-Markov process
%             is an integrated Ornstein-Uhlenbeck process. A second order
%             process is the Singer dynamic model used in tracking. It is
%             assumed that the a-priori mean of the highest-order moment is
%             zero.
%
%INPUTS: x The numDim*(order+1)X1 (aka xDimX1) dimensional target state.
%          For the given order, numDim is the dimensionality of the motion.
%          For example, 2D, 3D. The elements of the state are in the order
%          of position, velocity, acceleration, etc.
%      tau The time constant of the autocorrelation or the moment of the
%          given order. For example, if order=2, then tau is the time
%          constant of the decorrelation time of the acceleration in
%          seconds. The decorrelation time is approximately 2*tau. As tau
%          increases, the highest moment of the process process remains
%          correlated longer. A reasonable range for tau when order=2
%          (Singer's model) is between 5 and 20 seconds. The time constant
%          is assumed the same for all dimensions of motion, so this
%          parameter is scalar. If this parameter is omitted, the default
%          value of 20 is used.
%    order The order of the Gauss-Markov process. This is the number of
%          derivatives of position in the model. Thus, 0=position-only,
%          1=position and velocity, etc. If this parameter is omitted, the
%          default value of 2 (Singer's model) is used.
%
%OUTPUTS: aVal The time-derivative of the state under the Gauss-Markov model
%              having the given order and whose dimensionality is inferred
%              from x, whose elements are ordered
%              [position;velocity;acceleration;etc].
%       aJacob The xDimXxDim Jacobian of aVal. This is such that
%              aJacob(:,k) is the partial derivative of aVal with respect
%              to the kth element of x.
%        aHess The xDimXxDimXxDim hypermatrix such that aHess(:,k1,k2) is
%              the second partial derivative of aVal with respect to
%              elements k1 and k2 of x (all zero in this instance).
%         papt The xDimX1 partial derivative of aVal with respect to t (all
%              zero in this instance).
%
%The Singer model (order=2) is introduced for tracking in [1]. The first-
%order model, the integrated Ornstein-Uhlenbeck process, is discussed for
%tracking in Chapter 3.2.2 of [2]. In Chapter 8.2.3 of [3] it is stated
%that a typical value of tau for a slowly turning aircraft is 20s, and is
%5s for an evasive maneuver.
%
%This function simply generalized the concept from first and second orders
%to arbitrary orders.
%
%The dynamic model is written in each dimension
%dx/dt=-(1/tau)*x+sqrt(q)*dW where x is the highest order moment. The other
%moments are just integrated. dW is a differential Wiener process.
%
%The drift function corresponds to the state transition given in
%discrete-time by the function FGaussMarkov. To get the full dynamic model,
%the drift function should be used with the diffusion function
%DPoly(1,q,order,numDim), where q=(sqrt(2/tau)*sigmam)^2 and sigmam^2 is
%the instantaneous variance of the highest order moment and  numDim is
%the number of dimensions in the model. The autocorrelation of the
%highest-order moment is= sigmam^2*exp(-abs(deltaT)/tau) for some deltaT
%time interval. These relations come from Ch. 8.2 of [3], after relating
%the alpha term there to tau. The discrete-time process noise covariance
%matrix is given by QGaussMarkov.
%
%REFERENCES:
%[1] R. A. Singer,"Estimating optimal tracking filter performance for
%    manned maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
%[2] L. D. Stone, C. A. Barlow, and T. L. Corwin, Bayesian Multiple Target
%    Tracking. Boston: Artech House, 1999.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    order=2;
end

if(nargin<2)
	tau=20; 
end

xDim=size(x,1);
numDim=xDim/(order+1);

%Allocate space
aVal=zeros(xDim,1);

%All but the highest-order moment are just integrated.
aVal(1:(numDim*order))=x((numDim+1):end);

%The highest-order moment is set to -1/tau times itself. 
aVal((numDim*order+1):end)=-1/tau*x((numDim*order+1):end);

if(nargout>1)
    aJacob=[zeros(xDim,numDim),[eye(xDim-numDim,xDim-numDim);zeros(numDim,xDim-numDim)]];

    if(nargout>2)
        aHess=zeros(xDim,xDim,xDim);

        if(nargout>3)
            papt=zeros(xDim,1);
        end
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
