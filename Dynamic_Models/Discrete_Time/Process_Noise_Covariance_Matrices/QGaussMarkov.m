function Q=QGaussMarkov(T,xDim,q,tau,order)
%%QGAUSSMARKOV Get the process noise covariance matrix for an arbitrary
%             order Gauss-Markov process. Gauss-Markov processes have
%             exponentially-correlated noise in their highest moment. A
%             zeroth-order process is an Ornstein-Uhlenbeck model. A
%             first-order Gauss-Markov process is an integrated
%             Ornstein-Uhlenbeck process. A second order process is the
%             Singer dynamic model used in tracking. It is assumed that the
%             a-priori mean of the highest-order moment is zero. This
%             process noise covariance matrix is an exact discretization.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%     xDim The dimensionality of the target state. xDim equals
%          ((order+1)*numDim), where numDim is the number of position
%          dimensions of space (e.g. 2D or 3D).
%        q The power spectral density of the noise corrupting the
%          highest-order moment. q=(sqrt(2/tau)*sigmam)^2 where sigmam^2 is
%          the instantaneous variance of the highest order moment.
%      tau The time constant of the autocorrelation of the moment of the
%          given order. For example, if order=2, then tau is the time
%          constant of the decorrelation time of the acceleration in
%          seconds. The decorrelation time is approximately 2*tau. As tau
%          increases, the highest moment of the process remains correlated
%          longer. A reasonable range for tau when order=2 (Singer's model)
%          is between 5 and 20 seconds. The time constant is assumed the
%          same for all dimensions of motion, so this parameter is scalar.
%          If this parameter is omitted, the default value of 20 is used.
%          The use of large values of tau (typically tau/T ratios of
%          thousands or higher) will lead to finite precision effects that
%          can be sufficiently severe that the matrix ends up with negative
%          diagonal values for the low-order moments. Note that this
%          function will return the correct value if T=0 (which is just an
%          all-zero matrix) even though tau/T=Inf.
%    order The order of the Gauss-Markov process. This is the number of
%          derivatives of position in the model. Thus, 0=position-only,
%          1=position and velocity, etc. If this parameter is omitted, the
%          default value of 2 (Singer's model) is used.
%
%OUTPUTS: Q The process noise covariance matrix under a Gauss-Markov
%           dynamic model of the given order in numDim dimensions where the
%           state is stacked [position;velocity;acceleration;etc].
%
%A closed form-expression for the process noise covariance matrix for an
%order=2 model is in Equation 8.2.3-10 in Chapter 8.2.3 of [1]. In Chapter
%8.2.3 of the aforementioned book also states that a typical value of tau
%for a slowly turning aircraft is 20s, and is 5s for an evasive maneuver.
%The Singer model (order=2) is introduced with solution for tracking in
%[2]. The first-order model, the integrated Ornstein-Uhlenbeck process, is
%discussed for tracking in Chapter 3.2.2 of [3] with a discretization. In
%[4], a general solution is written out.
%
%This function simply generalizes the result to an arbitrary order.
%The dynamic model is written in continuous time in each dimension
%dx/dt=-(1/tau)*x+sqrt(q)*dW where x is the highest order moment. The other
%moments are just integrated. dw is a differential Wiener process. Note
%that in some of the texts, the process noise is expressed in terms of a
%variance such that sigma^2=q*tau/2.
%
%The state transition matrix associated with this model can be
%obtained using the function FGaussMarkov. The continuous-time drift
%function associated with the model is aGaussMarkov and the continuous-time
%diffusion matrix is given by DPoly(1,q,order,numDim) where q is the
%power spectral density of the noise.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] R. A. Singer,"Estimating optimal tracking filter performance for
%    manned maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
%[3] L. D. Stone, C. A. Barlow, and T. L. Corwin, Bayesian Multiple Target
%    Tracking. Boston: Artech House, 1999.
%[4] D. F. Crouse, "The tracker component library: Free routines for
%    rapid prototyping," IEEE Aerospace and Electronic Systems Magazine,
%    vol. 32, no. 5, pp. 18-27, May 2017.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    order=2;
end

if(nargin<4)
    tau=20;
end

if(T==0)
    Q=zeros(xDim,xDim);
    return;
end


numDim=xDim/(order+1);

%First, create the matrix for 1D motion.
Q1=zeros(order+1,order+1);

kSpan=0:(order-1);
innerSumVals=cumsum((T/tau).^kSpan./factorial(kSpan));

eTTau=exp(-T/tau);

constTerm=0.5*(1-exp(-2*T/tau));
%Rows
for r=0:order
    iSpanR=0:(order-r-1);
    rowSum=sum((-1).^iSpanR.*(1-eTTau*innerSumVals(iSpanR+1)));

    for c=0:order
        coeff=q*(-tau)^(2*order-r-c+1);
        iSpanC=0:(order-c-1);
        colSum=sum((-1).^iSpanC.*(1-eTTau*innerSumVals(iSpanC+1)));

        %Q1(r+1,c+1)=rowSum+colSum-constTerm;
        Q1(r+1,c+1)=coeff*(rowSum+colSum-constTerm);
        %Now, for the final term, which is a double sum. This will be done
        %without loops in Matlab.
        [jCGrid,iRGrid]=meshgrid(iSpanC,iSpanR);
        rVals=(-1).^iSpanR.*(T/tau).^iSpanR./factorial(iSpanR);
        cVals=(-1).^iSpanC.*(T/tau).^iSpanC./factorial(iSpanC);

        [CGridVals,rGridVals]=meshgrid(cVals,rVals);

        %Q1(r+1,c+1)=q*(-tau)^(2*order-r-c)*(-tau*Q1(r+1,c+1)+T*sum(sum(rGridVals.*CGridVals./(1+jCGrid+iRGrid))));
        Q1(r+1,c+1)=Q1(r+1,c+1)+q*(-tau)^(2*order-r-c)*T*sum(sum(rGridVals.*CGridVals./(1+jCGrid+iRGrid)));
    end
end

%Note that if order==2, then the following lines are equivalent to the
%above loops and are from Chapter 8.2.3 of [1], but are more subject to
%finite precision errors.
% alpha=1/tau;
% 
% %Equation  8.2.2-7 in [1].
% sigma2=q*tau/2;
% 
% eaT=exp(-alpha*T);
% e2aT=eaT*eaT;%exp(-2*alpha*T);
% 
% %Equation 8.2.3-12
% q11=tau^4*(1-e2aT+2*alpha*T+2*alpha^3*T^3/3-2*alpha^2*T^2-4*alpha*T*eaT);    
% %Equation 8.2.3-13
% q12=tau^3*(e2aT+1-2*eaT+2*alpha*T*eaT-2*alpha*T+alpha^2*T^2);
% %Equation 8.2.3-14
% q13=tau^2*(1-e2aT-2*alpha*T*eaT);
% %Equation 8.2.3-15
% q22=tau^2*(4*eaT-3-e2aT+2*alpha*T);
% %Equation 8.2.3-16
% q23=tau*(e2aT+1-2*eaT);
% %Equation 8.2.3-17
% q33=1-e2aT;
% 
% Q1=sigma2*[q11,q12,q13;
%            q12,q22,q23;
%            q13,q23,q33];

%Then repeat the matrix across multiple dimensions. This is done using a
%Kronecker product.
Q=kron(Q1,eye(numDim));

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
