function Q=QGaussMarkov(T,x,q,tau,order)
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
%INPUTS: T    The time-duration of the propagation interval in seconds.
%        x    The numDim*(order+1)X1 dimensional target state. For the
%             given order, numDim is the dimensionality of the motion. For
%             example, 2D, 3D. The elements of the state are in the order
%             of position, velocity, acceleration, etc. The state itself is
%             not used; this vector is just used to determine numDim. The
%             ability to pass x is for functions that expect the first two
%             parameters of a process noise covariance matrix function to
%             be T and x.
%         q The power spectral density of the noise corrupting the
%           highest-order moment. q=(sqrt(2/tau)*sigmam)^2 where sigmam^2
%           is the instantaneous variance of the highest order moment and
%           numDim is the number of dimensions in the model.
%       tau The time constant of the autocorrelation or the moment of the
%           given order. For example, if order=2, then tau is the
%           time constant of the decorrelation time of the acceleration in
%           seconds. The decorrelation time is approximately 2*tau.
%           As tau increases, the highest moment of the process process
%           remains correlated longer. A reasonable range for tau when
%           order=2 (Singer' model) is between 5 and 20 seconds. The time
%           constant is assumed the same for all dimensions of motion, so
%           this parameter is scalar. If this parameter is omitted, the
%           default value of 20 is used.
%     order The order of the Gauss-Markov process. This is the number of
%           derivatives of position in the model. Thus, 0=position-only,
%           1=position and velocity, etc. If this parameter is omitted, the
%           default value of 2 (Singer's model) is used.
%
%OUTPUTS: Q The process noise covariance matrix under a Gauss-Markov
%           dynamic model of the given order in numDim dimensions where the
%           state is stacked [position;velocity;acceleration;etc].
%
%A closed form-expression for the state transition matrix for an order=2
%model is in Equation 8.2.3-9 in Chapter 8.2.3 of [1]. In Chapter 8.2.3 of
%the aforementioned book also states that a typical value of tau for a
%slowly turning aircraft is 20s, and is 5s for an evasive maneuver.
%The Singer model (order=2) is introduced with solution for tracking in
%[2]. The first-order model, the integrated Ornstein-Uhlenbeck process, is
%discussed for tracking in Chapter 3.2.2 of [3] with a discretization.
%
%This function simply generalizes the redult to an arbitrary order.
%The dynamic model is written in continuous time in each dimension
%dx/dt=-(1/tau)*x+sqrt(q)*dW where x is the highest order moment. The other
%moments are just integrated. dw is a differential Wiener process.
%
%The state transition matrix associated with this model can be
%obtained using the function FGaussMarkov. The continuous-time drift
%function associated with the model is aGaussMarkov and the continuous-time
%diffusion matrix is given by DPoly(x,t,q,order,numDim) where q is the
%power spectral density of the noise.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] R.A.Singer,"Estimating optimal tracking filter performance for manned
%    maneuvering targets," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
%[3] L. D. Stone, C. A. Barlow, and T. L. Corwin, Bayesian Multiple Target
%    Tracking. Boston: Artech House, 1999.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    order=2;
end

if(nargin<4)
    tau=20;
end

xDim=size(x,1);
numDim=xDim/(order+1);

%First, create the matrix for 1D motion.
Q1=zeros(order+1,order+1);

kSpan=0:(order-1);
innerSumVals=cumsum((T/tau).^kSpan./factorial(kSpan));

constTerm=0.5*(1-exp(-2*T/tau));
%Rows
for r=0:order
    iSpanR=0:(order-r-1);
    rowSum=sum((-1).^iSpanR.*(1-exp(-T/tau)*innerSumVals(iSpanR+1)));
    
    for c=0:order
        coeff=q*(-tau)^(2*order-r-c+1);
        iSpanC=0:(order-c-1);
        colSum=sum((-1).^iSpanC.*(1-exp(-T/tau)*innerSumVals(iSpanC+1)));
        
        Q1(r+1,c+1)=coeff*(rowSum+colSum-constTerm);
        %Now, for the final term, which is a double sum. This will be done
        %without loops in matlab.
        [jCGrid,iRGrid]=meshgrid(iSpanC,iSpanR);
        rVals=(-1).^iSpanR.*(T/tau).^iSpanR./factorial(iSpanR);
        cVals=(-1).^iSpanC.*(T/tau).^iSpanC./factorial(iSpanC);
        
        [CGridVals,rGridVals]=meshgrid(cVals,rVals);
        
        Q1(r+1,c+1)=Q1(r+1,c+1)+q*(-tau)^(2*order-r-c)*T*sum(sum(rGridVals.*CGridVals./(1+jCGrid+iRGrid)));
    end
end

if(any(~isfinite(Q1(:))))
    %tau=Inf creates NaN. This is the result of limit as tau->Inf.
    switch order
        case 0
            Q1=q*T;
        case 1
            Q1=q*[T^3/3     T^2/2
                  T^2/2     T];
        case 2
            Q1=q*[T^5/20    T^4/8   T^3/6
                  T^4/8     T^3/3   T^2/2
                  T^3/6     T^2/2   T];
    end
end

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
