function [xHat, P]=DiscPriorPModel(k,xInit,param3,param4,mode)
%%DISCMODELPPRIOR  Obtain the mean and covariance of the prior
%                  distribution of a target under a discrete-time linear
%                  motion model at discrete-time k where at time k=0
%                  the state is given by a delta function at xInit.
%                  Standard polynomial motion models in 3 Cartesian
%                  dimensions are provided, or parameters for an arbitrary
%                  model can be specified.
%
%INPUTS: k The discrete-time-step at which the prior distribution is
%          desired.
%    xInit The xDimX1 target state at time k=0. The PDF at this time is
%          modeled as
%          a delta function. If the mode argument of this function is
%          omitted or is set to 0 then the dimensionality of xInit
%          determines the order of the model used. A 3X1 vector means a
%          stationary target model; a 6X1 vector means a nearly constant
%          velocity model (the discretized continuous white noise
%          acceleration model) and a 9X1 vector means a nearly constant
%          acceleration model (the discretized continuous white noise jerk
%          model).
%   param3 If mode=0, then param3 is T, the time interval between discrete
%          time steps. Otherwise, it is F, the state transition matrix
%          under a linear, discrete-time model.
%   param4 If mode=0, then param 4 is q0, the power spectral density of the
%          process noise. Otherwise, it is Q, the process noise covariance
%          matrix under a linear, discrete-time model.
%     mode If mode=0 or the mode parameter is omitted, then a polynomial
%          linear motion based on the dimensionality of xInit is used.
%          Otherwise, for mode=1, the given transition matrix and state
%          covariance matrix in param3 and param4 for an arbitrary
%          discrete-time model are used.
%
%OUTPUTS: xHat The xDimX1 mean of the distribution after the specified
%              number of steps.
%            P The xDimXxDim covariance matrix of the distribution after
%              the specified number of steps.
%
%The first and second moments of the prior (i.e. measurement-free)
%distribution that are computed are useful for evaluating the posterior
%Cramér-Rao bound (PCRLB) in simulations. The moments are computed at
%time-step k, where at time step 0, the distribution is given by a delta
%function at xInit. The assumes dynamic model has the form
%x(k)=F*x(k-1)+noise
%where noise has covariance matrix Q.
%
%The computation of a prior distribution of a target given its state at
%time k=0 as a delta function is discussed in Section VI of [1].
%
%The polynomial motion models are derived in Chapter 6.2 of [2].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%October 2013 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(mode))
   mode=0; 
end

if(mode~=0)
    F=param3;
    Q=param4;
    
    xHat=F^k*xInit;
    xDim=size(xInit,1);
    P=zeros(xDim,xDim);
    for n=0:(k-1)
        FPow=F^n;
        P=P+FPow*Q*FPow';
    end
    %Ensure symmetry is preserved.
    P=(P+P')/2;
    return;
end

T=param3;
q0=param4;

switch(length(xInit))
    case 3
        FSum=1;
        qt11=k*T*q0;
        P=blkdiag(qt11,qt11,qt11);
    case 9
        [FSum,P]=getFQSumDCWJA(k,T,q0);
    otherwise
        [FSum,P]=getFQSumDCWNA(k,T,q0);
end
%Ensure symmetry is preserved.
P=(P+P')/2;
xHat=FSum*xInit;
end
    
function [FSum,QSum]=getFQSumDCWNA(k,T,q0)
    %For the discretized continuous white noise acceleration model.
    FSumTilde=[1,k*T;
               0 1];

    qt11=k^3*T^3/3;
    qt12=k^2*T^2/2;
    qt22=k*T;

    QTilde=[qt11,   qt12;
            qt12,   qt22];
    
    FSum=blkdiag(FSumTilde,FSumTilde,FSumTilde);
    QSum=blkdiag(QTilde,QTilde,QTilde)*q0;
%The ordering of the elements in FSum and QSum now corresponds to a state
%of [x;xDot;y;yDot;z;zDot]. The elements shall be rearranged to make it
%correspond to a state of [x;y;z;xDot;yDot;zDot]
    sel=[1;3;5;2;4;6];
    FSum=FSum(sel,sel);
    QSum=QSum(sel,sel);
end

function [FSum,QSum]=getFQSumDCWJA(k,T,q0)
    %For the discretized continuous white noise jerk model.
    FSumTilde=[1,   k*T,    k^2*T^2/2;
               0,   1,      k*T;
               0,   0,      1];

    qt11=k^5*T^5/20;
    qt12=k^4*T^4/8;
    qt13=k^3*T^3/6;
    qt22=k^3*T^3/3;
    qt23=k^2*T^2/2;
    qt33=k*T;

    QTilde=[qt11,   qt12,   qt13;
            qt12,   qt22,   qt23;
            qt13,   qt23,   qt33];
    
    FSum=blkdiag(FSumTilde,FSumTilde,FSumTilde);
    QSum=blkdiag(QTilde,QTilde,QTilde)*q0;
%The ordering of the elements in FSum and QSum now corresponds to a state
%of [x;xDot;xDotDot;y;yDot;yDotDot;z;zDot;zDotDot]. The elements shall be
%rearranged to make it correspond to a state of
%[x;y;z;xDot;yDot;zDot;xDotDot;yDotDot;zDotDot]
    sel=[1;4;7;2;5;8;3;6;9];
    FSum=FSum(sel,sel);
    QSum=QSum(sel,sel);
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
