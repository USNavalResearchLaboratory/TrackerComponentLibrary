function CRLB=delayDopplerRateCRLBLFMApprox(B,T,P,TR,SNRIn,returnFIM,c,fc,useHalfRange)
%DELAYDOPPLERRATECRLBLFMAPPROX Determine an approximate Cramér-Rao lower
%      bound matrix for estimating the delay and Doppler rate when
%      receiving a chain of identical linearly frequency modulated
%      (LFM) chirps. Whether they are up or down chirps is
%      determined by the sign of the bandwidth. This uses a narrowband
%      approximation.
%
%INPUTS: B The _signed_ bandwidth of the chirp. If this is positive, then
%          this is an upchirp. If this is negative, then this is a
%          downchirp.
%        T The duration of the chirp in each pulse, usually in seconds.
%        P The integer number of pulses in the train; P>1.
%       TR The duration of each pulse in the pulse train, usually in
%          seconds. TR>=T.
%    SNRIn The SNR of the fully integrated signal (NOT in dB). How this
%          relates to other parameters is discussed below. This is a power
%          SNR (not a voltage SNR).
% returnFIM If this is true, the Fisher information matrix (FIM) is
%          returned instead of the CRLB. The default if omitted or an empty
%          matrix is passed is false.
%    c, fc These two optional parameters are only provided if one wishes to
%          have the results given in terms of round trip range and
%          approximate range rate. c is the propagation speed of the medium
%          (e.g. Constants.speedOfLight) and fc is the center frequency of
%          the transmission. The range rate conversion is only approximate
%          as a more precise conversion requires the true range rate as can
%          be inferred from [2]. If one doesn't wish to convert the units
%          of the output, these can be omitted or empty matrices can be
%          passed.
% useHalfRange If and fc and c are provided, so that the CRLB is for
%          outputs are range and approximate range rate instead of delay
%          and Doppler shift, then if this is true, the outputs are divided
%          by 2. This can be useful when dealing with a monostatic system
%          where one wants the one-way range and range rate. The default if
%          omitted or an empty matrix is passed is false.
%
%OUTPUTS: CRLB The 2X2 Cramér-Rao lower bound for a vector [delay;Doppler
%              rate] (or the FIM if returnFIM is true. The delay is
%              usually given in units of seconds and the Doppler rate in
%              Hertz.
%
%This function implements the approximation that is given in Equation 3.16
%of [1]. The sign of m12 has been flipped as compared ot the paper, as the
%sign of it appears to be wrong in the paper. The approximations assumes
%that the time-bandwidth prodict (B*T) of each chirp is much greater than 1
%so that the spectrum of each chirp is essentially a rectangle. It also
%assumes that the noise is not correlated over time. Also, Equation 3.15
%for the chirp in [1] is the complex envelope of the chirp, not the actual
%chirp, so if one plots the real part, it looks like a downchirp followed
%by an upchirp. However, that is not the case if one shifts the frequency
%to passband (multiplies by exp(1j*2*pi*fc) to get it up to a high center
%frequency of fc).
%
%Given a discrete baseband sampled complex signal model
%y[n]=A*exp(j*omegaD*n)*s[n]+w[n], where exp(j*omegaD*n) provides the
%effects of the carrier frequency at discrete sampled time n, s[n] is the
%complex envelope of the signal, A is the complex amplitude, and w[n] is
%additive noise with variance sigma2, the SNR is related to the signal
%amplitude as
%A=sqrt(SNR*NSigma/epsilon)
%Where epsilon is the integral of the magnitude squared of s(t) (when
%considering the contuuous time signal) over the entire integration period.
%If s(t) is 1 over the pulselength and the pulselength is T seconds and
%there are NB pulses, then epsilon=T*NB. The value NSigma is sigma2*T0,
%where T0 is the sampling period. These relations are from [1].
%
%REFERENCES:
%[1] A. Dogandzic and A. Nehorai, "Cramér Rao bounds for estimating range,
%    velocity, and direction with an active array," IEEE Transactions on
%    Signal Processing, vol. 49, no. 6, pp. 1122-1137, Jun. 2001.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<6||isempty(returnFIM))
    returnFIM=false;
end

if(TR<T)
    error('The pulse repetition interval TR cannot be less than than the duration of the chirp in each pulse T.')
end

SNR1=SNRIn/P;

%The terms in the matrix in Equation 3.16.
m11=(1/12)*T^2*(1+(TR/T)^2*(P^2-1));
%Corrected sign as compared to the paper.
m12=(1/6)*pi*B*T;
m22=(1/3)*pi^2*B^2;

if(returnFIM)
    %The inverse of Equation 3.16, but with the delay component in the
    %upper left. This is the Fisher Information matrix, not the CRLB.
    CRLB=(2*P*SNR1)*[m22, m12;
                     m12, m11];
    D=[1,0;
       0,2*pi];
    CRLB=D*CRLB*D;
else
    %Equation 3.16 with the matrix inverted with the delay accuracy
    %component in the upper left.
    denom=(2*P*SNR1)*(m11*m22-m12^2);
    CRLB=[m11,-m12;
         -m12, m22]/denom;
    D=[1,0;
       0,1/(2*pi)];
    CRLB=D*CRLB*D;
end

if(nargin>6&&~isempty(c))
    if(returnFIM)
        D=[1/c, 0;
             0,-(fc)/c];
        CRLB=D*CRLB*D;
        if(useHalfRange)
            D=[2,0;
                0,2];
            CRLB=D*CRLB*D;
        end
    else
        D=[c, 0;
           0,-c/(fc)];
        CRLB=D*CRLB*D;

        if(useHalfRange)
            D=[1/2,0;
                0,1/2];
            CRLB=D*CRLB*D;
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
