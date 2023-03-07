function CRLB=delayDopplerRateCRLBLFMApprox(B,T,P,TR,SNRIn,c,fc)
%DELAYDOPPLERRATECRLBLFMAPPROX Determine an approximate Cramér-Rao lower
%               bound matrix for estimating the delay and Doppler rate when
%               receiving a chain of identical linearly frequency modulated
%               (LFM) chirps. 
%
%INPUTS: B The bandwidth of each chirp, usually in Hertz. This is defined
%          to be the absolute value of the difference of the starting
%          instantaneous frequency and the ending instantaneous frequency
%          of the chirps, as used in the function LFMChirp.
%        T The duration of the chirp in each pulse, usually in seconds.
%        P The integer number of pulses in the train; P>1.
%       TR The duration of each pulse in the pulse train, usually in
%          seconds.
%    SNRIn The unitless ratio of the instantaneous signal power to average
%          noise power ratio of each sample. For a signal
%          x(t)=A*exp(1i*f(t))+w(t) where w(t) is noise, this is
%          A^2/E(w(t)^2). It is assumed that the noise is not correlated
%          over time.
%     c,fc These two optional parameters are only provided if one wishes to
%          have the results given in terms of range and range rate. c is
%          the propagation speed of the medium (e.g.
%          Constants.speedOfLight) and fc is the center frequency of the
%          transmission.
%
%OUTPUTS: CRLB The 2X2 Cramér-Rao lower bound for a vector [delay;Doppler
%              rate]. The delay is usually given in units of seconds and
%              the Doppler rate in radians per second. The example below
%              shows how to turn these quantities into range and velocity.
%              However, if one provides c and fc, to this function, then
%              the function will provide the outputs as range and range
%              rate.
%
%This function implements the approximation that is given in Equation 3.16
%of [1]. The approximations assumes that the time-bandwidth prodict (B*T)
%of each chirp is much greater than 1 so that the spectrum of each chirp is
%essentially a rectangle. It also assumes that the noise is not correlated
%over time.
%
%EXAMPLE:
% B=10e6;%10Mhz bandwidth.
% T=250e-6;%Chirp duration in seconds.
% P=16;%Number of pulses used.
% TR=1e-3;%The pulse repetition interval.
% %-10dB per pulse, so this is the equivalent SNRIn
% SNRIn=10^(-10/10)/T;
% CRLB=delayDopplerRateCRLBLFMApprox(B,T,P,TR,SNRIn);
% 
% %Rather than having a standard deviation for delay, we would rather use
% %range. Delay times the speed of light equals round-trip range. Thus, we
% %multiple the delay components in the CRLB by the speed of light (for the
% %cross term) or the square of the speed of light
% c=Constants.speedOfLight;
% %Cross terms
% CRLB(1,2)=CRLB(1,2)*c;
% CRLB(2,1)=CRLB(1,2);
% CRLB(1,1)=CRLB(1,1)*c^2;
% %Next, we want to convert the Doppler rate into a range rate with units
% %of meters per second. This requires the carrier frequency fc and the
% %propagation speed c. We will use 2GHz as the carrier frequency.
% fc=2e9;
% scal=-c/(2*pi*fc);
% CRLB(2,2)=CRLB(2,2)*scal^2;
% CRLB(1,2)=CRLB(1,2)*scal;
% CRLB(2,1)=CRLB(1,2);
% CRLB
% %One will see that the CRLB is around what one might expect. Ignoring the
% %cross terms, one can take the square roots of the entries to get a
% %standard deviation in range of about 9.2408 meters and in range rate of
% %about 2.8931 meters per second.
%
%REFERENCES:
%[1] A. Dogandzic and A. Nehorai, "Cramér Rao bounds for estimating range,
%    velocity, and direction with an active array," IEEE Transactions on
%    Signal Processing, vol. 49, no. 6, pp. 1122-1137, Jun. 2001.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 3.17.
SNR1=SNRIn*T;

%Equation 3.16 with the matrix inverted with the range accuracy in the
%upper left.
CRLB=1/(2*P*SNR1)*[3*(T^2+(P^2-1)*TR^2)/(B^2*(P^2-1)*pi^2*TR^2), 6*T/(B*(P^2-1)*pi*TR^2);
                   6*T/(B*(P^2-1)*pi*TR^2),                      12/(TR^2*(P^2-1))];
               
if(nargin>5&&~isempty(c))
    CRLB(1,2)=CRLB(1,2)*c;
    CRLB(2,1)=CRLB(1,2);
    CRLB(1,1)=CRLB(1,1)*c^2;
    
    scal=-c/(2*pi*fc);
    CRLB(2,2)=CRLB(2,2)*scal^2;
    CRLB(1,2)=CRLB(1,2)*scal;
    CRLB(2,1)=CRLB(1,2);
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
