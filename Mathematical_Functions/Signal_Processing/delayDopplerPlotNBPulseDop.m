function [delayDopPlot,Doppler,delay]=delayDopplerPlotNBPulseDop(x,y,M1,M2,wDoppler,T0)
%%DELAYDOPPLERPLOTNBPULSEDOP This function produces the complex
%               delay-Doppler plot for a narrowband signal assuming that
%               the same waveform is repeated every single pulse or that
%               each pulse ends with a region of non-broadcasting (>100%
%               duty cycle) and the maximum delay is less than that
%               duration. Range migration is not taken into account. Such
%               plots can be used for detection.
%
%INPUTS: x The NPulseX1 complex baseband waveform that is repeated every
%          pulse repetition interval (PRI). If NPulse<Ns, where Ns is the
%          number of samples in a PRI, then  it is assumed that the rest of
%          the signal is filled with zeros.
%        y The NsXNB set of complex baseband samples received over each of
%          the NB pulse repetition intervals. Ns is the number of samples
%          per interval. It is assumed that no gaps exist between
%          intervals.
%    M1,M2 These are the scaling factors to increase the resolution in
%          delay (M1) and Doppler (M2). These are integers. The number of
%          bins used in the plot increases by a factor of these. If omitted
%          or empty matrices are passed, a value of 1 is used.
% wDoppler This is an optional NBX1 or 1XNB window vector to lower the
%          Doppler sidelobes. One could use a window from, for example, the
%          windowFunSym function. If this parameter is omitted or an empty
%          matrix is passed, then the default of no windowing (the
%          rectangular window) is used. Range sidelobes should be lowered
%          either by windowing the pulse before broadcast or by other means
%          such as mismatched filtering.
%       T0 The sampling period. This parameter is only required if the
%          outputs Doppler and delay are desired.
%
%OUTPUTS: delayDopPlot This is an (M1*Ns)X(M2*NB) complex range-Doppler
%                      map. The values have been scaled so that peaks are
%                      approximate representations of the correct complex
%                      amplitude (not counting the effects of any
%                      windowing). However, as demonstrated in the example
%                      below, range migration will cause the complex
%                      amplitude of the peak to be biased from the true
%                      value.
%              Doppler An (M2*NB)X1 vector holding the values of the
%                      Doppler shifts used for each row on delayDopPlot.
%                      This requires that T0 be passed. These are positive
%                      and negative. These values are inverse seconds and
%                      when multiplied by the propagation speed (e.g. c=
%                      Constants.speedOfLight) divided by the carrier
%                      frequency fc will provide the range rate. If seconds
%                      are used for time and meters for distance, when
%                      multiplied by c/fc one gets range rate in meters per
%                      second. Positive is going away from the sensor.
%                delay An (M1*Ns)X1 vector holding the values of the delays
%                      used for each columns on delayDopPlot. This requires
%                      that T0 be passed. These are all non-negative.
%                      Multiplied by the propagation speed c this gives the
%                      range.
%
%Range-Doppler plots are discussed in many sources, but many aspects of the
%approximations are often omitted, so a brief derivation is given here.
%
%Let x(t) be the complex baseband waveform sent, concatenating all pulses
%together. The received passband signal is modeled as
%x_p(t)=real(exp(1j*2*pi*fc)*x(t))
%where x(t) is the complex baseband signal. To come to baseband, one
%would take the real part of exp(1j*2*pi*fc) times the complex baseband
%signal. In practice, one would use quadrature demodulation to take the
%real measured signal and obtain the in-phase (real) and quadrature
%(imaginary) parts of the complex baseband signal.
%
%The delayed and Doppler-shifted passband signal is
%y_p(t)=real(A*exp(1j*2*pi*fc*(t-tau-a*t)*x(t-tau-a*t))
%where tau is the delay measured from the time the pulse begins
%broadcasting, a is a Doppler shift, and A is a complex amplitude. This
%assumes a target with a constant range rate (constant Doppler frequency).
%For targets accelerating in range, additional terms are needed. The
%complex baseband received signal is thus
%y(t)=A*exp(-1j*2*pi*fc*(tau+a*t))*x(t-tau-a*t)
%
%The range-Doppler map is just a matched filter. A matched filter
%integrates the product of the complex conjugate of the ideal signal times
%the received signal. Here, the matched filter while observing over a time
%T is:
%integral_0^T exp(1j*2*pi*fc*(tau+a*t))*conj(x(t-tau-a*t))*y(t) dt
%In the problem at hand, the signal is approximated as a set of NB blocks
%or pulses of duration TB. Thus, T=NB*T. To simplify the signal
%processing, we assume either that all blocks have the same waveform or
%that each block ends with a period of no broadcasting (<100% duty cycle)
%and the maximum delay does not exceed that period. Additionally, the
%approximation is made that the delays due to the Doppler shift are
%constant over each block. All together, this makes the use of circular
%convolutions (and hence FFTs in the signal processing) possible.
%The matched filter under such a model is
%integral_0^T sum_{i=0}^{N_B-1} exp(1j*2*pi*fc*(tau+a*i*TB))*conj(x_i(t-tau-a*i*TB))*y_i(t) dt
%where x_i and y_i are the signals from the ith block. It is assumed that
%outside of a time windows of length T_B each of those signals is zero.
%To further simplify the signal processing, we assume that the Doppler
%shift delay in the x term is not significant (related to a narrowband
%approximation). Thus, the signal processing model for the matched filter
%is 
%exp(1j*2*pi*fc*tau)*sum_{i=0}^{N_B-1}exp(1j*2*pi*fc*a*i*TB)*integral_0^T conj(x_i(t-tau))*y_i(t) dt
%As the signal y_i is sampled, we will approximate the integral in time
%using a Riemann sum. We will also say that tauHat=tau*T0 (T0 is the sample
%period) where tauHat is an integer >=0 for the discrete delay. Also, we
%will throw out the phase constant of exp(1j*2*pi*fc*tau) in front. The
%result is thus:
%sum_{i=0}^{N_B-1}exp(1j*2*pi*fc*a*i*TB)*sum_{k=0}^{Ns-1} conj(x(k-tauHat))*y(i,k) 
%Thus, the integral has turned into a discrete sum over the samples in each
%block (the arguments of x and y access the discrete elements). There is no
%block marking on x, because it is assumed to be the same for all blocks.
%Looking at the form of the sum, it can be shown (for example, see Chapter
%5.7 of [1] for useful identities) that it is equivalent to the circular
%convolution of flipud(conj(x)) with y. This can be evaluated for all
%discrete delays tau and block i in parallel using ffts and an ifft. This
%is what is implemented below via the circConv function. The outer sum in
%the above can be seen to be NB times the ifft taken across all values of i
%when we replace a with a=aHat/(fc*TB*NB) where aHat is a discrete Doppler
%shift from 0 to (NB-1). Thus, the ifft function, as used below, evaluated
%all of the Doppler shifts in parallel for the matched filter, scaled by
%1/NB.
%
%Thus, as described above, we are evaluating 
%(1/NB)*sum_{i=0}^{N_B-1}exp(1j*2*pi*aHat*i/NB)*sum_{k=0}^{Ns-1} conj(x(k-tauHat,i))*y(i,k) 
%for discrete values of aHat and tauHat. However, it would be nice if the
%output could be related to the amplitude A. Well, if we assume that tauHat
%and aHat perfectly match the original signal AND the original signal was
%generated using the various narrowband simplications mentioned, then we
%have to multiply the result by 1/(x'*x). This is what is done below.
%However, opne must realize that the various approximations means that the
%amplitude estimates can still be biased. This is demonstrated in the
%example below.
%
%Often, the sidelobes of the signal in Doppler might be high. Thus, the
%input wDoppler can be used to apply tapering across each bin before the
%final ifft of the algorithm.
%
%The final thing to note is that the above method describes only processing
%for positive Doppler values. Due to the aliasing of large values to
%negative values with circular convolutions and FFTs, we use fftshift to
%get positive and negative values in the unambiguous Doppler region.
%
%Note that TB*c gives the maximum unambiguous range, where c is the speed
%of propagation in the medium, for example c=Constants.speedOfLight.
%Note that c/(2*fc*TB) is the magnitude of the maximum unambiguous range
%rate. One can get those as
% range=c*delay;
% rangeRate=Doppler*(c/fc);
%
%EXAMPLE:
%Our signal is a up-chirp. Here, we create the time-delayed Doppler shifted
%signal at baseband. The target signal is generated in two ways: Once
%including all of the phase shifts associated with range migration, and
%once with the ideal signal processing model that is used to derive the
%range-Doppler map using Fourier transforms to make it fast. We choose the
%target location to be EXACTLY on a delay-Doppler bin. It will be seen
%that the complex amplitude obtained using the ideal model matches the
%complex amplitude used in generating the signal (no noise is added).
%However, the complex amplitude obtained when using the more realistic
%model for the signal, including all range migration, is biased. This
%highlights how the lack of migration compensation biases complex amplitude
%determination. The range-Doppler plot of the realistic signal is
%displayed.
% fc=1e9;%1GHz carrier frequency.
% B=2e6;%2Mhz bandwidth.
% %Baseband start and end frequencies.
% fStart=-B/2;
% fEnd=B/2;
% %Sampling rate is two times the Nyquist rate.
% T0=1/(2*2*fEnd);%Sampling period in seconds.
% T=2e-5;%Chirp duration in seconds.
% 
% PRF=2000;%Pulse repetition frequency (Hertz)
% TB=1/PRF;%The pulse repetition period.
% %The number of samples per PRI. The above parameters were chosen so that
% %this is an integer. Fix just deals with finite precision errors.
% Ns=fix(TB/T0);
% 
% %Generate the reference signal. This is an up-chirp.
% x=LFMChirp(T,fStart,fEnd,T0);
% x=x(:);
% 
% %We will use 64 pulse repetition intervals.
% NB=64;
% 
% %True target parameters
% c=Constants.speedOfLight;
% rTrue=512*c*T0;%Target placed on
% tau=rTrue/c;%The true delay (s).
% 
% %The true range rate (m/s).
% rrTrue=6/NB*(c/(fc*TB));
% a=rrTrue/c;
% 
% %Allocate space for the received signal. The first dimensions is "fast
% %time"; the second dimension if "slow time".
% yReal=zeros(Ns,NB);
% yIdeal=zeros(Ns,NB);
% 
% %Create the received signal, properly delayed and Doppler shifted for
% %each PRI. The same waveform is used in each PRI.
% t=0:T0:((Ns-1)*T0);%Sample times
% for i=0:(NB-1)
%     %The subtraction of i*TB deals with the start time of this pulse if
%     %it is not time zero.
%     tCur=t-a*t-tau-i*TB;
%     yReal(:,i+1)=512*exp(-1j*2*pi*fc*(tau+a*t)).*LFMChirp(T,fStart,fEnd,tCur);
%     
%     tCur=t-tau-i*TB;
%     yIdeal(:,i+1)=512*exp(-1j*2*pi*fc*(tau+a*i*TB)).*LFMChirp(T,fStart,fEnd,tCur);
%     
%     t=t+TB;%Increment to the next time step.
% end
% 
% M1=1;
% M2=1;
% %If desired, one could use a window in the Doppler domain. For example,
% %by using wDoppler=windowFunSym(NB,'Blackman',1); Here, we choose to use
% %no window (same as a rectangular window).
% wDoppler=[];
% valsReal=delayDopplerPlotNBPulseDop(x,yReal,M1,M2,wDoppler,T0);
% [valsIdeal,Doppler,delay]=delayDopplerPlotNBPulseDop(x,yIdeal,M1,M2,wDoppler,T0);
% 
% %Note that the amplitude of the ideal signal matches the true amplitude
% %used, but that of the real signal (with range migration) is biased.
% [magReal,idxReal]=max(abs(valsReal(:)));
% [magIdeal,idxIdeal]=max(abs(valsIdeal(:)));
% amplitudeMagReal=magReal
% amplitudeMagIdeal=magIdeal
% 
% range=c*delay;
% rangeRate=Doppler*(c/fc);
% 
% %We will make sure that the maximum point on the plots matches the true
% %inputs, since the true inputs were made to align with range-Doppler cells.
% [idxR,idxRR]=ind2sub(size(valsReal),idxReal);
% rDiffReal=range(idxR)-rTrue
% rrDiffReal=rangeRate(idxRR)-rrTrue
% 
% [idxR,idxRR]=ind2sub(size(valsIdeal),idxIdeal);
% rDiffIdeal=range(idxR)-rTrue
% rrDiffIdeal=rangeRate(idxRR)-rrTrue
% 
% figure(1)
% clf
% imagesc([rangeRate(1),rangeRate(end)],[range(1), range(end)]/1e3,10*log10(abs(valsReal)));
% set(gca,'YDir','normal')
% caxis([-30 27])
% colormap(jet(256))
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach,
%    3rd ed. Boston: McGraw Hill, 2006.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(M1))
    M1=1;
end

if(nargin<4||isempty(M2))
   M2=1; 
end

NB=size(y,2);
Ns=size(y,1);

%The extra circshift aligns the zero-Doppler point to the position assumed
%in the delay output.
rangePlots=circshift(circConv(flipud([conj(x);zeros(Ns-length(x),1)]),y,Ns),1);

if(nargin>4&&~isempty(wDoppler))
    %Perform windowing to lower Doppler sidelobes.
    rangePlots=bsxfun(@times,wDoppler(:).',rangePlots);
end

delayDopPlot=ifftshift(ifft(rangePlots,NB*M2,2),2);

%Remove the scaling at the matched points so that we can get the true
%complex amplitude back.
delayDopPlot=M1*M2*delayDopPlot/(x'*x);

if(nargout>1)
    %T0 must be provided for these outputs.
    
    delay=(0:(T0/M1):((T0/M1)*(M1*Ns-1))).';
    %The pulse repetition interval.
    TB=Ns*T0;

    numNeg=(NB*M2-1)-ceil((NB*M2-1)/2)+1;
    Doppler=([(-numNeg):1:-1,0:(floor((NB*M2-1)/2)-1)].'/(NB*M2))*(1/TB);
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
