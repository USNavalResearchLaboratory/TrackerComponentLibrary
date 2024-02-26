function [s,t,dsdt,d2sdt2]=LFMChirp(T,fStart,fEnd,tSamp,phi)
%%LFMCHIRP Generate a complex linearly frequency modulated (LFM) signal. It
%          can be an upchirp or a downchirp. LFM chirps have a role in
%          radar and sonar signal processing. Letting fStart=-fEnd will
%          generate a chirp at baseband with bandwidth 2*fEnd.
%
%INPUTS:  T The duration of the chirp (seconds).
%    fStart The frequency at the start of the chirp (Hertz). For a signal
%           at baseband, this should be -B/2, where B is the bandwidth.
%      fEnd The frequency at the end of the chirp (Hertz). For a signal
%           at baseband, this should be B/2, where B is the bandwidth.
%     tSamp If this is scalar value in a cell array, it is the sampling
%           rate of the chirp and the number of samples depends on T and is
%           given by round(T/tSamp{1}). Otherwise, this is a numSampX1 or
%           1XnumSamp array of times at which samples are desired. Samples
%           before time 0 and after time T are set to zero. If this
%           parameter is omitted or an empty matrix is passed, the Nyquist
%           sampling rate is used.
%       phi The phase offset of the chirp. If this parameter is omitted or
%           an empty matrix is passed, then the default of pi/2 is used,
%           which makes the real part of the signal start at zero.
%
%OUTPUT: s The complex LFM chirp signal. This is either 1XnumSamp if tSamp
%          was a cell array or this has the dimensions of tSamp.
%        t The sample times used. If tSamp was not a cell array, then this
%          is just tSamp.
%     dsdt The first derivative of s with respect to time t. This has the
%          same dimensionality as s.
%   d2sdt2 The 1XnumSamp second derivative of s with respect to time t.
%          This has the same dimensionality as s.
%
%LFM chirps are common in radar. They are discussed in Section 8.2 of [1].
%In [2], biases relevant to tracking from LFM signals processed without
%Doppler filtering are discussed. The real component of the complex
%waveform is cos(2*pi*fStart*t+pi*alpha*t.^2+phi), where
%alpha=(fEnd-fStart)/T.
%
%EXAMPLE 1:
%Here is an example of an upchirp:
% T=1e-3;
% fStart=0;
% fEnd=10e3;
% tSamp=1/(100*2*fEnd);%100 times the Nyquist frequency.
% [s,t]=LFMChirp(T,fStart,fEnd,{tSamp});
% figure(1)
% clf
% plot(t,real(s),'linewidth',2)
% h1=xlabel('t');
% h2=ylabel('s');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here is an example of an downchirp:
% T=1e-3;
% fStart=10e3;
% fEnd=0;
% tSamp=1/(100*2*fStart);%100 times the Nyquist frequency.
% [s,t]=LFMChirp(T,fStart,fEnd,{tSamp});
% figure(2)
% clf
% plot(t,real(s),'linewidth',2)
% h1=xlabel('t');
% h2=ylabel('s');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 3:
%This verifies that the time derivative is consistent with numerical
%differentiation. The results agree to over 10 digits as implied by the
%relative error.
% B=2e6;%2Mhz bandwidth.
% T=2e-5;%Chirp duration in seconds.
% phi=pi/2;
% tPt=T/4.1;%Arbitrary time at which to evaluate the gradient.
% [~,~,dsdt]=LFMChirp(T,-B/2,B/2,tPt,phi);
% f=@(t)LFMChirp(T,-B/2,B/2,t,phi);
% dsdtNumDiff=numDiff(tPt,f,1,2,1e-10);
% RelErr=abs((dsdt-dsdtNumDiff)./abs(dsdtNumDiff))
%
%EXAMPLE 4:
%This verifies that the second time derivative is consistent with numerical
%differentiation. The results agree to over 10 digits as implied by the
%relative error.
% B=2e6;%2Mhz bandwidth.
% T=2e-5;%Chirp duration in seconds.
% phi=pi/2;
% tPt=T/4.1;%Arbitrary time (with a nonzero signal).
% [~,~,~,d2sdt2]=LFMChirp(T,-B/2,B/2,tPt,phi);
% f=@(t)getOutputN(@(t)LFMChirp(T,-B/2,B/2,t,phi),3,t);
% d2sdt2NumDiff=numDiff(tPt,f,1,2,1e-10);
% RelErr=abs((d2sdt2-d2sdt2NumDiff)./abs(d2sdt2NumDiff))
%
%REFERENCES:
%[1] M. I. Skolnik, Introduction to Radar Systems, 3rd ed. Boston:
%    McGraw Hill, 2001.
%[2] R. J. Fitzgerald, "Effects of range-Doppler coupling on chirp radar
%    tracking accuracy," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 10, no. 4, pp. 528-532, Jul. 1974.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(tSamp))
   %Default to the Nyquist sampling rate.
   tSamp={1/(2*(fEnd-fStart))};
end

if(nargin<5||isempty(phi))
    phi=pi/2;%Phase offset so the signal starts at zero, by default.
end

alpha=(fEnd-fStart)/T;

if(iscell(tSamp))
    %If a scalar sampling rate is given.
    NPulse=round(T/tSamp{1});
    t=(0:(NPulse-1))*tSamp{1};
else
    %If the samples themselves are given
    t=tSamp;
end
%The real signal corresponding to the complex waveform given below is
%cos(2*pi*fStart*t+pi*alpha*t.^2+phi);
s=exp(1j*(2*pi*fStart*t+pi*alpha*t.^2+phi));

sel=(t>T)|(t<0);
s(sel)=0;
if(nargout>2)
    dsdt=2*1j*pi*s.*(fStart+alpha*t);
    dsdt(sel)=0;

    if(nargout>3)
        d2sdt2=-2*pi*s.*(-1j*alpha+2*pi*(fStart+alpha*t).^2);
        d2sdt2(sel)=0;
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
