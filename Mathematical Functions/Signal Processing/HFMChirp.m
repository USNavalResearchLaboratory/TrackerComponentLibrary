function [s,t]=HFMChirp(T,fStart,fEnd,tSamp,phi)
%%LFMCHIRP Generate a complex hyperbolic frequency modulated (HFM) signal.
%          It can be an upchirp or a downchirp. HFM chirps have a role in
%          radar and sonar signal processing. They are more Doppler
%          insensitive than linear FM chirps.
%
%INPUTS: T The duration of the chirp (seconds).
%   fStart The frequency at the start of the chirp (Hertz). This must be
%          >0.
%     fEnd The frequency at the end of the chirp (Hertz). This must be >0.
%    tSamp If this is scalar, it is the sampling rate of the chirp and the
%          number of samples depends on T. Otherwise, this is a numSampX1
%          array of times at which samples are desired. Samples before
%          time 0 and after time T are set to zero. If this parameter is
%          omitted or an empty matrix is passed, the Nyquist sampling
%          rate is used.
%      phi The phase offset of the chirp. If this parameter is omitted or
%          an empty matrix is passed, then the default of pi/2 is used,
%          which makes the real part of the signal start at zero.
%
%OUTPUT: s The numSampX1 complex chirp signal.
%        t The sample times used. If tSamp was not a scalar, then this is
%          just tSamp.
%
%HFM chirps are described in [1] along with a range bias model for target
%tracking that is relevant when Doppler filtering is not performed.
%
%EXAMPLE 1:
%Here is an example of an upchirp:
% T=1e-3;
% fStart=1e3;
% fEnd=50e3;
% tSamp=1/(100*2*fEnd);%100 times the Nyquist frequency.
% [s,t]=HFMChirp(T,fStart,fEnd,tSamp);
% figure(1)
% clf
% plot(t,real(s),'linewidth',2)
% h1=xlabel('t');
% h2=ylabel('s');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')

%EXAMPLE 2:
%Here is an example of a downchirp:
% T=1e-3;
% fStart=50e3;
% fEnd=1e3;
% tSamp=1/(100*2*fStart);%100 times the Nyquist frequency.
% [s,t]=HFMChirp(T,fStart,fEnd,tSamp);
% figure(2)
% clf
% plot(t,real(s),'linewidth',2)
% h1=xlabel('t');
% h2=ylabel('s');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] X. Song, P. Willett, and S. Zhou, "Range bias modeling for hyperbolic-
%    frequency-modulated waveforms in target tracking," IEEE Journal of
%    Oceanic Engineering, vol. 37, no. 4, pp. 670-679, Oct. 2012.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(tSamp))
   %Default to the Nyquist sampling rate.
   tSamp=1/(2*max(fStart,fEnd));
end

if(nargin<5||isempty(phi))
   phi=pi/2;%Phase offset so the signal starts at zero, by default.
end

%Equation 1
b=(fStart-fEnd)/(fStart*fEnd*T);

if(isscalar(tSamp))
    %If a scalar sampling rate is given.
    NPulse=fix(T/tSamp);
    t=(0:(NPulse-1))*tSamp;
else
    %If the samples themselves are given
    t=tSamp;
end

s=exp(1j*(phi+(2*pi/b)*log(1+b*fStart*t)));

sel=(t>T)&(t<0);
s(sel)=0;

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
