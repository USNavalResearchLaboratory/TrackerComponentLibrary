function [s,t,dsdt,d2sdt2]=HFMChirp(T,fStart,fEnd,tSamp,phi,makeBaseband)
%%HFMCHIRP Generate a complex hyperbolic frequency modulated (HFM) signal.
%          It can be an upchirp or a downchirp. HFM chirps have a role in
%          radar and sonar signal processing. They are more Doppler
%          insensitive than linear FM chirps. If fStart<fEnd, this is an
%          upchirp; if fStart>fEnd, this is a downchirp.
%
%INPUTS: T The duration of the chirp (seconds).
%   fStart The passband frequency at the start of the chirp (Hertz). This
%          must be >0. Use the makeBaseband input to shift to baseband if
%          desired.
%     fEnd The passband frequency at the end of the chirp (Hertz). This
%          must be >0. Use the makeBaseband input to shift to baseband if
%          desired.
%    tSamp If this is scalar value in a cell array, it is the sampling
%          rate of the chirp and the number of samples depends on T and is
%          given by round(T/tSamp{1}). Otherwise, this is a numSampX1 or
%          1XnumSamp array of times at which samples are desired. Samples
%          before time 0 and after time T are set to zero. If this
%          parameter is omitted or an empty matrix is passed, the Nyquist
%          sampling rate is used.
%      phi The phase offset of the chirp. If this parameter is omitted or
%          an empty matrix is passed, then the default of pi/2 is used,
%          which makes the real part of the signal start at zero.
% makeBaseband If this is true, a baseband signal is generated. Otherwise,
%          a passband signal is generated. The default if omitted or an
%          empty matrix is passed is false. The shift to baseband is
%          accomplished by multiplying by exp(-1j*2*pi*fc*t) where
%          fc=(fStart+fEnd)/2.
%
%OUTPUT: s The complex HFM chirp signal. This is either 1XnumSamp if tSamp
%          was a cell array or this has the dimensions of tSamp.
%        t The sample times used. If tSamp was not a cell array, then this
%          is just tSamp.
%     dsdt The first derivative of s with respect to time t. This has the
%          same dimensionality as s.
%   d2sdt2 The 1XnumSamp second derivative of s with respect to time t.
%          This has the same dimensionality as s.
%
%HFM chirps are described in [1] along with a range bias model for target
%tracking that is relevant when Doppler filtering is not performed.
%
%EXAMPLE 1:
%Here is an example of an upchirp:
% T=1e-3;
% fStart=1e3;
% fEnd=50e3;
% tSamp=1/(100*2*(fEnd-fStart));%100 times the Nyquist frequency.
% [s,t]=HFMChirp(T,fStart,fEnd,{tSamp});
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
%Here is an example of a downchirp:
% T=1e-3;
% fStart=50e3;
% fEnd=1e3;
% tSamp=1/(100*2*abs(fEnd-fStart));%100 times the Nyquist frequency.
% [s,t]=HFMChirp(T,fStart,fEnd,{tSamp});
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
%Here, we compare finite differencing to the explicit derivative. The
%results agree to more than 8 digits.
% fc=1e8;
% B=2e6;%2Mhz bandwidth.
% fStart=fc-B/2;
% fEnd=fc+B/2;
% T=2e-5;%Chirp duration in seconds.
% phi=pi/2;
% tPt=T/4.1;%Arbitrary time at which to evaluate the gradient.
% makeBaseband=true;
% [~,~,dsdt]=HFMChirp(T,fStart,fEnd,tPt,phi,makeBaseband);
% f=@(t)HFMChirp(T,fStart,fEnd,t,phi,makeBaseband);
% dsdtNumDiff=numDiff(tPt,f,1,2,1e-9);
% RelErr=abs((dsdt-dsdtNumDiff)./abs(dsdtNumDiff))
%
%EXAMPLE 4:
%This verifies that the second time derivative is consistent with numerical
%differentiation. The results agree to over 8 digits as implied by the
%relative error.
% fc=1e8;
% B=2e6;%2Mhz bandwidth.
% fStart=fc-B/2;
% fEnd=fc+B/2;
% T=2e-5;%Chirp duration in seconds.
% phi=pi/2;
% tPt=T/4.1;%Arbitrary time (with a nonzero signal).
% makeBaseband=true;
% [~,~,~,d2sdt2]=HFMChirp(T,fStart,fEnd,tPt,phi,makeBaseband);
% f=@(t)getOutputN(@(t)HFMChirp(T,fStart,fEnd,t,phi,makeBaseband),3,t);
% d2sdt2NumDiff=numDiff(tPt,f,1,2,1e-9);
% RelErr=abs((d2sdt2-d2sdt2NumDiff)./abs(d2sdt2NumDiff))
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
    tSamp={1/(2*(fEnd-fStart))};
end

if(nargin<5||isempty(phi))
    phi=pi/2;%Phase offset so the signal starts at zero, by default.
end

if(nargin<6||isempty(makeBaseband))
    makeBaseband=false;
end

%Equation 1
b=(fStart-fEnd)/(fStart*fEnd*T);

if(iscell(tSamp))
    %If a scalar sampling rate is given.
    NPulse=round(T/tSamp{1});
    t=(0:(NPulse-1))*tSamp{1};
else
    %If the samples themselves are given.
    t=tSamp;
end

if(makeBaseband)
    fc=(fEnd+fStart)/2;
    s=exp(1j*(-2*pi*fc*t+phi+(2*pi/b)*log(1+b*fStart*t)));
else
    s=exp(1j*(phi+(2*pi/b)*log(1+b*fStart*t)));
end

sel=(t>T)|(t<0);
s(sel)=0;
if(nargout>2)
    if(makeBaseband)
        dsdt=1j*(-2*fc*pi+(2*fStart*pi)./(1+b*fStart*t)).*s;
    else
        dsdt=((2*1j*fStart*pi)./(1+b*fStart*t)).*s;
    end
    dsdt(sel)=0;

    if(nargout>3)
        if(makeBaseband)
            d2sdt2=1j*(-2*fc*pi+(2*fStart*pi)./(1+b*fStart*t)).*dsdt-((2*1j*b*fStart^2*pi)./(1+b*fStart*t).^2).*s;
        else
            d2sdt2=((2*1j*fStart*pi)./(1+b*fStart*t)).*dsdt-((2*1j*b*fStart^2*pi)./(1+b*fStart*t).^2).*s;
        end
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
