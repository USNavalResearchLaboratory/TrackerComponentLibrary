function val=HFMChirpEnergy(T,fStart,fEnd,phi)
%%HFMCHIRPENERGY Given the parameters for a hyperbolic frequency modulated
%                (HFM) chirp signal, which can be an upchirp or a
%                downchirp, determine the scaled energy of the continuous-
%                time signal. For a signal x(t) going from time 0 to T, the
%                energy is the integral_0^T of x(t)^2 dt divided by the
%                magnitude of the load driven. Here, we omit the dividing
%                term, so the result is the scaled energy. For a signal in
%                Volts and the time of integration in sections, the output
%                is in Volt^2-seconds. Dividing that by a resistance in
%                Ohms leads to an energy in Joules.
%
%INPUTS:  T The duration of the chirp (seconds).
%    fStart The frequency at the start of the chirp (Hertz). This must be
%           >0.
%      fEnd The frequency at the end of the chirp (Hertz). This must be >0.
%       phi The phase offset of the chirp. If this parameter is omitted or
%           an empty matrix is passed, then the default of pi/2 is used,
%           which makes the real part of the signal start at zero.
%
%OUTPUT: val The scaled energy of the chirp.
%
%The time variable for integration is taken to be from 0 to T. The real
%passband chirp signal is assumed to have the form
% cos(phi+(2*pi*log(1+b*fStart*t])/b)
%where b=(fStart-fEnd)/(fStart*fEnd*T). This is consistent with the
%function HFMChirp. The solution was obtained by performing symbolic
%integration.
%
%EXAMPLE:
%Here, we find the energy in a 1ms chirp going from 50kHz down to 1kHz with
%a phase offset of pi/2.
% T=1e-3;
% fStart=50e3;
% fEnd=1e3;
% phi=pi/2;
% val=HFMChirpEnergy(T,fStart,fEnd,phi)
%One will get a result of about 5.0103e-04 Volt^2*seconds assuming the
%signal is given in Volts.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(phi))
   phi=pi/2;%Phase offset so the signal starts at zero, by default.
end

%Equation 1
b=(fStart-fEnd)/(fStart*fEnd*T);

temp=2*phi+4*pi*log(1+b*fStart*T)/b;
sinVal1=sin(temp);
cosVal1=cos(temp);

temp=2*phi;
sinVal2=sin(temp);
cosVal2=cos(temp);

temp=fStart*(b^2+16*pi^2);

val=(1/(2*temp))*(temp*T-b*cosVal2+b*(1+b*fStart*T)*cosVal1-4*pi*sinVal2+4*pi*(1+b*fStart*T)*sinVal1);

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
