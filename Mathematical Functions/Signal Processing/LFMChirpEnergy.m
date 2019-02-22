function val=LFMChirpEnergy(T,fStart,fEnd,phi)
%%LFMCHIRPENERGY Given the parameters for a linearly frequency modulated
%                (LFM) chirp signal, which can be an upchirp or a
%                downchirp, determine the scaled energy of the continuous-
%                time signal. For a signal x(t) going from time 0 to T, the
%                energy is the integral_0^T of x(t)^2 dt divided by the
%                magnitude of the load driven. Here, we omit the dividing
%                term, so the result is the scaled energy. For a signal in
%                Volts and the time of integration in sections, the output
%                is in Volt^2-seconds. Dividing that by resistance in Ohms
%                leads to an energy in Joules.
%
%INPUTS:  T The duration of the chirp (seconds).
%    fStart The positive frequency at the start of the chirp (Hertz).
%      fEnd The positive frequency at the end of the chirp (Hertz).
%       phi The phase offset of the chirp. If this parameter is omitted or
%           an empty matrix is passed, then the default of pi/2 is used,
%           which makes the real part of the signal start at zero.
%
%OUTPUT: val The scaled energy of the chirp.
%
%The time variable for integration is taken to be from 0 to T. The real
%passband chirp signal is assumed to have the form
%cos(2*pi*fStart*t+pi*alpha*t.^2+phi)
%where alpha=(fEnd-fStart)/T. This is consistent with the function
%LFMChirp. The solution was obtained by performing symbolic integration and
%simplifying the results to be in terms of Fresnel integrals. A special
%case is used for when alpha is negative so as to avoid taking the square
%roots of negative numbers.
%
%EXAMPLE:
%Here, we find the energy in a 1ms chirp going from 1MHz to 10MHz with a
%phase offset of pi/2.
% T=1e-3;
% fStart=1e6;
% fEnd=10e6;
% phi=pi/2;
% val=LFMChirpEnergy(T,fStart,fEnd,phi)
%One will get a result of about 0.005 Volt^2*seconds assuming the signal is
%given in Volts.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(phi))
    phi=pi/2;
end

alpha=(fEnd-fStart)/T;

if(alpha>0)
    sqrtAlpha=sqrt(alpha);

    [S1,C1]=FresnelInt(2*fStart/sqrtAlpha);
    [S2,C2]=FresnelInt(2*(fStart+T*alpha)/sqrtAlpha);
    temp=2*(fStart^2*pi/alpha-phi);
    sinVal=sin(temp);
    cosVal=cos(temp);

    val=1/(4*sqrtAlpha)*(2*T*sqrtAlpha+cosVal*(C2-C1)+sinVal*(S2-S1));
else
    alpha=abs(alpha);
    sqrtAlpha=sqrt(alpha);
    
    [S1,C1]=FresnelInt(2*fStart/sqrtAlpha);
    [S2,C2]=FresnelInt(2*(-fStart+T*alpha)/sqrtAlpha);
    temp=2*(fStart^2*pi/alpha+phi);
    sinVal=sin(temp);
    cosVal=cos(temp);
    
    val=1/(4*sqrtAlpha)*(2*T*sqrtAlpha+cosVal*(C1+C2)+sinVal*(S1+S2));
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
