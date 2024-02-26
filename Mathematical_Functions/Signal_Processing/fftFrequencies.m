function freqVals=fftFrequencies(N,fSamp,isFFTShifted)
%%FFTFREQUENCIES Return the frequency values that correspond to each entry
%       of the discrete Fourier trasnform (DFT) of a vector, so in fft(x)
%       or fftshift(fft(x)) if isFFTShifted is true.
%
%INPUTS: N The length of the vector whose DFT is taken.
%       fs The sample frequency, typically in Hertz. The default if omitted
%          or an empty matrix is passed is 1.
% isFFTShifted This is true if one want the frequencies for an fftshifted
%          vector. The default if omitted or an empty matrix is passed is
%          false.
%
%OUTPUTS: freqVals A 1XN vector where each entry is the frequency
%                  represented by the correspnding bin of the DFT of a
%                  length N vector.
%
%The relationship between the bins of a DFT and the corresponding frequency
%values are in many signal processing text books, such as in Chapter 15.1
%of [1], where here, we use the aliasing property of the DFT to determine
%the negative frequencies when one uses fftshift.
%
%EXAMPLE:
%Here, we create a 10Hz  sinusoid. The sinusoid is sampled above the
%Nyquist rate. Exactly the correct number of samples is taken so that is
%the sinusoid were to continue indefinitely, the next point after the end
%would be the first point here again --so what is FFt'd here is band
%limited to a single frequency. We see that the peak lines up exactly with
%the specified frequency and the frequency is reportede correctly by the
%function fftFrequencies.
% f=10;%10 Hertz.
% fN=2*f;%Nyquist rate.
% fSamp=4*fN;%Sample above the Nyquist rate.
% Ts=1/fSamp;%Sample period.
% wavePeriod=1/f;
% NSamp=10*wavePeriod/Ts;
% 
% t=(0:(NSamp-1))*Ts;
% x=sin(2*pi*f*t);
% figure(1)
% clf
% plot(t,x)
% 
% X=fft(x);
% f=fftFrequencies(NSamp,fSamp,false);
% figure(2)
% clf
% plot(f,abs(X))
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach,
%    3rd ed. Boston: McGraw Hill, 2006.
%
%September 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(isFFTShifted))
    isFFTShifted=false;
end

if(nargin<2||isempty(fSamp))
    fSamp=1;
end

freqVals=(0:1:(N-1))/N;

if(isFFTShifted)
    shiftAmount=floor(N/2);
    freqVals=circshift(freqVals,shiftAmount);
    sel=1:shiftAmount;%Negative frequency indices.
    freqVals(sel)=freqVals(sel)-1;
end

freqVals=freqVals*fSamp;

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
