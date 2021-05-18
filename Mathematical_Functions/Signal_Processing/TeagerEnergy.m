function E=TeagerEnergy(x,dim)
%%TEAGERENERGY Compute the Teager energy of a discrete signal as defined in
%              [1]. The Teager energy is a type of instantaneous "energy"
%              signal that produces a result related to the rate of change
%              of the amplitude and frequency of the original signal. It
%              has found use in a number of applications, primarily related
%              to acoustic processing, such as in the wavelet-based speech
%              enhancement of [2].
%
%INPUTS: x A vector signal or a matrix of signals whose Teager energy is
%          desired.
%      dim The dimension of x over which the Teager energy is to be found.
%          For example, for a NXnumSigs set of signals, one would use
%          dims=1 as each column is a separate signal.
%
%OUTPUTS: E The Teager energy. This has the same size as x.
%
%EXAMPLE 1:
%Here we compute the Teager energy signal of an up-chirp followeb by a down
%chirp.
% fStart=10;
% fEnd=100;
% T=1;
% 
% %We want the chirp to be sampled at 100 times the Nyquist frequency.
% T0=1/(100*2*fEnd);%The sampling period.
% %We want to sample just the duration of the signal. Thus, time goes from
% %T=0 to T=1. This means that we need N=T/T0 samples.
% N=T/T0;
% 
% %Make a chirp up and then a chirp down
% [x1,t]=LFMChirp(T,fStart,fEnd,T0);
% x2=LFMChirp(T,fEnd,fStart,T0);
% x=real([x1(:);x2(:)]);
% t=[t(:);t(end)+t(2)-t(1)+t(:)];
% 
% E=TeagerEnergy(x);
% 
% figure(1)
% clf
% subplot(2,1,1)
% hold on
% plot(t,x,'-r','linewidth',2)
% title('Original Signal')
% subplot(2,1,2)
% plot(t,E,'-r','linewidth',2)
% title('Transformed Signal')
%
%EXAMPLE 2:
%Here we compute the Teager energy of a sum of sinusoids. It essentially
%just follows the beat pattern.
% n=(0:100).';
% x=sin(pi/6*n)+sin(pi/4*n);
% 
% E=TeagerEnergy(x);
% 
% figure(2)
% clf
% subplot(2,1,1)
% hold on
% plot(n,x,'-r','linewidth',2)
% title('Original Signal')
% subplot(2,1,2)
% plot(n,E,'-r','linewidth',2)
% title('Transformed Signal')
%
%REFERENCES:
%[1] J. F. Kaiser, "On a simple algorithm to calculate the 'energy' of a
%    signal," in Proceedings of the International Conference on Acoustics,
%    Speech and Signal Processing, Albuquerque, NM, 3-6 Apr. 1990, pp.
%    381-384.
%[2] M. Bahoura and J. Rouat, "Wavelet speech enhancement based on the
%    Teager energy operator," IEEE Signal Processing Letters, vol. 8, no.
%    1, pp. 10-12, Jan. 2001.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(dim))
    E=x.^2-circshift(x,-1).*circshift(x,1);
else
    E=x.^2-circshift(x,-1,dim).*circshift(x,1,dim);    
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
