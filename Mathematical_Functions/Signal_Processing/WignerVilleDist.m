function [val,t,f]=WignerVilleDist(s,a,freqSamps)
%%WIGNERVILLEDIST Compute the discrete Wigner-Ville distribution. This is a
%    time-frequency distribution. Multiple discretizations of the
%    continuous-time version of the Wigner-Ville distribution exist. This
%    function implements the discrete Wigner-Ville distribution C (DWVD-C),
%    which is described in [1] with a detailed analysis given in [2]. As
%    noted in [1], this specific definition has a number of nice
%    properties. The Wigner-Ville distribution is related to the wideband
%    ambiguity function. An approximate algorithm is also implemented,
%    which allows one to obtain the transform for specific frequencies
%    while decimating the time-scale.
%
%INPUTS: s An NX1 or 1XN signal. s can be real (in which case the complex
%          analytic signal will be found) or s can be complex (representing
%          an analytic signal).
%        a If a is provided, then decimation in time is performed. The time
%          decimation factor a is an integer >=1. It is required that
%          (2*N/a) be an integer. If this parameter is omitted or an empty
%          matrix is provided, then a=1 is used.
% freqSamps If only specific frequency values are desired (rather than all
%          frequencies), then this input should be provided and is a set of
%          indices from 1 to N of which frequencies in the full transform
%          are desired. If omitted or an empty matrix is passed, then
%          freqSamps=1:N is used.
%
%OUTPUTS: val The ((2*N)/a)XnumFreqSamps discrete Wigner-Ville
%             distribution. The first index is time; the second index is
%             frequency.
%           t The values of the time-samples of the time samples
%             corresponding to the rows of the distribution. These are in
%             units of the sample period.
%           f The normalized values of the frequency samples corresponding
%             to the columns of val. Multiply these by 1/(2*T0) to get the
%             actual frequency values, where T0 is the sampling period.
%
%The non-decimated function is based on Algorithm 4 in Section 5.3.2 of
%[2]. The "kernel" used is just 1. The decimated algorithm is algorithm 6
%of [2] and is only used if more than one inputs are provided to this
%function. The function real2ComplexSignalW is used to obtain the analytic
%signal when s is real, because the algorithm implemented there was shown
%in [2] and [3] to have less aliasing issues then Hough-transform based
%methods.
%
%EXAMPLE 1:
%Here, we consider an HFM chirp.
% fStart=3e3;
% fEnd=10e3;%10kHz
% 
% %Sampling rate is two times the Nyquist rate for the highest frequency.
% T0=1/(2*2*fEnd);%Sampling period in seconds.
% T=10e-3;%Chirp duration in seconds.
% 
% %Generate the real signal, here a chirp.
% xSig=real(HFMChirp(T,fStart,fEnd,T0));
% 
% [G,t,f]=WignerVilleDist(xSig);
% t=t*T0;
% f=f./(2*T0);
% 
% figure(1)
% clf
% imagesc([t(1),t(end)]*1e3,[f(1),f(end)]/1e3,20*log10(abs(G.'/max(G(:)))));
% set(gca,'YDir','normal')
% colormap(jet(256))
% colorbar()
% caxis([-50, 0])
% title('Wigner-Ville Distribution (dB from Peak)')
% h1=xlabel('Time (ms)');
% h2=ylabel('Frequency kHz');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here we look at obtaining the energy of an LFM chirp.
% fStart=3e3;
% fEnd=10e3;%10kHz
% 
% %Sampling rate is two times the Nyquist rate for the highest frequency.
% T0=1/(2*2*fEnd);%Sampling period in seconds.
% T=10e-3;%Chirp duration in seconds.
% 
% %Generate the analytic signal, here a chirp.
% xSig=LFMChirp(T,fStart,fEnd,T0);
% 
% G=WignerVilleDist(xSig);
% signalEnergy=2*T0*sum(sum(abs(G).^2))
% analyticEnergy=LFMChirpEnergy(T,fStart,fEnd)
%One will see that the signal energy from the Wigner-Ville transform and
%the analytic solution are aboth about 0.0050 Volt^2*seconds. That is, they
%agree. If one were to pass the real signal, the function used in
%WignerVilleDist to get the analytic signal reduces aliasing effects, but
%it also makes the power computation slightly less accurate.
%
%REFERENCES:
%[1] J. M. O'Toole, M. Mesbah, and B. Boashash, "Improved discrete
%    definition of quadratic time-frequency distributions," IEEE
%    Transactions on Signal Processing, vol. 58, no. 2, pp. 906-911, Feb.
%    2010.
%[2] J. M. O'Toole, "Discrete quadratic time-frequency distributions:
%    Definition, computation, and a newborn electroencephalogram
%    application," Ph.D. dissertation, School of Medicine, The University
%    of Queensland, Australia, Mar. 2009.
%[3] J. M. O'Toole, M. Mesbah, and B. Boashash, "A new discrete analytic
%    signal for reducing aliasing in the discrete Wigner-Ville
%    distribution," IEEE Transactions on Signal Processing, vol. 56, no.
%    11, pp. 5427-5434, Nov. 2008.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(s);

if(nargin>1)
    %If the approximately decimated algorthm is used.
    if(isempty(a))
        a=1;
    end
    
    if(nargin<3||isempty(freqSamps))
        freqSamps=1:N;
    end
    freqSamps=freqSamps-1;%Reduce indices to the 0:(N-1) range.
end

%First, get the analytic signal padded with zeros at the end. If an
%analytic signal is provided, then just zero-pad it.
if(~isreal(s))
    %Just zero-pad z.
    z=[s(:);zeros(N,1)];
else
    %The real2ComplexSignalW is designed to reduce aliasing in Wigner-Ville
    %transforms as compared to the standard approach, which is implemented
    %in the real2ComplexSignal function.
    z=real2ComplexSignalW(s,true);
end
%z is length 2*N; s is length N.

if(nargin>1)%If decimation is to be performed.
    %The decimated algorithm approximation.
    val=DWVDDec(z,N,a,freqSamps);
    
    if(nargout>1)
        t=((0:(2*N/a-1)))*(1/2);
        f=freqSamps*(1/(N-1));
    end
else
    %No decimation.
    val=DWVD(z,N);
    if(nargout>1)
        t=((0:(2*N-1)))*(1/2);
        
        f=linspace(0,1,N);
    end
end

end

function val=DWVD(z,N)
%%DWVD The non-decimated Wigner-Ville transform of Algorithm 4 of [1].
%
%REFERENCES:
%[1] J. M. O'Toole, "Discrete quadratic time-frequency distributions:
%    Definition, computation, and a newborn electroencephalogram
%    application," Ph.D. dissertation, School of Medicine, The University
%    of Queensland, Australia, Mar. 2009.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

N2=length(z);

%Step 1 in Algorithm 4 in [1]. The zero-padding in z interpolates the
%half-frequency bins.
Z=fft(z); 

%Steps 2 and 3 in Algorithm 4 in [1] obtain the Doppler-frequency function.
%It is assumed that the lag-independent kernel is the length of the signal
%(N) and that all of its values in the frequency domain are 1 (in other
%words, no kernel is used).
RC=zeros(N2,N);
span=0:N;
span2=1:(N-1);
for k=0:(N-1)
    %The +1 values are because Matlab indexes things from 1 instead of from
    %0.
    idx1=mod(k+span,N2)+1;
    idx2=mod(k-span,N2)+1; 
  
    %The +1 values in the RC indices are because Matlab indexes things from
    %1, not 0.
    RC(span+1,k+1)=Z(idx1).*conj(Z(idx2));%Step 2
    RC(N2-span2+1,k+1)=conj(RC(span2+1,k+1));%Step 3.
end

%Step 4, transform the Doppler-frequency function Rc to the time-frequency
%domain.
%The ifft value should be real. We are just using the real comment here in
%case any finite precision limitations left the result complex.
val=ifft(RC)/N2;

end

function val=DWVDDec(z,N,a,freqSamps)
%DWVDDEC The decimated Wigner-Ville transform of Algorithm 6 of [1].
%
%REFERENCES:
%[1] J. M. O'Toole, "Discrete quadratic time-frequency distributions:
%    Definition, computation, and a newborn electroencephalogram
%    application," Ph.D. dissertation, School of Medicine, The University
%    of Queensland, Australia, Mar. 2009.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

%The time-decimated DWVD
N2=length(z);
J=length(freqSamps);
L=N2/a;

if(fix(L)~=L)
    error('2*length(x)/a must be an integer.')
end

%Step 1 in Algorithm 6 of [1].
Lh=ceil(L/2);
Z=fft(z);

span0=0:N;
span1=1:(N-1);
spanh0=0:Lh;
spanh1=1:(Lh-1);

%Steps 2, 3  and 4 in Algorithm6 of [1] to form the Doppler-frequency
%function and fold the samples according to the time-decimation.
K=zeros(L,J);
Ktmp=zeros(N2,1);
for k=1:J
    idx1=mod(freqSamps(k)+span0,N2)+1;
    idx2=mod(freqSamps(k)-span0,N2)+1; 

    %Step 2.
    Ktmp(span0+1) = Z(idx1).*conj(Z(idx2));  
    Ktmp(N2-span1+1)=conj(Ktmp(span1+1)); 

    if(a~=1)%If there is decimation in time.
        %THis part is step 3.
        %Fold the Ktmp values.
        x=zeros(Lh+1,1);

        for k1=0:(a-1)
            x(spanh0+1)=x(spanh0+1)+Ktmp(k1*L+spanh0+1);    
        end
        K(spanh0+1,k)=(1/a)*x;

        %Get the negative Doppler values from the positive ones.
        K(L-spanh1+1,k)=conj(K(spanh1+1,k));          
    else
        K(:,k)=Ktmp;
    end
end

%Step 5 in Algorithm 6 of [1].
val=real(ifft(K))/N2;
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
