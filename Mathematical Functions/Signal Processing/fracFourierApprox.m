function fracVal=fracFourierApprox(f,a)
%%FRACTFOURIERAPPROX Compute an approximate discrete fractional Fourier
%           transformation using a fast approximation that is valid for
%           fractional parameter values from -1<=p<=1.
%
%INPUTS: f The NXnumEls matrix of real or complex sequences over which the
%          fractional Fourier transform is desired. The transform is
%          performed over each column.
%        a The order of the fractional Fourier transform.
%
%OUTPUTS: fracVal The NXnumEls approximate values of the fractional Fourier
%                 transform over each column of f. 
%
%This function implements the algorithm described at the end of [1].
%
%The true discrete fractional Fourier transform is given by the function
%fracFourier, which this approximates. We will note that for a column
%vector x,
%fft(x)==ifftshift(fracFourier(fftshift(x),1))*sqrt(length(x))
%ifft(x)==ifftshift(fracFourier(fftshift(x),3))/sqrt(length(x))
%and
%fracFourier(x,1)=fftshift(fft(ifftshift(x)))/sqrt(length(x))
%fracFourier(y,3)=fftshift(ifft(ifftshift(y)))*sqrt(length(y))
%to within finite precision limits. That is, the discrete fractional
%Fourier transform is centered differently than the fft.
%
%EXAMPLE:
%Here, we consider how well the approximation compres to the exact
%algorithm when considering a chirp.
% T=1e-3;
% fStart=0;
% fEnd=10e3;
% tSamp=1/(50*2*fEnd);%50 times the Nyquist frequency.
% [s,t]=LFMChirp(T,fStart,fEnd,tSamp);
% s=s(:);
% a=0.5;
% fApprox=fracFourierApprox(s,a);
% fDisc=fracFourier(s,a,length(s)-1);
% 
% figure(1)
% clf
% hold on
% plot(abs(fDisc),'-k','linewidth',4)
% plot(abs(fApprox),'-r','linewidth',1)
% h2=ylabel('Magnitude');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Discrete Transform','Approximate Transform')
% 
% figure(2)
% clf
% hold on
% plot(angle(fDisc),'-k','linewidth',4)
% plot(angle(fApprox),'-r','linewidth',1)
% h2=ylabel('Phase (radians)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Discrete Transform','Approximate Transform')
%
%REFERENCES:
%[1] J. García, D. Mas, and R. G. Dorsch, "Fractional-Fourier-transform
%   calculation through the fast-Fourier-transform algorithm," Applied
%   Optics, vol. 35, no. 35, pp. 7013-7018, 10 Dec. 1996.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(f,1);
numEls=size(f,2);

m=((-N/2):(N/2-1)).';

fracVal=zeros(N,numEls);
for curEl=1:numEls
    x=f(:,curEl);

    %Step 1
    x=x.*exp(-1j*pi/N*tan(a*pi/4)*m.^2);

    %Step 2a. The definition of the fft used in the paper is not the same
    %as in the fft function. Thus, we must convert when doing the fft for
    %step 2a.
    x=fftshift(fft(ifftshift(x)))/sqrt(N);

    %Step 2b
    x=x.*exp(-1j*pi/N*sin(a*pi/2)*m.^2);

    %Step 2c
    x=fftshift(ifft(ifftshift(x)))*sqrt(N);

    %Step 3
    x=x.*exp(-1j*pi/N*tan(a*pi/4)*m.^2);

    %Equation 18. In the book, the i term before the pi/4 was missing.
    Mp=exp(-1i*pi*sign(sin(a*pi/2))/4+1i*a*pi/4+1i*pi/4);
    fracVal(:,curEl)=Mp*x;
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
