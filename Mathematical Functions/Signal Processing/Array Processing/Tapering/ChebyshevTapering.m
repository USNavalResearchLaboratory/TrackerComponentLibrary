function w=ChebyshevTapering(N,sidelobedB)
%%CHEBYSHEVTAPERING The Dolph-Chebyshev tapering is a set of amplitude
%                   weights for a discrete linear array where the elements
%                   are spaced a half wavelength apart (narrowband). This
%                   function provides the tapering at a set of values. The
%                   Dolph-Chebyshev tapering makes ALL sidelobes have a
%                   given value. The weights are symmetric about the
%                   origin, which is taken to be the center of the array.
%
%INPUTS:    N The integer number of elements in the array.
%  sidelobedB The number of decibels of the ratio of any of the sidelobe
%             voltages to the main lobe voltage. This must be a negative
%             number. A typical value is -30.
%
%OUTPUTS: w The NX1 set of Dolph-Chebyshev tapering values for each of the
%           elements in the array. w is symmetric. All tapering values are
%           positive and real. The coefficients are normalized so that the
%           largest value is 1.
%
%It is assumed that the positions of the elements in the array are
%1/2 wavelengths apart.
%
%The Dolph-Chebyshev comes from [1]. For odd N, we implement the ifft
%method given in problem 3.4.18 of [2]. For even N, we implement Equation
%3.156 (the method of Chapter 3.4.2).
%
%EXAMPLE:
%Here, we plot how the sidelobes using Chebyshev tapering have uniform
%height.
% nVal=20;
% xPoints=(-(nVal-1)/2):1/2:((nVal-1)/2);
% N=length(xPoints);
% sidelobedB=-20;
% w=ChebyshevTapering(N,sidelobedB);
% T=diag(w);
% [Rsp,U]=standardUVBeamPattern(T,xPoints,[],[],[],[],500);
% 
% figure(1)
% clf
% plot(U,10*log10(Rsp))
% axis([-1 1 -20 20])
% axis square
% h1=xlabel('u');
% h2=ylabel('Array Gain');
% title('Array Power Gain in Decibels')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] C. L. Dolph, "A current distribution for broadside arrays which
%    optimizes the relationship between beam width and side-lobe level,"
%    Proceedings of I.R.E. and Waves and Electrons, vol. 34, no. 6, pp.
%    335-348, Jun. 1946.
%[2] H. L. Van Trees, Optimum Array Processing. New York: Wiley-
%    Interscience, 2002.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(sidelobedB))
   sidelobedB=-30; 
end

R=10^(-sidelobedB/20);

%Equation 3.145
x0=cosh(1/(N-1)*acosh(R));

%For N being odd
if(mod(N,2)==1)
    m=ifftshift((-(N-1)/2):((N-1)/2))';
    p=ChebyshevPoly1(x0*cos(m*pi/N),N-1,-Inf,Inf);
    w=fftshift(ifft(p));
else%N is even
    p=1:(N-1);
    psip=2*acos((1/x0)*cos((2*p-1)*pi/((2*(N-1)))));
    
    V=zeros(N,N);
    V(:,1)=exp(-1j*(N-1)/2)*ones(N,1);
    for k=1:(N-1)
        n=0:(N-1);
        V(:,k+1)=exp(-1j*(N-1)/2)*(exp(1j*n*psip(k)));
    end
    e1=zeros(N,1);
    e1(1)=1;
    w=(V')\e1;
    
    w=(w+flip(w))/2;
end

w=real(w/max(w));
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
