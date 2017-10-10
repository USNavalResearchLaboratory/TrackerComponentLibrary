function w=windowFunSym(N,algorithm,alpha,beta)
%%WINDOWFUNSYM This produces one of a large number of symmetric window
%              functions. Windows functions are used to adjust the
%              sidelobes and Gibbs phenomenon associated with Fourier
%              transforms of signals.
%
%INPUTS:  N The integer length of the transformation.
% algorithm The algorithm to use to generate the window. Some algorithms
%           take an additional parameter, which can be given in alpha.
%           Possible values are:
%           'rectangular' Section VA of [1]. The first sidelobe of a
%                  rectangular window is around -13dB. This is just all
%                  ones.
%           'triangular' Section VB. The alpha parameter specifies which
%                  type of triangular window to use. For alpha=0, the
%                  endpoints are zero just as in [1]; this is equivalent to
%                  the first-order Fejer window of Chapter 2.1.5 of [2].
%                  For alpha=1, the denominator has been modified so that
%                  the endpoints are not zero. Finally, in Section VF of
%                  [1], it is mentioned how Bartlett made this window from
%                  the convolution of two rectangular windows. For odd
%                  values, one can do
%                  conv(ones((N+1)/2,1),ones((N+1)/2,1)). For even values,
%                  one can do conv(ones(N/2+1,1),ones(N/2,1)). For alpha=3,
%                  the window is a scaled version of those convolution
%                  windows and the endpoints are not zero. 
%           'cosAlpha' Section VC of [1]. This is the cosine function
%                  raised to a power. This takes a parameter alpha, which
%                  is typically a small value (e.g. 1-4). For alpha=2, this
%                  is the Hann window. The parameter beta specifies the
%                  type of algorithm. For beta=0, the algorithm in Section
%                  VC is given. The endpoints are zero. For beta=1, the
%                  denominator has been changed so that the endpoints are
%                  not zero. For beta=2, the denominator has been
%                  increased again so that the endpoints are even further
%                  above zero, which is consistent with some common
%                  implementations.
%           'Hamming' Section VD of [1]. This is the Hamming window. The
%                  parameter alpha selects whether the approximate Hamming
%                  window of Equation 30b, which is the most common is
%                  used, or the exact Hamming window of Equation 30a with a
%                  parameter of a=25/46 is used. The sidelobe level is
%                  about -43dB. In both instances, the endpoints are not
%                  zero.
%           'GenHamming' Section VD of [1]. This is the general Hamming
%                  window of Equation 30a using a user-supplied version of
%                  the alpha parameter. The endpoints are not zero.
%           'Blackman' Section VE. This is the Blackman window. The input
%                  alpha selects whether exact or approximate coefficients
%                  are used. For alpha=0, the approximate coefficients are
%                  used. For alpha=1, the exact coefficients are used. If
%                  approximate coefficients are used (the ones most
%                  commonly used as "the" Blackman window), then the
%                  endpoints are zero; otherwise, they are nonzero. The
%                  first sidelobe is about -51dB down. This window uses
%                  three terms in the sums of the optimization problem in
%                  Equation 31.
%           'Blackman-Harris-Min-Approx' Section VE of [1]. These are the
%                  Blackman- Harris windows with coefficients from the
%                  Table on page 186. The parameter alpha selects which
%                  parameters to use. For alpha=0, the approximate
%                  three-term minimum sidelobe level value of -67dB is
%                  used; alpha=1 is a 3-term -61dB formula. For alpha=1,
%                  the  approximate four-term minimum sidelobe level
%                  formula is used for -92dB sidelobes, and alpha=3 gives a
%                  four-term -74dB formula. 
%           'Riesz' Section VF1 of [1]. This is the Riesz (Parzen, Bochner)
%                  window. The first sidelobe is -22dB. If alpha=0, then
%                  the method of Section VF1 is used and the endpoints are
%                  zero. For alpha=1, the denominator is modified so that
%                  the endpoints are not zero.
%           'Riemann' Section VF2 of [1]. This is the Riemann window. The
%                  input alpha selects which version to use. This window is
%                  the central lobe of a sinc kernel. For alpha=0, the
%                  version as in Section VF2 is used and the endpoints are
%                  zero. For alpha=1, a version with a modified denominator
%                  is used so that the endpoints are not zero.
%           'Poussin' Section VF3. This is the de la Vallé-Poussin window.
%                  For alpha=0, the algorithm of Section VF3 is sued and
%                  the endpoints are zero. For alpha=1, the denominators
%                  are modified so that the endpoints are not zero.
%           'Tukey' Section VF4 of [1]. This is the Tukey window, also
%                  known as the cosine tapered window. It takes a scalar
%                  parameter alpha, which can vary from 0 to 1. The value
%                  pi in Equation 38 should be multiplied by 2. The input
%                  alpha specifies the fraction of either half of the pulse
%                  that is 1. This is essentially a rectangular pulse with
%                  modified edges. The option beta specified the algorithm
%                  to use. For beta=0, the algorithm of Section VF4 is used
%                  and the endpoints are zero. For beta=1, the denominators
%                  are modified so that the endpoints are not zero.
%           'Bohman' Section VF5. This is the Bohman window. For alpha=0,
%                  the algorithm of Section VF5 is used and the endpoints
%                  are zero. For alpha=1, the denominators are modified so
%                  that the endpoints are not zero.
%           'Poisson' Section VF6 of [1]. This is the Poisson window. It is
%                  a two- sided exponential and it takes a positive
%                  parameter alpha. The endpoints are not zero.
%           'Hanning-Poisson' Section VF7 of [1]. The Hanning-Poisson
%                  window. It takes a parameter alpha that is similar to
%                  that used in the Poisson window. The input beta selects
%                  the type of algorithm. For beta=0, the algorithm as
%                  given in Section VF7 with zeros at the endpoints is
%                  used. For beta=1, the method is modified so that the
%                  endpoints are not zero.
%           'Cauchy' Section VF8 of [1]. The Cauchy window. It takes the
%                  parameter alpha. The endpoints are not zero.
%           'Gaussian' Section VG of [1]. The Gaussian (Weierstrass)
%                  window. This takes the parameter alpha. The parameter
%                  alpha is related to the "standard deviation" of the
%                  Gaussian as sigma=(N-1)/(2*alpha). The endpoints are not
%                  zero. Increasing alpha decreases the width of the
%                  window.
%           'Kaiser' Section VI. The Kaiser-Bessel window. This takes an
%                  input parameter alpha. The end points are not zero.
%           'Fejer-2' %This is the second-order Fejer filter of Chapter
%                   2.1.5 of [2]. The endpoints are not zero.
%           'Nuttall'Section C7 of [3]. These are Nuttall windows. The
%                   parameter alpha selects the window. Possible values
%                   are:
%                   0 -46.7dB, Nuttall3
%                   1 -64.2dB, Nuttall3a
%                   2 -71.5dB, Nuttall3b
%                   3 -60.9dB, Nuttall4
%                   4 -82.6dB, Nuttall4a
%                   5 -93.3dB, Nuttall4b
%                   6 -98.1dB, Nuttall4c
%           'FlatTop' Sections D1-D3 of [3]. The symmetric flat-top window.
%                   The parameter alpha selects the window. Possible values
%                   are:
%                   0 -31.7dB, SFT3F
%                   1 -44.7dB, SFT4F
%                   2 -57.3dB, SFT5F
%                   3 -44.2dB, SFT3M
%                   4 -66.5dB, SFT4M
%                   5 -89.9dB, SFT5M
%                   6 -44dB, FTNI
%                   7 -70.4dB, FTHP
%                   8 -76.6dB, FTSRS
%                   9 -70.4dB, HFT70
%                   10 -95dB, HFT95
%                   11 -90.2dB, HFT90D
%                   12 -116.8dB, HFT116D
%                   13 -144.1dB, HFT144D
%                   14 -169.5dB, HFT169D
%                   15 -196.2dB,HFT196D
%                   16 -223dB, HFT223D
%                   17 -248.4dB, HFT248D
%            'Chebyshev' This is the Chebyshev window. Alpha is the number
%                   of decibels down the sidelobe should be and must be
%                   negative. This just calls the ChebyshevTapering
%                   function.
%            'Taylor' The Taylor window. This type of window is most
%                   commonly used with antenna arrays. The input alpha is
%                   the number of sidelobes to push down and the input beta
%                   is a negative number indicating the number of decibels
%                   they should try to be pushed down. Taylor tapering
%                   samples a continuous tapering. This just calls the
%                   function TaylorTapering.
% alpha, beta These inputs are parameters for the various algorithms. The
%             descriptions of the algorithms specifies when these are
%             necessary.          
%
%OUTPUTS: w This is an NX1 real vector containing the weights of the
%           window.
%
%EXAMPLE:
%Here we show how the FFT of a Fejer-2 window has no sidelobes, whereas
%that of a rectangular windows has many sidelobes. We zero-pad the end of
%the signal before the FFT to interpolate in the Fourier domain. The
%zero-freuqnecy is in the center.
% N=100;
% padding=1000;
% wRect=windowFunSym(N,'rectangular');
% wFejer=windowFunSym(N,'Fejer-2');
%
% WRectdB=fftshift(20*log10(abs(fft(wRect,N+padding))));
% WFejerdB=fftshift(20*log10(abs(fft(wFejer,N+padding))));
% %Normalize
% WRectdB=WRectdB-max(WRectdB);
% WFejerdB=WFejerdB-max(WFejerdB);
%
% figure(1)
% clf
% hold on
% plot(WRectdB,'-k','linewidth',2)
% plot(WFejerdB,'-r','linewidth',2)
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Rectangular','Fejer')
% axis([0, N+padding, -50, 0])
%
%REFERENCES:
%[1] F. J. Harris, "On the use of windows for harmonic analysis with the
%    discrete Fourier transform," Proceedings of the IEEE, vol. 66, no. 1,
%    pp. 172-204, Jan. 1978.
%[2] W. Wasylkiwskyj, Signal and Transforms in Linear Systems Analysis.
%    Heidelberg: Springer, 2013.
%[3] G. Heinzel, A. Rüdinger, and R. Schilliung, "Spectrum and spectral
%    density estimation by the discrete fourier transform (DFT), including
%    a comprehensive list of window functions and some new flat-top
%    windows," Max-Planck-Institut für Gravitationsphysik (Albert Einstein
%    Institut)), Hanover, Germany, Tech. Rep., 15 Feb. 2002.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(alpha))
    alpha=2;%0.75;
end

if(nargin<2||isempty(algorithm))
   algorithm='triangular'; 
end

if(mod(N,2))%If odd
    n=(-(N-1)/2):((N-1)/2);
else%If even
    n=(-N/2):(N/2-1);
    n=n+1/2;
end

switch(algorithm)
    case 'rectangular'%Section VA of [1].
        w=ones(N,1);
    case 'triangular'%Section VB of [1].
        switch(alpha)
            case 0
                w=1-abs(n)/((N-1)/2);
            case 1
                w=1-abs(n)/(N/2);
            case 2
                if(mod(N,2))
                    w=1-abs(n)/((N+1)/2); 
                else
                    w=1-abs(n)/(N/2);
                end
            otherwise
                error('Unknown option specified')
        end
    case 'cosAlpha'%Section VC of [1].
        switch(beta)
            case 0
                w=cos(pi*n/(N-1)).^alpha;
            case 1
                w=cos(pi*n/N).^alpha;
            case 2
                w=cos(pi*n/(N+1)).^alpha;
            otherwise
                error('Unknown option specified')
        end
    case 'Hamming'%Section VD of [1].
        switch(alpha)
            case 0
                w=0.54+0.46*cos(2*pi*n/(N-1));
            case 1
                a=25/46;
                w=a+(1-a)*cos(2*pi*n/(N-1));
            otherwise
                error('Unknown option specified')
        end        
    case 'GenHamming'%Section VD of [1].
        w=alpha+(1-alpha)*cos(2*pi*n/(N-1));
    case 'Blackman'%Section VE of [1].
                   switch(alpha)
                       case 0%Approximate coefficients
                            w=0.42+0.5*cos(2*pi/(N-1)*n)+0.08*cos(2*pi/(N-1)*2*n);
                       case 1%Exact coefficients
                            a0=7938/18608;
                            a1=9240/18608;
                            a2=1430/18608;
                            w=a0+a1*cos(2*pi/(N-1)*n)+a2*cos(2*pi/(N-1)*2*n);
                       otherwise
                           error('Unknown option specified')
                   end
    case 'Blackman-Harris-Min-Approx'%Section VE of [1].
        switch(alpha)
            case 0%-67dB
                a0=0.42323;
                a1=0.49755;
                a2=0.07922;
                a3=0;
            case 1%-61dB
                a0=0.44959;
                a1=0.49364;
                a2=0.05677;
                a3=0;
            case 2%-92dB
                a0=0.35875;
                a1=0.48829;
                a2=0.14128;
                a3=0.01168;
            case 3
                a0=0.40217;
                a1=0.49703;
                a2=0.09392;
                a3=0.00183;
            otherwise
                error('The number of terms specified is not supported')
        end
        
        w=a0+a1*cos(2*pi/(N-1)*n)+a2*cos(2*pi/(N-1)*2*n)+a3*cos(2*pi/(N-1)*3*n);
    case 'Riesz'%Section VF1 of [1].
        switch(alpha)
            case 0
                w=1-(n/((N-1)/2)).^2;
            case 1
                w=1-(n/(N/2)).^2;
            otherwise
                error('Unknown option specified')
        end
    case 'Riemann'%Section VF2 of [1].
        switch(alpha)
            case 0
                w=sincFun(2*n/(N-1));
            case 1
                w=sincFun(2*n/(N)); 
            otherwise
                error('Unknown option specified')
        end  
    case 'Poussin'%Section VF3 of [1].
        switch(alpha)
            case 0
                sel=abs(n)<=((N-1)/4);
                w=zeros(N,1);
                w(sel)=1-6*(n(sel)/((N-1)/2)).^2.*(1-abs(n(sel))/((N-1)/2));
                w(~sel)=2*(1-abs(n(~sel))/((N-1)/2)).^3;
            case 1
                sel=abs(n)<=((N)/4);
                w=zeros(N,1);
                w(sel)=1-6*(n(sel)/(N/2)).^2.*(1-abs(n(sel))/(N/2));
                w(~sel)=2*(1-abs(n(~sel))/(N/2)).^3;
            otherwise
                error('Unknown option specified')
        end
    case 'Tukey'%Section VF4 of [1].
         switch(beta)
             case 0
                w=ones(N,1);
                sel=abs(n)>alpha*(N-1)/2;

                w(sel)=(1/2)*(1+cos(2*pi*(abs(n(sel))-alpha*(N-1)/2)/(2*(1-alpha)*(N-1)/2)));
             case 1
                w=ones(N,1);
                sel=abs(n)>alpha*N/2;

                w(sel)=(1/2)*(1+cos(2*pi*(abs(n(sel))-alpha*N/2)/(2*(1-alpha)*N/2)));
             otherwise
                error('Unknown option specified')
         end
    case 'Bohman'%Section VF5 of [1].
        switch(alpha)
            case 0
                w=(1-abs(n)/((N-1)/2)).*cos(pi*abs(n)/((N-1)/2))+(1/pi)*sin(pi*abs(n)/((N-1)/2));
            case 1
                w=(1-abs(n)/(N/2)).*cos(pi*abs(n)/(N/2))+(1/pi)*sin(pi*abs(n)/(N/2));
            otherwise
                error('Unknown option specified')
        end
    case 'Poisson'%Section VF6 of [1].
        w=exp(-alpha*abs(n)/((N-1)/2));
    case 'Hanning-Poisson'%Section VF7 of [1].
        switch(beta)
            case 0
                w=(1/2)*(1+cos(pi*n/((N-1)/2))).*exp(-alpha*abs(n)/((N-1)/2));
            case 1
                w=(1/2)*(1+cos(pi*n/(N/2))).*exp(-alpha*abs(n)/(N/2));
            otherwise
                error('Unknown option specified')
        end
    case 'Cauchy'%Section VF8 of [1].
        w=1./(1+(alpha*n/((N-1)/2)).^2);
    case 'Gaussian'%Section VG of [1].
        w=exp(-(1/2)*(alpha*n/((N-1)/2)).^2);
    case 'Kaiser'%Section VI of [1].
        val=alpha*sqrt(1-(n/((N-1)/2)).^2);
        w=besseli(0,val)/besseli(0,alpha);
    case 'Fejer-2'%This is the second-order Fejer filter of Chapter 2.1.5
        %of [2]. The endpoints are not zero.

        M=(N-1)/2;

        w=zeros(N,1);
        for idx=1:N
            k=n(idx);
            
            nList=1:fix(M-abs(k)+1);
            w(idx)=(1/(M+1))*sum(nList./(abs(k)+nList));
        end
    case 'Nuttall'%Section C7 of [3]. The Nuttall windows. The signs of the
        %terms have been corrected.
        switch(alpha)
            case 0%-46.7dB, Nuttall3
                c0=0.375;
                c1=0.5;
                c2=0.125;
                c3=0;
            case 1%-64.2dB, Nuttall3a
                c0=0.40897;
                c1=0.5;
                c2=0.09103;
                c3=0;
            case 2%-71.5dB, Nuttall3b
                c0=0.4243801;
                c1=0.4973406;
                c2=0.0782793;
                c3=0;
            case 3%-60.9dB, Nuttall4
                c0=0.3125;
                c1=0.46875;
                c2=0.1875;
                c3=0.03125;
            case 4%-82.6dB, Nuttall4a
                c0=0.338946;
                c1=0.481973;
                c2=0.161054;
                c3=0.018027;
            case 5%-93.3dB, Nuttall4b
                c0=0.355768;
                c1=0.487396;
                c2=0.144232;
                c3=0.012604;
            case 6%-98.1dB, Nuttall4c
                c0=0.3635819;
                c1=0.4891775;
                c2=0.1365995;
                c3=0.0106411;
            otherwise
                error('Unknown type of Nuttall window selected')
        end
        w=c0+c1*cos(2*pi*n/N)+c2*cos(4*pi*n/N)+c3*cos(6*pi*n/N);
    case 'FlatTop'%Section D1-D3 of [3]. The symmetric flat-top window. The
        %signs of the coefficients in the report were not all correct; they
        %have been corrected.
        c3=0;
        c4=0;
        c5=0;
        c6=0;
        c7=0;
        c8=0;
        c9=0;
        c10=0;
        switch(alpha)
            case 0%-31.7dB, SFT3F
                c0=0.26526;
                c1=0.5;
                c2=0.23474;
            case 1%-44.7dB, SFT4F
                c0=0.21706;
                c1=0.42103;
                c2=0.28294;
                c3=0.07897;
            case 2%-57.3dB, SFT5F
                c0=0.1881;
                c1=0.36923;
                c2=0.28702;
                c3=0.13077;
                c4=0.02488;
            case 3%-44.2dB, SFT3M
                c0=0.28235;
                c1=0.52105;
                c2=0.19659;
            case 4%-66.5dB, SFT4M
                c0=0.241906;
                c1=0.460841;
                c2=0.255381;
                c3=0.041872;
            case 5%-89.9dB, SFT5M
                c0=0.209671;
                c1=0.407331;
                c2=0.281225;
                c3=0.092669;
                c4=0.0091036;
            case 6%-44dB, FTNI
                c0=0.2810639;
                c1=0.5208972;
                c2=0.1980399;
            case 7%-70.4dB, FTHP
                c0=1;
                c1=1.912510941;
                c2=1.079173272;
                c3=0.1832630879;
            case 8%-76.6dB, FTSRS
                c0=1;
                c1=1.93;
                c2=1.29;
                c3=0.388;
                c4=0.028;
            case 9%-70.4dB, HFT70
                c0=1;
                c1=1.90796;
                c2=1.07349;
                c3=0.18199;
            case 10%-95dB, HFT95
                c0=1;
                c1=1.9383379;
                c2=1.3045202;
                c3=0.4028270;
                c4=0.0350665;
            case 11%-90.2dB, HFT90D
                c0=1;
                c1=1.942604;
                c2=1.340318;
                c3=0.440811;
                c4=0.043097;
            case 12%-116.8dB, HFT116D
                c0=1;
                c1=1.9575375;
                c2=1.4780705;
                c3=0.6367431;
                c4=0.1228389;
                c5=0.0066288;
            case 13%-144.1dB, HFT144D
                c0=1;
                c1=1.96760033;
                c2=1.57983607;
                c3=0.81123644;
                c4=0.22583558;
                c5=0.02773848;
                c6=0.00090360;
            case 14%-169.5dB, HFT169D
                c0=1;
                c1=1.97441842;
                c2=1.65409888;
                c3=0.95788186;
                c4=0.33673420;
                c5=0.06364621;
                c6=0.00521942;
                c7=0.00010599;
            case 15%-196.2dB,HFT196D
                c0=1;
                c1=1.979280420;
                c2=1.710288951;
                c3=1.081629853;
                c4=0.448734314;
                c5=0.112376628;
                c6=0.015122992;
                c7=0.000871252;
                c8=0.000011896;
            case 16%-223dB, HFT223D
                c0=1;
                c1=1.98298997309;
                c2=1.75556083063;
                c3=1.19037717712;
                c4=0.56155440797;
                c5=0.17296769663;
                c6=0.03233247087;
                c7=0.00324954578;
                c8=0.00013801040;
                c9=0.00000132725;
            case 17%-248.4dB, HFT248D
                c0=1;
                c1=1.985844164102;
                c2=1.791176438506;
                c3=1.282075284005;
                c4=0.667777530266;
                c5=0.240160796576;
                c6=0.056656381764;
                c7=0.008134974479;
                c8=0.000624544650;
                c9=0.000019808998;
                c10=0.000000132974;
            otherwise
                error('Unknown type of flat top window selected')
        end
        
        x=n/N;
        w=c0+c1*cos(2*pi*x)+c2*cos(4*pi*x)+c3*cos(6*pi*x)+c4*cos(8*pi*x)+c5*cos(10*pi*x)+c6*cos(12*pi*x)+c7*cos(14*pi*x)+c8*cos(16*pi*x)+c9*cos(18*pi*x)+c10*cos(20*pi*x);
    case 'Chebyshev'
        w=ChebyshevTapering(N,alpha);
    case 'Taylor'
        w=TaylorTapering(alpha,beta,n/2);
    otherwise
        error('Unknown window selected')
end

w=w(:);%Make it a column vector.

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
