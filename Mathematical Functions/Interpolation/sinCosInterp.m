function [fInterp,c]=sinCosInterp(f,tList,deltaT,algorithm,numDeriv)
%%SINCOSINTERP Given uninformely spaced samples of a function, interpolate
%              the values of a function or the values of the first
%              derivative of a function over a desired region.
%              Trigonometric interpolation using sinusoid or cosine
%              formulas is used.
%
%INPUT: f An NX1 or 1XN vector of points that are uniformly spaced in the
%         independent variable.
%   tList The NX1 or 1XN set of values of the independent variable at which
%         the function or its derivative should be interpolated. This
%         function assumes that the first element in f is evaluated at t=0
%         and subsequent elements are spaced deltaT apart. Interpolated
%         points outside of the values given as references tend to be bad.
%  deltaT The spacing of the samples in f in terms of the independent
%         variable t. If this parameter is omitted or an empty matrix is
%         passed, the default of 1 is used.
% algorithm Select the interpolation algorithm. Possible variables are:
%         0 Use sinusoid interpolation on Chapter 2.1.8 of [1]. This tends
%           to be oscillatory at the edges of the interpolation region.
%         1 Use the cosine algorithm of Chapter 2.1.8 of [1]. This is the
%           default if this parameter is omitted or an empty matrix is
%           passed.
% numDeriv This can be 0 or 1. This is the number of derivatives to take.
%
%OUTPUTS: fInterp The interpolated function values of the derivatives
%                 corresponding to the values in tList. This has the same
%                 dimensions as tList.
%
%This implements the sinusoid interpolation algorithm in Equation 2.74 of
%Chapter 2.1.8 and the cosine interpolation algorithm in Equation 2.83 of
%Chapter 2.1.8.
%
%EXAMPLE 1:
%Here, we interpolate the values of exp(-t) from 0 to 1 using 11 samples,
%similar to the example in Chapter 2.1.8 of [1]. One will see that the
%interpolation due to algorithm 1 is much better than algorithm 0.
%Algorithm 0 oscillates much more.
% N=11;
% tSamp=linspace(0,1,N);
% deltaT=tSamp(2)-tSamp(1);
% tDes=linspace(0,1,1000*N);
% 
% fSamp=exp(-tSamp);
% fTrue=exp(-tDes);
% 
% fInterp0=sinCosInterp(fSamp,tDes,deltaT,0);
% fInterp1=sinCosInterp(fSamp,tDes,deltaT,1);
% 
% figure(1)
% clf
% hold on
% plot(tDes,fTrue,'-k','linewidth',4)
% plot(tDes,fInterp0,'--r','linewidth',2)
% plot(tDes,fInterp1,'--g','linewidth',2)
% h1=xlabel('t');
% h2=ylabel('f');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Truth', 'Algorithm 0', 'Algorithm 1')
%
%EXAMPLE 2:
%Here, we interpolate the derivative of the function. One will see that
%algorithm 0 still has terrible oscillations. Algorithm 1 is better, but
%still is not as stable as one might desire. If the signal were upsampled
%so that none of the samples hit any of the training samples, then
%algorithm 1 performs better.
% N=100;
% tSamp=linspace(0,2*pi-0.1,N);
% deltaT=tSamp(2)-tSamp(1);
% tDes=linspace(0,2*pi-0.1,N);
% 
% fSamp=sin(tSamp);
% fTrueDeriv=cos(tDes);
% 
% fInterpDeriv0=sinCosInterp(fSamp,tDes,deltaT,0,1);
% fInterpDeriv1=sinCosInterp(fSamp,tDes,deltaT,1,1);
% 
% figure(1)
% clf
% hold on
% plot(tDes,fTrueDeriv,'-k','linewidth',4)
% plot(tDes,fInterpDeriv0,'-r','linewidth',2)
% plot(tDes,fInterpDeriv1,'--g','linewidth',2)
% h1=xlabel('t');
% h2=ylabel('f');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Truth', 'Algorithm 0', 'Algorithm 1')
%
%REFERENCES:
%[1] W. Wasylkiwskyj, Signal and Transforms in Linear Systems Analysis.
%    Heidelberg: Springer, 2013.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(f);

if(nargin<3||isempty(deltaT))
   deltaT=1; 
end

if(nargin<4||isempty(algorithm))
   algorithm=1;%Cosine interpolation is the default. 
end

if(nargin<5||isempty(numDeriv))
   numDeriv=0; 
end

numT=length(tList);
fInterp=zeros(size(tList));

switch(numDeriv)
    case 0%Evaluate the functions directly.
        switch(algorithm)
            case 0%The sinusoid interpolation in Equation 2.74.
                l=0:(N-1);
                for tCur=1:numT
                    t=tList(tCur);
                    c=sin(pi*(t/deltaT-l))./(N*sin(pi/N*(t/deltaT-l)));
                    %Assume that any NaN values are due to t/deltaT=l. In
                    %such an instance, we insert the limit:
                    c(isnan(c))=1;
                    fInterp(tCur)=sum(c(:).*f(:));
                end
            case 1
                %Cosine interpolation as Equation 2.82
                m=(0:(N-1)).';
                %epsilon_m, eM, is defined before Equation 2.61 in Section
                %2.1.7. This is the epsilon_m/2 term from Equation 2.82.
                eM2=ones(N,1);
                eM2(1)=1/2;
                coeff=1/(N-1/2);
                for tCur=1:numT
                    t=tList(tCur);

                    fInterp(tCur)=sum(coeff*eM2.*f(:).*(1+kM(t/deltaT-m,N)+kM(t/deltaT+m,N)));
                end
            otherwise
                error('Unknown algorithm selected')
        end
    case 1%Interpolate the first derivative.
        switch(algorithm)
            case 0%The sinusoid interpolation in Equation 2.74.
                l=0:(N-1);
                for tCur=1:numT
                    t=tList(tCur);
                    x=t/deltaT-l;
                    c=(pi*csc((pi*x)/N).*(N*cos(pi*(-x))+cot((pi*x)/N).*sin(pi*(-x))))/(N^2*deltaT);
                    %Assume that any non-finite values are due to
                    %t/deltaT=l. In such an instance, we insert the limit:
                    c(~isfinite(c))=0;
                    fInterp(tCur)=sum(c(:).*f(:));
                end
            case 1
                %Cosine interpolation as Equation 2.82
                m=(0:(N-1)).';
                %epsilon_m, eM, is defined before Equation 2.61 in Section
                %2.1.7. This is the epsilon_m/2 term from Equation 2.82.
                eM2=ones(N,1);
                eM2(1)=1/2;
                coeff=1/(deltaT*(N-1/2));
                for tCur=1:numT
                    t=tList(tCur);

                    fInterp(tCur)=sum(coeff*eM2.*f(:).*(kMDeriv(t/deltaT-m,N)+kMDeriv(t/deltaT+m,N)));
                end
            otherwise
                error('Unknown algorithm selected')
        end
    otherwise
        error('Unsupported derivative requested')
end
end

function val=kM(t,M)
%This function implements Equation 2.83 for cosine interpolation.

    val=cos(pi*M/(2*M-1)*t).*(sin(pi*((M-1)/(2*(M-1/2)))*t)./sin((pi/(2*(M-1/2)))*t));
    
    %Assume that any NaN values are due to t=0. In such an instance, we
    %insert the limit:
    val(isnan(val))=M-1;
end

function val=kMDeriv(t,M)
%This function implements the derivative of Equation 2.83 for cosine
%interpolation.
    val=(pi*csc((pi*t)/(-1+2*M)).*((-1+2*M)*cos(pi*t)-cos((pi*t)/(-1+2*M))-2*cos((M*pi*t)/(-1+2*M)).*cot((pi*t)/(-1+2*M)).*sin(((-1+M)*pi*t)/(-1+2*M))))/(-2+4*M);

    %Assume that any non-finite values are due to t=0. In such an instance,
    %we insert the limit:
    val(~isfinite(val))=0;
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
