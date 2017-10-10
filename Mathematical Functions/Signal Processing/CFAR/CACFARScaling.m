function data2DSqrScaled=CACFARScaling(data2D,numGuardCells,numAvgCells,PFA,method)
%%CACFARSCALING In cell-averaging CFAR, the detection threshold at one cell
%           depends on a sum of values in the surrounding cells. Here, the
%           magnitude squared of data2D is taken and then scaled based upon
%           the CACFAR detection threshold such that detection scan now be
%           performed using a threshold of 1. This is often done prior to
%           noncoherently adding together multiple range-Doppler maps as
%           one tends to get performance doing CFAR scaling before finding
%           detections rather than performing CFAR on the output of the
%           non-coherent integrator.
%
%INPUTS: data2D The delay-Doppler plot or a set of delay-Doppler plots.
%               This is a numRowsXnumColsXnumPlots set of numPlots 2D
%               range-Doppler maps (or similar matrices on which CFAR
%               should be performed). This can contain complex values.
% numGuardCells The integer number of cells in each dimension around each
%               test cell that are not considered in the average used for
%               determining the changing detection threshold. This is a 2X1
%               vector. If the same value is used in both directions, then
%               a scalar can be passed.
%   numAvgCells The width of the region in cells in each dimension after
%               the guard cell region that define the average used for
%               determinig the threshold. This is a 2X1 vector. If the same
%               value is used in both directions, then a scalar can be
%               passed.
%           PFA The scalar probability of false alarm between 0 and 1 that
%               determines the threshold for detection.
%        method The method used. This is an optional value from 0 to 2 that
%               does not change the output of the function, only how the
%               results are computed. Value 0 uses the most intuitive
%               method involving two-dimensional convolutions. Methods 1
%               and 2 are more complicated and try to break the convolution
%               regions into parts to speed up the problem.
%
%OUTPUTS: data2DSqrScaled This is the magnitude squared of data2D, scaled
%               such that one can use a detection threshold of 1 at every
%               point.
%
%Cell averaging CFAR is discussed in [1]. To determine whether a value in a
%particular cell should be flagged as a detection, rather than comparing
%its value to a fixed constant depending on the signal noise power, the
%threshold is based on a sum of cells surrounding the test cell, outside of
%a guard interval. A rectangle of cells is used here.
%
%EXAMPLE:
%In this example, we make and dispay a range-Doppler map, perform CA-CFAR
%scaling and display the result, then do detection and then display the
%results.
% fc=1e9;%1GHz carrier frequency.
% B=2e6;%2Mhz bandwidth.
% %Baseband start and end frequencies.
% fStart=-B/2;
% fEnd=B/2;
% %Sampling rate is two times the Nyquist rate.
% T0=1/(2*2*fEnd);%Sampling period in seconds.
% T=2e-5;%Chirp duration in seconds.
% 
% PRF=2000;%Pulse repetition frequency (Hertz)
% TB=1/PRF;%The pulse repetition period.
% %The number of samples per PRI. The above parameters were chosen so that
% %this is an integer. Fix just deals with finite precision errors.
% Ns=fix(TB/T0);
% 
% %Generate the reference signal. This is an up-chirp.
% x=LFMChirp(T,fStart,fEnd,T0);
% x=x(:);
% 
% %We will use 64 pulse repetition intervals.
% NB=64;
% 
% %True target parameters.
% c=Constants.speedOfLight;
% rTrue=[40e3;50e3;60e3];
% tau=rTrue/c;%The true delay (s).
% 
% %The true range rate (m/s).
% rrTrue=[100;60;-100];
% a=rrTrue/c;
% 
% %Complex amplitudes
% A=[48;512*exp(1j*2*pi*rand(1));32*exp(1j*2*pi*rand(1))];
% 
% numTargets=length(rTrue);
% 
% %Allocate space for the received signal. The first dimensions is "fast
% %time"; the second dimension if "slow time".
% y=zeros(Ns,NB);
% 
% %Create the received signal, properly delayed and Doppler shifted for
% %each PRI. The same waveform is used in each PRI.
% t=0:T0:((Ns-1)*T0);%Sample times
% for i=0:(NB-1)
%     
%     for curTar=1:numTargets
%         tCur=t-a(curTar)*t-tau(curTar)-i*TB;
%         %The signal simulated with range migration.
%         y(:,i+1)=y(:,i+1)+(A(curTar)*exp(-1j*2*pi*fc*(tau(curTar)+a(curTar)*t)).*LFMChirp(T,fStart,fEnd,tCur)).';
%     end
% 
%     y(:,i+1)=y(:,i+1)+ComplexGaussianD.rand(Ns).';
%     
%     t=t+TB;%Increment to the next time step.
% end
% 
% M1=1;
% M2=1;
% 
% %Windowing in range and Doppler to lower sidelobes.
% Nx=length(x);
% wRange=windowFunSym(Nx,'Blackman',0);
% x=wRange.*x;
% wDoppler=windowFunSym(NB,'Nuttall',1);
% 
% [delayDopPlot,Doppler,delay]=delayDopplerPlotNBPulseDop(x,y,M1,M2,wDoppler,T0);
% 
% numGuardCells=[2;4];
% numAvgCells=[5;3];
% PFA=1e-7;
% 
% range=c*delay;
% rangeRate=Doppler*(c/fc);
% 
% %Display the range-Doppler plot
% figure(1)
% clf
% imagesc([rangeRate(1),rangeRate(end)],[range(1), range(end)]/1e3,20*log10(abs(delayDopPlot)));
% set(gca,'YDir','normal')
% caxis([-30 27])
% colormap(jet(256))
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% delayDopPlotMag2Scaled=CACFARScaling(delayDopPlot,numGuardCells,numAvgCells,PFA);
% 
% %Display the scaled range-Doppler plot
% figure(2)
% clf
% imagesc([rangeRate(1),rangeRate(end)],[range(1), range(end)]/1e3,10*log10(delayDopPlotMag2Scaled));
% set(gca,'YDir','normal')
% caxis([-30 27])
% colormap(jet(256))
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% [indexR,indexD]=find(delayDopPlotMag2Scaled>1);
% DetectionList(1).Index=[indexR.';indexD.'];
% numRows=size(delayDopPlotMag2Scaled,1);
% numCols=size(delayDopPlotMag2Scaled,2);
% idx=sub2ind([numRows,numCols],indexR,indexD);
% DetectionList(1).Value=delayDopPlotMag2Scaled(idx);
% 
% idx=DetectionList(1).Index;
% 
% rVals=range(idx(1,:));
% DopVals=rangeRate(idx(2,:));
% 
% %Display the CFAR detections.
% figure(3)
% clf
% hold on
% scatter(DopVals,rVals)
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% axis([rangeRate(1), rangeRate(end), range(1), range(end)])
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%as a modification of the function CACFAR by Tom Higgins and Daniel
%Scholnick.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(method))
    method=0; 
end

if(length(numGuardCells)==1)
    numGuardCells=[numGuardCells;numGuardCells];
end

if(length(numAvgCells)==1)
    numAvgCells=[numAvgCells;numAvgCells];
end

NFilter=numGuardCells+numAvgCells;

numRows=size(data2D,1);
numCols=size(data2D,2);
numPlots=size(data2D,3);
data2DSqrScaled=zeros(numRows,numCols,numPlots);

%The number of elements in the mask.
NMask=2*numGuardCells+2*numAvgCells+1;

%%Create the parts of the mask. The mask itself could be constructed using 
%Mask1=ones(NMask(:)');
%Mask2=zeros(NMask(:)');
%Mask2((numAvgCells(1)+1):(numAvgCells(1)+2*numGuardCells(1)+1),(numAvgCells(2)+1):(numAvgCells(2)+2*numGuardCells(2)+1))=1;
%Mask=Mask1-Mask2;
%However, we create it in parts to speed things up below.

Mask1R=ones(NMask(1),1);
Mask1D=ones(1,NMask(2));
%Mask1=Mask1R*Mask1D;
Mask2R=zeros(NMask(1),1);
Mask2R((numAvgCells(1)+1):(numAvgCells(1)+2*numGuardCells(1)+1),1)=1;
Mask2D=zeros(1,NMask(2));
Mask2D(1,(numAvgCells(2)+1):(numAvgCells(2)+2*numGuardCells(2)+1))=1;
%Mask2=Mask2R*Mask2D;
%and again, one could make the mask as
%Mask=Mask1-Mask2;

%The number of cells in the CFAR average.
NCFAR=prod(2*NFilter+1)-prod(2*numGuardCells+1); 

%This is the CFAR threshold from Equation 14 of [1].
T=PFA^(-1/NCFAR)-1;

%For each of the range-Doppler plots that was passed.
for curPlot=1:numPlots
    if(numPlots==1)
        plotMag2=abs(data2D).^2; %Somewhat faster without the indexing
    else
        plotMag2=abs(data2D(:,:,curPlot)).^2;
    end
    
    if(method==0)%Do everything at once by padding plotMag2        
        plotMag2=addAliasedPadding(plotMag2,NFilter);

        NseEst=conv2(T*Mask1R,Mask1D,plotMag2,'valid')-conv2(T*Mask2R,Mask2D,plotMag2,'valid');

        %[indexR,indexD] = find(plotMag2((NFilter(1)+1):(numRows+NFilter(1)),(NFilter(2)+1):(numCols+NFilter(2))) > NseEst);
        
        %These two lines are faster than the above line.
        data2DSqrScaled(:,:,curPlot)=plotMag2((NFilter(1)+1):(numRows+NFilter(1)),(NFilter(2)+1):(numCols+NFilter(2)))./NseEst;
    elseif(method==1||method==2)%Compute in segments
        %Do central region
        NseEst_cent=conv2(T*Mask1R,Mask1D,plotMag2,'valid')-conv2(T*Mask2R,Mask2D,plotMag2,'valid');
        
        RDRNumRows=4*NFilter(1);
        RDRNumCols=numCols+2*NFilter(2);
        
        %Extract the outer edges in range and periodically wrap.
        RD2_edgeR=zeros(RDRNumRows,RDRNumCols);
        RD2_edgeR(:,(NFilter(2)+1):(NFilter(2)+numCols))=[plotMag2((numRows-2*NFilter(1)+1):numRows,:);
                                                           plotMag2(1:2*NFilter(1),:)];                             
        RD2_edgeR(:,1:NFilter(2))=RD2_edgeR(:,(NFilter(2)+numCols-NFilter(2)+1):(NFilter(2)+numCols));
        RD2_edgeR(:,(numCols+NFilter(2)+1):RDRNumCols)=RD2_edgeR(:,(NFilter(2)+1):(2*NFilter(2)));

        NseEst_edgeR=conv2(T*Mask1R,Mask1D,RD2_edgeR,'valid')-conv2(T*Mask2R,Mask2D,RD2_edgeR,'valid');
        
        %The outer edges in Doppler
        RD2_edgeD=[plotMag2(:,(numCols-2*NFilter(2)+1):numCols), plotMag2(:,1:(2*NFilter(2)))];
        NseEst_edgeD=conv2(T*Mask1R,Mask1D,RD2_edgeD,'valid')-conv2(T*Mask2R,Mask2D,RD2_edgeD,'valid');
        if(method==1)
            %Combine into a single map: slow
            NseEst=[NseEst_edgeR((NFilter(1)+1):end,:);...
                    NseEst_edgeD(:,(NFilter(2)+1):end), NseEst_cent, NseEst_edgeD(:,1:NFilter(2));...
                    NseEst_edgeR(1:NFilter(1),:)];
            data2DSqrScaled(:,:,curPlot)=plotMag2./NseEst;
        else%Method 2
            %Scale section by section
            plotMag2(1:NFilter(1),:)=plotMag2(1:NFilter(1),:)./NseEst_edgeR((NFilter(1)+1):end,:);
            plotMag2((NFilter(1)+1):(numRows-NFilter(1)),1:NFilter(2))=plotMag2((NFilter(1)+1):(numRows-NFilter(1)),1:NFilter(2))./NseEst_edgeD(:,(NFilter(2)+1):end);
            plotMag2((NFilter(1)+1):numRows-NFilter(1),(NFilter(2)+1):(numCols-NFilter(2)))=plotMag2((NFilter(1)+1):numRows-NFilter(1),(NFilter(2)+1):(numCols-NFilter(2)))./NseEst_cent;
            plotMag2((NFilter(1)+1):numRows-NFilter(1),(numCols-NFilter(2)+1):numCols)=plotMag2((NFilter(1)+1):numRows-NFilter(1),(numCols-NFilter(2)+1):numCols)./NseEst_edgeD(:,1:NFilter(2));
            plotMag2((numRows-NFilter(1)+1):numRows,:)=plotMag2((numRows-NFilter(1)+1):numRows,:)./NseEst_edgeR(1:NFilter(1),:);
            data2DSqrScaled(:,:,curPlot)=plotMag2;
        end
    else
        error('Unknown method specified.')
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
