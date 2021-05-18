function DetectionList=OSCFAR(data2D,numGuardCells,numAvgCells,k,PFA)
%%OSCFAR Perform order-statistic constant false alarm rate (OS-CFAR)
%        detection on a two-dimensional grid. This can be, for example, a
%        range-Doppler plot. This function just performs the detection, it
%        does not centroid the detections.
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
%             k The integer order to use k. That is, the kth largest sample
%               in the test region is used as the test statistic.
%           PFA The scalar probability of false alarm between 0 and 1 that
%               determines the threshold for detection.
%
%OUTPUTS: DetectionList A numPlotsX1 collection of structures.
%               DetectionList(i).Index provides a 2XnumDetect set of the
%               row and column indices of each detection in the ith plot.
%               The values from data2D of the detections are given in
%               DetectionList(i).Value
%
%This implements the algorithm described in Section 4 of [1].
%
%EXAMPLE:
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
% imagesc([rangeRate(1),rangeRate(end)],[range(1), range(end)]/1e3,10*log10(abs(delayDopPlot)));
% set(gca,'YDir','normal')
% caxis([-30 27])
% colormap(jet(256))
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% k=165;
% DetectionList=OSCFAR(delayDopPlot,numGuardCells,numAvgCells,k,PFA);
% idx=DetectionList(1).Index;
% 
% rVals=range(idx(1,:));
% DopVals=rangeRate(idx(2,:));
% 
% %Display the CFAR detections.
% figure(2)
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
%REFERENCES
%[1] P. P. Gandhi and S. A. Kassam, "Analysis of CFAR processors in
%    nonhomogeneous background," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 24, no. 4, pp. 427-445, Jul. 1988.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(length(numGuardCells)==1)
    numGuardCells=[numGuardCells;numGuardCells];
end

if(length(numAvgCells)==1)
    numAvgCells=[numAvgCells;numAvgCells];
end

NFilter=numGuardCells+numAvgCells;

numRows=size(data2D,1);
numCols=size(data2D,2);
numEls=numRows*numCols;

numPlots=size(data2D,3);

%The number of elements in the mask.
NMask=2*numGuardCells+2*numAvgCells+1;

%Create the mask. The ones in the mask select which cells below the mask
%are considered  when performing the sorting and then averaging.
Mask1=ones(NMask(:)');
Mask2=zeros(NMask(:)');
Mask2((numAvgCells(1)+1):(numAvgCells(1)+2*numGuardCells(1)+1),(numAvgCells(2)+1):(numAvgCells(2)+2*numGuardCells(2)+1))=1;
Mask=logical(Mask1-Mask2);

%The number of cells in the mask
NCFAR=prod(2*NFilter+1)-prod(2*numGuardCells+1); 

%The threshold for detection.
T=OSCFARThreshold4PFA(PFA,NCFAR,k);

%Allocate space for the return variables.
DetectionList(numPlots).Index=[];
DetectionList(numPlots).Value=[];

%For each of the range-Doppler plots that was passed.
for curPlot=1:numPlots
    plotMag2=abs(data2D(:,:,curPlot)).^2;
    
    %Preallocate the maximum amount of space for detections. It can get
    %shrunk on return.
    indexR=zeros(numEls,1);
    indexD=zeros(numEls,1);
    numDet=0;
    
    %Extend
    plotMag2=addAliasedPadding(plotMag2,NFilter);
    
    %The initial span of the rows covered by the mask.
    rowSpan=1:NMask(1);
    for curRow=1:numRows
        %The initial span of the columns covered by the mask.
        colSpan=1:NMask(2);
        for curCol=1:numCols
            curVals=plotMag2(rowSpan,colSpan);
            sortedMaskedVals=sort(curVals(Mask),'ascend');
            
            if(sortedMaskedVals(k)*T<plotMag2(NFilter(1)+curRow,NFilter(2)+curCol))
                numDet=numDet+1;
                indexR(numDet)=curRow;
                indexD(numDet)=curCol;
            end
            colSpan=colSpan+1;
        end
        rowSpan=rowSpan+1;
    end
    
    %Shink to fit the actual detections.
    indexR=indexR(1:numDet);
    indexD=indexD(1:numDet);
    
    DetectionList(curPlot).Index=[indexR.';indexD.'];
    idx=sub2ind([numRows,numCols],indexR,indexD)+(curPlot-1)*numRows*numCols;
    
    DetectionList(curPlot).Value=data2D(idx);
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
