function centIdxVals=windowedGridCentroiding2D(index2D,weights,centroidWinLen)
%%WINDOWEDGRIDCENTROIDING2D Given a set of detections in two dimensions on
%                       an integer grid, cluster the detections
%                       together using a sliding window method described
%                       below.
%
%INPUTS: index2D A 2XN set of the index pairs of each of the detections.
%               The indices can be positive and negative, but must be
%               integers.
%       weights An NX1 or 1XN vector of all-positive weights associated
%               with each of the detections.
% centroidWinLen A 2X1 vector of the length of the window to use for
%               centroiding in each dimension. The length is given in
%               pixels. If a scalar is passed, then the same length is used
%               in both directions.
%
%OUTPUTS: centIdxVals A 2XnumCent matrix of the centroided index values of
%              the detections. These will generally not be integers.
%
%This function can be useful for centroiding detections from a range-
%Doppler plot.
%
%The algorithm creates a 2D mask of ondes of size 2*centroidWinLen(:).'+1.
%The mask is then convoluted with a sparse array containing all of the
%weights in the detection locations, as well as arrays containing the
%weights times the first and second indices. This means that each point is
%summed up with its neighbors in all directions at a distance determined by
%centroidWinLen. The nonzero elements of the convolution from the matrices
%times the row/ columns are divided by the one with just the weights. Thus,
%each point of the result is a weighted average of its nearest neighbors.
%Weighted averaged values that are less than 0.5 pixels from a nonzero
%convolution point in both directions are declared detections.
%
%EXAMPLE:
%Here we make a range-Doppler plot, detect three targets and then centroid
%the detections using this function.
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
% %Windowing in range and Doppler to lower sidelobes.
% Nx=length(x);
% wRange=windowFunSym(Nx,'Blackman',0);
% x=wRange.*x;
% wDoppler=windowFunSym(NB,'Nuttall',1);
% 
% [delayDopPlot,Doppler,delay]=delayDopplerPlotNBPulseDop(x,y,1,1,wDoppler,T0);
% 
% numGuardCells=[2;4];
% numAvgCells=[5;3];
% PFA=1e-7;
% 
% range=c*delay;
% rangeRate=Doppler*(c/fc);
% 
% DetectionList=CACFAR(delayDopPlot,numGuardCells,numAvgCells,PFA,[],0);
% idx=DetectionList(1).Index;
% 
% rVals=range(idx(1,:));
% RRVals=rangeRate(idx(2,:));
% 
% weights=abs(DetectionList(1).Value);
% centroidWinLen=4;
% centIdxVals=windowedGridCentroiding2D(idx,weights,centroidWinLen);
% 
% %Turn the centroided index values into range and range rate values.
% deltaR=range(2)-range(1);
% rValsCent=(centIdxVals(1,:)-1)*deltaR;
% 
% deltaRR=rangeRate(2)-rangeRate(1);
% minVal=rangeRate(1);
% RRValsCent=(centIdxVals(2,:)-1)*deltaRR+minVal;
% 
% %Display the CFAR detections and the centroided detections.
% figure(1)
% clf
% hold on
% scatter(RRVals,rVals,'ob')
% scatter(RRValsCent,rValsCent,'xr','linewidth',2)
% h1=xlabel('Range Rate (m/s)');
% h2=ylabel('Range (km)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% axis([rangeRate(1), rangeRate(end), range(1), range(end)])
%
%January 2017 Thomas Higgins, modified extensively by Daniel Scholnick and
%David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

indexR=index2D(1,:).';
indexD=index2D(2,:).';
weights=weights(:);

%This shifts the indices to be all positive and greater than zero.
minR=min(indexR);
minD=min(indexD);
indexR=indexR-minR+1;
indexD=indexD-minD+1;

%The dimensions of the sparse matrices below are minimized.
numRows=max(indexR);
numCols=max(indexD);

if(length(centroidWinLen)==1)
   centroidWinLen=[centroidWinLen;centroidWinLen];
end

%The matrix is sparse as it only contains detections.
Detection=sparse(indexR,indexD,weights,numRows,numCols);
%Weights multiplied by range bin indices
Detection_R=sparse(indexR,indexD,weights.*indexR,numRows,numCols); 
%Weights multiplied by Doppler bin indices
Detection_D=sparse(indexR,indexD,weights.*indexD,numRows,numCols); 

CM_mask=ones(2*centroidWinLen(:).'+1);
%This filtering does an implicit zero-padding, which means that a target
%that wraps around the edges of the viewing areas will potentially be
%split.
CM_sum=conv2Sparse(Detection,CM_mask,'same');
CM_R=conv2Sparse(Detection_R,CM_mask,'same');
CM_D=conv2Sparse(Detection_D,CM_mask,'same');

if(isa(CM_sum,'gpuArray'))
    %One can't use subsref on a gpuArray, so grab the nonzero
    %values directly.
    non_zero=find(CM_sum);
    nz_R=find(CM_R);
    nz_D=find(CM_D);
    if((length(non_zero)==length(nz_R)) && (length(non_zero)==length(nz_D)) && all(non_zero==nz_R) && all(non_zero==nz_D))
        CM_sum=nonzeros(CM_sum);
        CM_R=nonzeros(CM_R)./CM_sum;
        CM_D=nonzeros(CM_D)./CM_sum;
    else %Somehow we have different nonzeros for sum/R/D; this should be
         %very rare, but we cannot rule it out.
        error('Inconsistent Clustering values encountered. This dataset cannot be processed on the GPU.');
    end
else
    non_zero=find(CM_sum);
    CM_sum=full(CM_sum(non_zero));
    CM_R=full(CM_R(non_zero))./CM_sum;
    CM_D=full(CM_D(non_zero))./CM_sum;
end
[R, D]=ind2sub([numRows, numCols],non_zero);
R_diff=CM_R-R;
D_diff=CM_D-D;
centers=find((abs(R_diff)<=0.5)&(abs(D_diff)<=0.5));

%Undo the shifts in the indices that were applied at the beginning.
Rc=CM_R(centers)+minR-1;
Dc=CM_D(centers)+minD-1;

centIdxVals=[Rc.';Dc.'];
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
