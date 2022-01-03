function FIM=computePolyMeasFIM(x,sigma2List,fTx,measTypes,sensorIdxLists,sensorStates,c,xi,w)
%%COMPUTEPOLYMEASFIM Given simultaneous measurements consisting of
%           time-delay of arrival, bistatic range, range rate (Doppler)
%           from an assumed stationary emitter (given moving sensors),
%           and/ or the received frequency of measurements from
%           multiple sensors (assuming a stationary target and moving
%           sensors), determine the inverse Cramér-Rao lower bound (CRLB)
%           (the Fisher information matrix) for an estimate of the location
%           of the target assuming measurements are corrupted with Gaussian
%           noise. Cubature integration is used for evaluation of the
%           necessary integral.
%
%INPUTS: x The numDimX1 true target location in Cartesian coordinates. This
%          will generally be 2X1 or 3X1.
% sigma2List A numMeasX1 or 1XnumMeas list of the variance of each 
%          measurement. The measurement errors are assumed independent. The
%          variance must be positive.
%      fTx If emitter frequency measurements are used, then this is the
%          true (non-Doppler-shifted) frequency of the emitter. Otherwise,
%          an empty matrix can be passed. It is assumed that all sensors
%          measuring a frequency measure the same one (i.e. not that one
%          sensor measurements emissions on one band and another emissions
%          on a disjoint band for something emitting on multiple bands).
% measTypes A numMeasX1 or 1XnumMeas list specifying which type each of the
%          numMeas measurements is. There can be more measurements than the
%          dimensionality of the state. Possible values for each of the
%          entries are:
%          0 TDOA measurement.
%          1 Bistatic range (collocate the transmitter and receiver for
%            monostatic).
%          2 Range rate measurements of an emitter (sensors must be
%            moving). This is proportional to the Doppler shift.
%          3 Frequency ratio measurements (sensors must be moving).
% sensorIdxLists A 2XnumMeas matrix of indices of which sensor locations
%          in the collection sensorLocs are associated with each
%          measurement. sensorIdxLists(:,i) concerns the ith measurement.
%          For TDOA measurements sensorIdxLists(1,i) is the reference
%          sensor from which the measured delay is taken with respect to
%          sensorIdxLists(2,i). For range measurements, the two indices
%          indicate which sensors are the transmitter and receiver (order
%          does not matter). For range-rate  and frequency measurements
%          from an emitter only sensorIdxLists(1,i) is used; the second
%          value can be set to zero.
% sensorStates A numDimXnumSensors or (2*numDim)XnumSensors list of the
%          Cartesian sensor locations associated with the measurements and,
%          if range rate or frequency measurements are used, the velocities
%          of the sensors (where the 2*numDim dimensionality are
%          required and position components all come before velocity). For
%          sensors only used for TDOA or bistatic range measurements, any
%          provided velocity components are ignored and can be set to zero.
%          These sensors are selected via the sensorIdxLists input.
%        c The speed of signal propagation. This is needed when using TDOA
%          or frequency measurements. If this parameter is omitted or an
%          empty matrix is passed, then c=299792458 (the speed of light in
%          a vacuum in meters per second) is used.
%    xi, w Cubature points are required for the expected value in the
%          computation of the CRLB (assuming Gaussian measurement noise).
%          xi and w are the cubature points and weights. The cubature
%          points must be numMeas-dimensional. If these parameters are
%          omitted, then the default set of points from the function
%          fifthOrderCubPoints is used.
%
%OUTPUTS: FIM If frequency measurements are not used, then this is a
%             numDimXnumDim Fisher information matrix (inverse CRLB matrix)
%             for the target position. If frequency measurements are used,
%             then this is a (numDim+1)X(numDim+1) Fisher information
%             matrix where the final row/ column element deals with the
%             estimate of fTx.
%
%This function implements the equations for the Fisher-information matrix/
%CRLB given in [1] assuming that measurements are corrupted with Gaussian
%noise. Cubature integration is used to solve the integral.
%
%EXAMPLE I: Two bistatic ranges and one TDOA measurement
%This example is given in [1]. the FIM is evaluated over a coarse grid.
% numPoints=50;
% measTypes=[0;0;1];
% c=299792458;%The speed of light in m/s.
% sigmaR=10;%10 meter range standard deviation.
% sigmaTau=10/c;%Standard deviation in TDOA, equivalent to 10m range.
% 
% sigma2List=[sigmaTau;sigmaTau;sigmaR].^2;
% %sigma2List=[sigmaTau;sigmaTau;].^2;
% 
% %Latitude/longitudes of the sensors
% sensorLatLon=[[20.265901;-155.857544],...
% [19.878939; -155.107727],...
% [19.661825; -156.091003],...
% [20.069960; -155.434570]]*(pi/180);
% sensorAlt=[8000,7500,6000,5000];
% 
% cartSensorLocs=ellips2Cart([sensorLatLon;sensorAlt]);
% 
% targetAlt=4205;%meters
% 
% lonMin=-156*(pi/180);
% lonMax=-154.75*(pi/180);
% latMin=18.5*(pi/180);
% latMax=19.75*(pi/180);
% 
% lonPoints=linspace(lonMin,lonMax,numPoints);
% latPoints=linspace(latMin,latMax,numPoints);
% 
% [lon,lat]=meshgrid(lonPoints,latPoints);
% numGridPoints=length(lat(:));
% 
% latLonAlt=[lat(:)';lon(:)';targetAlt*ones(1,numGridPoints)];
% cartTarLocs=ellips2Cart(latLonAlt);
% 
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%TDOA between sensors 1 and 2.
% sensorIdxLists(:,2)=[1;3];%TDOA between sensors 1 and 3.
% sensorIdxLists(:,3)=[3;4];%Bistatic range, sensors 3 and 4.
% 
% [xi,w]=fifthOrderCubPoints(3);
% %[xi,w]=seventhOrderCubPoints(3);
% %[xi,w]=arbOrderGaussCubPoints(3,10);%ninth-order
% RMSEMin=zeros(numGridPoints,1);
% for curGridPoint=1:numGridPoints
%     x=cartTarLocs(:,curGridPoint);
% 
%     FIM=computePolyMeasFIM(x,sigma2List,[],measTypes,sensorIdxLists,cartSensorLocs,c,xi,w);
%     RMSEMin(curGridPoint)=sqrt(trace(pinv(FIM)));
% end
% 
% markerSize=128*10;
% 
% figure()
% clf
% hold on
% set(gca,'YDir','normal')
% 
% %Plot the RMSE
% imagesc([lonMin lonMax]*(180/pi),[latMin latMax]*(180/pi),reshape(RMSEMin,numPoints,numPoints))
% caxis([0;2e4])
% colorbar()
% axis tight
% h1=xlabel('Longitude');
% h2=ylabel('Latitude');
% 
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%EXAMPLE II:
%This uses frequency-ratio-only measurements. This example is given in [1].
% numPoints=50;
% measTypes=[3;3;3;3];
% 
% c=299792458;%The speed of light in m/s.
% %c=331.45;%Approximate speed of sound in m/s.
% fTx=8e9;%8Ghz (X-Band) transmit frequency
% %fTx=24e3;%24kHz
% sigmaf=1;%1Hz frequency standard deviation.
% sigma2List=[sigmaf;sigmaf;sigmaf;sigmaf].^2;
% 
% %Latitude/longitudes of the sensors
% sensorLatLon=[[20.265901;-155.857544],...
% [19.878939; -155.107727],...
% [19.661825; -156.091003],...
% [20.069960; -155.434570]]*(pi/180);
% sensorAlt=[8000,7500,6000,5000];
% 
% cartSensorLocs=ellips2Cart([sensorLatLon;sensorAlt]);
% 
% %Determine sensor velocities.
% sensorVel=zeros(3,2);%Allocate space
% uENU=getENUAxes([sensorLatLon(:,1);sensorAlt(1)]);
% sensorVel(:,1)=(uENU(:,1)-uENU(:,2))/2*300;
% 
% uENU=getENUAxes([sensorLatLon(:,2);sensorAlt(2)]);
% sensorVel(:,2)=(uENU(:,1)+uENU(:,2))/2*300;
% 
% uENU=getENUAxes([sensorLatLon(:,3);sensorAlt(3)]);
% sensorVel(:,3)=uENU(:,1)*250;
% 
% uENU=getENUAxes([sensorLatLon(:,4);sensorAlt(4)]);
% sensorVel(:,4)=uENU(:,2)*250;
% 
% sensorStates=[cartSensorLocs;sensorVel];
% 
% targetAlt=4205;%meters
% 
% lonMin=-156*(pi/180);
% lonMax=-154.75*(pi/180);
% latMin=18.5*(pi/180);
% latMax=19.75*(pi/180);
% 
% lonPoints=linspace(lonMin,lonMax,numPoints);
% latPoints=linspace(latMin,latMax,numPoints);
% 
% [lon,lat]=meshgrid(lonPoints,latPoints);
% numGridPoints=length(lat(:));
% 
% latLonAlt=[lat(:)';lon(:)';targetAlt*ones(1,numGridPoints)];
% cartTarLocs=ellips2Cart(latLonAlt);
% 
% sensorIdxLists=zeros(2,4);
% sensorIdxLists(:,1)=[1;0];%Frequency.
% sensorIdxLists(:,2)=[2;0];%Frequency.
% sensorIdxLists(:,3)=[3;0];%Frequency.
% sensorIdxLists(:,4)=[4;0];%Frequency.
% 
% [xi,w]=fifthOrderCubPoints(4);
% RMSEMin=zeros(numGridPoints,1);
% for curGridPoint=1:numGridPoints
%     x=cartTarLocs(:,curGridPoint);
%     FIM=computePolyMeasFIM(x,sigma2List,fTx,measTypes,sensorIdxLists,sensorStates,c,xi,w);
%     CRLB=inv(FIM);
%     RMSEMin(curGridPoint)=sqrt(trace(CRLB(1:3,1:3)));
% end
% 
% markerSize=128*10;
% 
% figure()
% clf
% hold on
% set(gca,'YDir','normal')
% 
% %Plot the RMSE
% imagesc([lonMin lonMax]*(180/pi),[latMin latMax]*(180/pi),reshape(RMSEMin,numPoints,numPoints))
% caxis([0;5e3])
% colormap(parula(256))
% colorbar();
% axis tight
% 
% h1=xlabel('Longitude');
% h2=ylabel('Latitude');
% 
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] D. F. Crouse, "General Multivariate Polynomial Target Localization and
%    Initial Estimation," Journal of Advances in Information Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%April 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=length(x);
numMeas=length(measTypes);
hasfTx=any(measTypes==3);

if(nargin<8||isempty(xi))
    [xi,w]=fifthOrderCubPoints(numMeas);
end
numCubPoints=length(w);

if(nargin<7||isempty(c))
   c=299792458;
end

%We will compute the noise-free measurement values.
measValsTrue=zeros(numMeas,1);
for curMeas=1:numMeas
    l1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
    
    switch(measTypes(curMeas))
        case 0%TDOA measurement.
            l2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
            measValsTrue(curMeas)=(norm(x-l1)-norm(x-l2))/c;
        case 1%Bistatic range measurement.
            l2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
            measValsTrue(curMeas)=norm(x-l1)+norm(x-l2);
        case 2%Range rate measurement of an emitter.
            l1Dot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));
            measValsTrue(curMeas)=-(x-l1)'*l1Dot/norm(x-l1);
        case 3%Frequency measurement.
            l1Dot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));
            rr=-(x-l1)'*l1Dot/norm(x-l1);
            measValsTrue(curMeas)=(1-rr/c)*fTx;
        otherwise
            error('Unknown measurement type specified')
    end
end

%Transform the cubature points to match the given distribution.
xi=transformCubPoints(xi,measValsTrue,diag(sqrt(sigma2List)));

%Allocate space
FIM=zeros(numDim+hasfTx,numDim+hasfTx);
for curCub=1:numCubPoints
    gradVec=zeros(numDim+hasfTx,1);
    for curMeas=1:numMeas
        l1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));

        switch(measTypes(curMeas))
            case 0%TDOA measurement.
                l2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                sigmaTDOA2=sigma2List(curMeas);
                TDOA=xi(curMeas,curCub);

                norm1=norm(x-l1);
                norm2=norm(x-l2);

                gradVec(1:numDim)=gradVec(1:numDim)-(1/(c*sigmaTDOA2))*(TDOA-(1/c)*(norm1-norm2))*((l1-x)/norm1-(l2-x)/norm2);
            case 1%Bistatic range measurement.
                l2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                sigmaR2=sigma2List(curMeas);
                rB=xi(curMeas,curCub);

                norm1=norm(x-l1);
                norm2=norm(x-l2);

                gradVec(1:numDim)=gradVec(1:numDim)-(1/sigmaR2)*(rB-(norm1+norm2))*((l1-x)/norm1+(l2-x)/norm2);
            case 2%Range rate measurement of an emitter.
                sigmarDot2=sigma2List(curMeas);
                rDot=xi(curMeas,curCub);
                l1Dot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));

                norm1=norm(x-l1);

                gradVec(1:numDim)=gradVec(1:numDim)+(1/sigmarDot2)*(rDot+(x-l1)'*l1Dot/norm1)*(l1Dot/norm1+(x-l1)'*l1Dot/norm1^3*(l1-x));
            case 3%Frequency measurement.
                sigmaf2=sigma2List(curMeas);
                f=xi(curMeas,curCub);
                l1Dot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));

                norm1=norm(x-l1);
                fDiffTerm=f-(1+(1/c)*(x-l1)'*l1Dot/norm(x-l1))*fTx;

                %The position components paet of the gradient.
                gradVec(1:numDim)=gradVec(1:numDim)+fTx/(c*sigmaf2)*fDiffTerm*(l1Dot/norm1+(x-l1)'*l1Dot/norm1^3*(l1-x));
                %The frequency part of the gradient.
                gradVec(end)=gradVec(end)+(1/sigmaf2)*fDiffTerm*(1+(1/c)*(x-l1)'*l1Dot/norm1);
            otherwise
                error('Unknown measurement type specified')
        end
    end

    FIM=FIM+w(curCub)*(gradVec*gradVec');
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
