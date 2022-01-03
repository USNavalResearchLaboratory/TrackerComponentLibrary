function tLocEst=polyMeasConvert(z,measTypes,sensorIdxLists,sensorStates,c,algorithm,opts,AbsTol,RelTol,scratchFolderPath,execPath)
%%POLYMEASCONVERT Given simultaneous measurements consisting of time-delay
%           of arrival, bistatic range, range rate (Doppler) from an
%           assumed stationary emitter (given moving sensors), and/ or the
%           ratio of the received frequency of measurements across sensors
%           (assuming a stationary target and moving sensors), determine
%           the location of the target. Unlike the function
%           polyMeasConvertAsync, this function only estimates a position
%           and thus cannot handle non-simultaneous measurements. This
%           function calls solvePolySysWithExtProg to solve a system of
%           equations. Some solvers can fail; Bertini appears to be the
%           most reliable.
%
%INPUTS: z A numDimX1 vector of the measurements. Note that each
%          measurement is a scalar and there must be exactly as many
%          measurements as the dimensionality of the problem being solved.
%          Thus, numDim=2 for 2D or numDim=3 for 3D.
%  measTypes A numDimX1 or 1XnumDim list specifying which type each of the
%          numDim measurements is. Possible values for each of the numDim
%          entries are
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
%          does not matter). For range rate measurements from an emitter,
%          only one receiver is used, so only sensorIdxLists(1,i) is used;
%          the second value can be set to zero. For frequency ratio
%          measurements, both components are used and sensorIdxLists(1,i)
%          is the numerator.
% sensorStates A numDimXnumSensors or (2*numDim)XnumSensors list of the
%          Cartesian sensor locations associated with the measurements and,
%          if frequency ratio or Doppler measurements are used, the
%          velocities of the sensors (where the 2*numDim dimensionality are
%          required and position components all come before velocity). For
%          sensors only used for TDOA or bistatic range measurements, any
%          provided velocity components are ignored and can be set to zero.
%          These sensors are selected via the sensorIdxLists input.
%        c The speed of signal propagation. This is needed when using TDOA
%          or frequency ratio measurements. If this parameter is omitted or
%          an empty matrix is passed, then c=299792458 (the speed of light
%          in a vacuum in meters per second) is used.
% algorithm An optional parameter specifying which algorithm to use. The
%          solvers for simultaneous multivariate polynomials are external
%          programs called via solvePolySysWithExtProg. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Bertini.
%          1 Use PHCpack.
%          2 Use the certified homotopy algorithm that is built into
%            Macaulay2 (in NAG4M2). This only uses the normalized solver
%            with the default options.
%     opts An optional input specifying options for the solver in the
%          function solvePolySysWithExtProg. This is described in more
%          detail in solvePolySysWithExtProg. Omitting this parameter or
%          passing am empty matrices uses the default values, except if
%          Bertini is used, then SecurityMaxNorm, EndpointFiniteThreshold
%          and PathTruncationThreshold are set to 1e9 by default and
%          RatioTolerance and FinalTol are respectively set to 0.01 and
%          1e-14 by default. If Bertini has been compiled with support for
%          the message passing interface and the mpirun command exists on
%          this computer, then one can include the parameter MPIRunProcs
%          set to an integer above 0 and that many processes will be used 
%          to parallelize Bertini. Setting the number too high will result
%          in slower performance. This parameter can be omitted or an empty
%          matrix passed if default values should not be changed.
% AbsTol, RelTol Absolute and relative tolerances for discarding false
%          solutions. Due to the squaring of equations, false solutions can
%          arise. Such solutions will produce, for example, TDOA
%          measurments with a sign opposite that of the input. ALso,
%          complex solutions are false, but one might want to keep complex
%          solutions that are very close to real solutions. Thus, the
%          algorithm tests real solutions and the complex solutions
%          (discarding the imaginary part) and converts them into the
%          measurement domain. If absDiff<AbsTol or
%          absDiff<RelTol*abs(z(curMeas)), where absDiff is the abslute
%          difference of the converted component versus the actual
%          component, then the measurement is kept. Default values are
%          respectively 1e-11 and 1e-8.
%  scratchFolderPath An optional parameter specifying the folder in which
%           temporary files will be written. This is needed, because
%           all solvers only works through files. If this parameter is
%           omitted or an empty matrix is passed, then a folder named temp
%           in the folder enclosing this function is used. Note that if an
%           error occurs while executing this function, then temporary
%           files might be left in the temp folder after this function
%           ends.
%  execPath The command line command to use to execute the solver. The 
%           default if this parameter is omitted or an empty matrix is 
%           passed is just the standard name on the command line: bertini,
%           phc and M2 for each of the algorithms. The default assumes that
%           the program is already in the default search path. For example,
%           on *NIX systems, one can usually put the executable in
%           /usr/local/bin for it to be in the default search path.
%           However, for some reason, that directory is often not
%           included in Matlab's search path. Matlab 2016a's help under
%           "Run External Commands, Scripts, and Programs" specifically
%           says how to add that folder to the default search path.
%
%OUTPUTS: tLocEst A numDimXnumEst matrix of the solutions to the problem.
%
%All of the problems addressed by this function are types of simultaneous
%multivariate polynomials, as discussed in [1]. The polynomials are solved
%using the specified algorithm.
%
%Note that if the solver fails, then this function will usually have an
%error when trying to read the solutions and the temporary files will not
%be deleted. Note that one can place calls to this function in
%try-catch statements to keep it from having an error if the solver fails.
%
%Bertini can solve all of the examples below. PHCpack will fail on Example
%9, causing an error and PHCpack produces inaccurate results on example 8.
%The Macaulay2 solver is less reliable.
%
%EXAMPLE 1:
%The first example is TDOA-only estimation in 2D.
% measTypes=[0;0];%Both TDOA
% sensorStates=zeros(2,3);
% sensorStates(:,1)=[9e3;39e3];
% sensorStates(:,2)=[65e3;10e3];
% sensorStates(:,3)=[64e3;71e3];
% tTrue=[27e3;42e3];%The true target location.
% c=299792458;%m/s The speed of light.
% 
% TDOA=zeros(2,1);
% sensorIdxLists=zeros(2,2);
% sensorIdxLists(:,1)=[1;2];%First TDOA pair.
% TDOA(1)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,2)))/c;
% sensorIdxLists(:,2)=[1;3];%Second TDOA pair.
% TDOA(2)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,3)))/c;
% tLocEst=polyMeasConvert(TDOA,measTypes,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 2:
%The second example is TDOA-only estimation in 3D.
% measTypes=[0;0;0];%All three TDOA
% sensorStates=zeros(3,4);
% sensorStates(:,1)=[9e3;9e3;100];
% sensorStates(:,2)=[65e3;-10e3;1e3];
% sensorStates(:,3)=[64e3;-30e3;8e3];
% sensorStates(:,4)=[-12e3;6e3;12];
% tTrue=[27e3;42e3;5e3];%The true target location.
% c=299792458;%m/s The speed of light.
% 
% TDOA=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%First TDOA pair.
% TDOA(1)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,2)))/c;
% sensorIdxLists(:,2)=[1;3];%Second TDOA pair.
% TDOA(2)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,3)))/c;
% sensorIdxLists(:,3)=[4;3];%Third TDOA pair.
% TDOA(3)=(norm(tTrue-sensorStates(:,4))-norm(tTrue-sensorStates(:,3)))/c;
% tLocEst=polyMeasConvert(TDOA,measTypes,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 3:
%The third example is range-only estimation in 2D
% measTypes=[1;1];%Both range
% sensorStates=zeros(2,3);
% sensorStates(:,1)=[9e3;39e3];
% sensorStates(:,2)=[65e3;10e3];
% sensorStates(:,3)=[64e3;71e3];
% tTrue=[27e3;42e3];%The true target location.
% 
% r=zeros(2,1);
% sensorIdxLists=zeros(2,2);
% sensorIdxLists(:,1)=[1;2];%First range pair.
% r(1)=norm(tTrue-sensorStates(:,1))+norm(tTrue-sensorStates(:,2));
% sensorIdxLists(:,2)=[1;3];%Second range pair.
% r(2)=norm(tTrue-sensorStates(:,1))+norm(tTrue-sensorStates(:,3));
% tLocEst=polyMeasConvert(r,measTypes,sensorIdxLists,sensorStates)
%
%EXAMPLE 4:
%The fourth example is two range and one TDOA measurement in 3D.
% measTypes=[1;1;0];%Two range, one TDOA
% sensorStates=zeros(3,4);
% sensorStates(:,1)=[9e3;9e3;100];
% sensorStates(:,2)=[65e3;-10e3;1e3];
% sensorStates(:,3)=[64e3;-30e3;8e3];
% sensorStates(:,4)=[-12e3;6e3;12];
% tTrue=[27e3;42e3;5e3];%The true target location.
% c=331.45;%m/s Approximate speed of sound in air.
% 
% z=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%First range pair.
% z(1)=norm(tTrue-sensorStates(:,1))+norm(tTrue-sensorStates(:,2));
% sensorIdxLists(:,2)=[1;3];%Second range pair.
% z(2)=norm(tTrue-sensorStates(:,1))+norm(tTrue-sensorStates(:,3));
% sensorIdxLists(:,3)=[4;3];%Third range pair.
% z(3)=(norm(tTrue-sensorStates(:,4))-norm(tTrue-sensorStates(:,3)))/c;
% tLocEst=polyMeasConvert(z,measTypes,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 5:
%The fifth-example is range rate-only (Doppler) estimation of an emitter in
%2D.
% measTypes=[2;2];
% sensorStates=zeros(4,2);
% sensorStates(1:2,:)=[500,1100;
%                      2500, 2500];%Sensor locations.
% sensorStates(3:4,:)=[300,300;
%                      0,  100];%Sensor velocities.
% tTrue=[1e3;5e3];%Stationary emitter location.
% 
% rr=zeros(2,1);
% sensorIdxLists=zeros(2,2);
% sensorIdxLists(1,1)=1;%Generate measurement from sensor 1.
% rr(1)=-sensorStates(3:4,1)'*(tTrue-sensorStates(1:2,1))/norm(tTrue-sensorStates(1:2,1));
% sensorIdxLists(1,2)=2;%Generate measurement from sensor 2.
% rr(2)=-sensorStates(3:4,2)'*(tTrue-sensorStates(1:2,2))/norm(tTrue-sensorStates(1:2,2));
% tLocEst=polyMeasConvert(rr,measTypes,sensorIdxLists,sensorStates)
%
%EXAMPLE 6:
%The sixth example is range rate-only (Doppler) estimation of an emitter in
%3D.
% measTypes=[2;2;2];
% sensorStates=zeros(6,3);
% sensorStates(1:3,:)=[500, 1100,5000;
%                     3000, 1000,-1000;
%                     0,    2000,1000];%Sensor locations.
% sensorStates(4:6,:)=[300, 100,0;
%                      0,   100,-300;
%                      0,   0,  0];%Sensor velocities.
% tTrue=[1e3;5e3;2e3];%Stationary emitter location.
% 
% rr=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(1,1)=1;%Generate measurement from sensor 1.
% rr(1)=-sensorStates(4:6,1)'*(tTrue-sensorStates(1:3,1))/norm(tTrue-sensorStates(1:3,1));
% sensorIdxLists(1,2)=2;%Generate measurement from sensor 2.
% rr(2)=-sensorStates(4:6,2)'*(tTrue-sensorStates(1:3,2))/norm(tTrue-sensorStates(1:3,2));
% sensorIdxLists(1,3)=3;%Generate measurement from sensor 3.
% rr(3)=-sensorStates(4:6,3)'*(tTrue-sensorStates(1:3,3))/norm(tTrue-sensorStates(1:3,3));
% tLocEst=polyMeasConvert(rr,measTypes,sensorIdxLists,sensorStates)
%
%EXAMPLE 7:
%The seventh example is frequency ratio-only estimation in 2D
% measTypes=[3;3];
% sensorStates=zeros(4,3);
% sensorStates(1:2,:)=[1000, 500,  1100;
%                      3000, 2500, 2500];%Sensor locations
% sensorStates(3:4,:)=[150, 300, 300;
%                     -150, 0,   0];%Sensor velocities.
% tTrue=[1e3;5e3];%True emitter location.
% c=299792458;%m/s The speed of light.
% 
% %Compue all three range rates. The frequency ratio can be derived from the
% %range rates
% rr=zeros(3,1);
% for curRx=1:3
%     rr(curRx)=-sensorStates(3:4,curRx)'*(tTrue-sensorStates(1:2,curRx))/norm(tTrue-sensorStates(1:2,curRx));
% end
% 
% fRat=zeros(2,1);
% sensorIdxLists=zeros(2,2);
% sensorIdxLists(:,1)=[1;2];%Ratio of sensor 1 to 2.
% fRat(1)=(1-rr(1)/c)./(1-rr(2)/c);
% sensorIdxLists(:,2)=[2;3];%Ratio of sensor 2 to 3.
% fRat(2)=(1-rr(2)/c)./(1-rr(3)/c);
% tLocEst=polyMeasConvert(fRat,measTypes,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 8:
%The eighth example is frequency ratio-only estimation in 3D
% measTypes=[3;3;3];
% sensorStates=zeros(3,4);
% sensorStates(1:3,:)=[1000, 500, 1100,2500;
%                      3000, 2500, 2500, 0
%                         0,    0,    0,    400];%Sensor locations
% sensorStates(4:6,:)=[150, 300, 300, 0;
%                     -150, 0,   0,   200;
%                        0, 0,  -100, 0;];%Sensor velocities.
% tTrue=[1e3;5e3;4e3];%True emitter location.
% c=299792458;%m/s The speed of light.
% 
% %Compue all four range rates. The frequency ratio can be derived from the
% %range rates
% rr=zeros(4,1);
% for curRx=1:4
%     rr(curRx)=-sensorStates(4:6,curRx)'*(tTrue-sensorStates(1:3,curRx))/norm(tTrue-sensorStates(1:3,curRx));
% end
% 
% fRat=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%Ratio of sensor 1 to 2.
% fRat(1)=(1-rr(1)/c)./(1-rr(2)/c);
% sensorIdxLists(:,2)=[1;3];%Ratio of sensor 1 to 3.
% fRat(2)=(1-rr(1)/c)./(1-rr(3)/c);
% sensorIdxLists(:,3)=[1;4];%Ratio of sensor 1 to 4.
% fRat(3)=(1-rr(1)/c)./(1-rr(4)/c);
% tLocEst=polyMeasConvert(fRat,measTypes,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 9:
%The ninth example is two TDOA and one frequency ratio measurement.
% measTypes=[0;3;0];
% sensorStates=zeros(3,4);
% sensorStates(1:3,:)=[1000, 500, 1100,2500;
%                      3000, 2500, 2500, 0
%                         0,    0,    0,    400];%Sensor locations
% sensorStates(4:6,:)=[150, 300, 300, 0;
%                     -150, 0,   0,   200;
%                        0, 0,  -100, 0;];%Sensor velocities.
% tTrue=[1e3;5e3;4e3];%True emitter location.
% c=299792458%331.45;%m/s Approximate speed of sound in air.
% 
% %Compue all four range rates. The frequency ratio can be derived from the
% %range rates
% rr=zeros(4,1);
% for curRx=1:4
%     rr(curRx)=-sensorStates(4:6,curRx)'*(tTrue-sensorStates(1:3,curRx))/norm(tTrue-sensorStates(1:3,curRx));
% end
% 
% z=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%TDOA between sensors 1 and 2.
% z(1)=(norm(tTrue-sensorStates(1:3,1))-norm(tTrue-sensorStates(1:3,2)))/c;
% sensorIdxLists(:,2)=[1;2];%Frequency ratio between sensors 1 and 2.
% z(2)=(1-rr(1)/c)./(1-rr(2)/c);
% sensorIdxLists(:,3)=[3;4];%TDOA between sensors 3 and 4.
% z(3)=(norm(tTrue-sensorStates(1:3,3))-norm(tTrue-sensorStates(1:3,4)))/c;
% tLocEst=polyMeasConvert(z,measTypes,sensorIdxLists,sensorStates,c,0)
% %Note that Bertini will correctly solve this system, but PHCpack will
% %fail.
%
%REFERENCES:
%[1] D. F. Crouse, "General Multivariate Polynomial Target Localization and
%    Initial Estimation," Journal of Advances in Information Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(c))
   c=299792458;%m/s The speed of light in a vacuum.
end

if(nargin<6||isempty(algorithm))
    algorithm=0;
end

if(nargin<7)
    opts=[];
end

if(nargin<8||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<9||isempty(RelTol))
    RelTol=1e-8;
end 

if(nargin<10)
    scratchFolderPath=[];
end

if(nargin<11)
    execPath=[];
end

if(algorithm==0&&isempty(opts))%If using Bertini, set default tolerances
    %Default options when using Bertini
    opts=struct('SecurityMaxNorm',1e9,'EndpointFiniteThreshold',1e9,'PathTruncationThreshold',1e9,'RatioTolerance',0.01,'FinalTol',1e-14);
elseif(algorithm==0)
    if(~isfield(opts,'SecurityMaxNorm'))
       opts.SecurityMaxNorm=1e9;
    end
    
    if(~isfield(opts,'EndpointFiniteThreshold'))
        opts.EndpointFiniteThreshold=1e9;
    end
    
    if(~isfield(opts,'PathTruncationThreshold'))
        opts.PathTruncationThreshold=1e9;
    end
    
    if(~isfield(opts,'RatioTolerance'))
        opts.RatioTolerance=0.01;
    end

    if(~isfield(opts,'FinalTol'))
        opts.FinalTol=1e-14;
    end
end

%The number of dimensions of the problem equals the number of measurements.
numDim=length(z);

%We subtract out the locations of the sensors and then put them back in at
%the end. This should help with numerical sensitivity issues when
%localizing things on the ground but using ECEF coordinates. The offset is
%removed in the end. We then also scale everything with regard to the most
%distant sensor.
centerOfRegion=mean(sensorStates(1:numDim,:),2);
sensorStates(1:numDim,:)=bsxfun(@minus,sensorStates(1:numDim,:),centerOfRegion);

scale=max(sqrt(sum(sensorStates.*sensorStates,1)));
sensorStates=sensorStates/scale;
sel=measTypes~=0&measTypes~=3;
z(sel)=z(sel)/scale;
c=c/scale;

numSensors=length(unique(sensorIdxLists(:)));

%These mark when an equation for a range from the target to a particular
%sensor has been added. When a range  measurement has been added, this
%gives the index of the variable used for it.
rangeEqAdded=zeros(numSensors,1);
numRangeMeasAdded=0;

%Allocate space for the maximum possible number of equations
xPolys=cell(9,1);

numPolysAdded=0;

switch(numDim)
    case 2
        xList=cell(2,1);
        rVars=cell(4,1);

        xList{1}='xa';
        xList{2}='xb';

        rVars{1}='ra';
        rVars{2}='rb';
        rVars{3}='rc';
        rVars{4}='rd';

        for curMeas=1:numDim
            switch(measTypes(curMeas))
                case 0%TDOA measurement
                    lRx1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lRx2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                    TDOA=z(curMeas);
                    
                    xPolyCur=TDOAMeas2DEq(TDOA,lRx1,lRx2,xList,c);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 1%Bistatic range measurement
                    lTx=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lRx=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                    r=z(curMeas);
                    
                    xPolyCur=rangeMeas2DEq(r,lTx,lRx,xList);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 2%Doppler measurement of an emitter
                    l=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lDot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));
                    rDot=z(curMeas);
                    
                    xPolyCur=emitterDopplerMeas2DEq(rDot,l,lDot,xList);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 3%Frequency ratio measurement.
                    l=cell(numDim,1);
                    lDot=cell(numDim,1);
                    selRVars=cell(numDim,1);
                    
                    %Get the appropriate sensor and range variables. If the
                    %range variables have not already been added to the
                    %estimation problem, then add them.
                    for curIdx=1:2
                        idx=sensorIdxLists(curIdx,curMeas);
                        l{curIdx}=sensorStates(1:numDim,idx);
                        lDot{curIdx}=sensorStates((numDim+1):(2*numDim),idx);
                    
                        %If no equation for the current range has been
                        %added.
                        if(rangeEqAdded(idx)==0)
                            numRangeMeasAdded=numRangeMeasAdded+1;
                            rangeEqAdded(idx)=numRangeMeasAdded;
                            selRVars{curIdx}=rVars{numRangeMeasAdded};

                            rPoly=create2DRangeEq(l{curIdx},xList,selRVars{curIdx});
                            numPolysAdded=numPolysAdded+1;
                            xPolys{numPolysAdded}=rPoly;
                        else
                            selRVars{curIdx}=rVars{rangeEqAdded(idx)};
                        end
                    end

                    xPolyCur=freqRatMeas2DEq(z(curMeas),l,lDot,xList,selRVars,c);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                otherwise
                    error('An invalid measurement type was specified')
            end
        end
    case 3
        xList=cell(3,1);
        rVars=cell(6,1);
        
        xList{1}='xa';
        xList{2}='xb';
        xList{3}='xc';

        rVars{1}='ra';
        rVars{2}='rb';
        rVars{3}='rc';
        rVars{4}='rd';
        rVars{5}='re';
        rVars{6}='rf';

        for curMeas=1:3
            switch(measTypes(curMeas))
                case 0%TDOA measurement
                    lRx1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lRx2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                    TDOA=z(curMeas);
                    
                    xPolyCur=TDOAMeas3DEq(TDOA,lRx1,lRx2,xList,c);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 1%Bistatic range measurement
                    lTx=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lRx=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                    r=z(curMeas);
                    
                    xPolyCur=rangeMeas3DEq(r,lTx,lRx,xList);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 2%Doppler measurement of an emitter
                    l=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                    lDot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));
                    rDot=z(curMeas);
                    
                    xPolyCur=emitterDopplerMeas3DEq(rDot,l,lDot,xList);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                case 3%Frequency ratio measurement.
                    l=cell(numDim,1);
                    lDot=cell(numDim,1);
                    selRVars=cell(numDim,1);
                    
                    %Get the appropriate sensor and range variables. If the
                    %range variables have not already been added to the
                    %estimation problem, then add them.
                    for curIdx=1:2
                        idx=sensorIdxLists(curIdx,curMeas);
                        l{curIdx}=sensorStates(1:numDim,idx);
                        lDot{curIdx}=sensorStates((numDim+1):(2*numDim),idx);
                    
                        %If no equation for the current range has been
                        %added.
                        if(rangeEqAdded(idx)==0)
                            numRangeMeasAdded=numRangeMeasAdded+1;
                            rangeEqAdded(idx)=numRangeMeasAdded;
                            selRVars{curIdx}=rVars{numRangeMeasAdded};

                            rPoly=create3DRangeEq(l{curIdx},xList,selRVars{curIdx});
                            numPolysAdded=numPolysAdded+1;
                            xPolys{numPolysAdded}=rPoly;
                        else
                            selRVars{curIdx}=rVars{rangeEqAdded(idx)};
                        end
                    end

                    xPolyCur=freqRatMeas3DEq(z(curMeas),l,lDot,xList,selRVars,c);
                    numPolysAdded=numPolysAdded+1;
                    xPolys{numPolysAdded}=xPolyCur;
                otherwise
                    error('An invalid measurement type was specified')
            end
        end
    otherwise
        error('Invalid number of measurements/ dimensionality')
end

%Shrink to fit.
xPolys=xPolys(1:numPolysAdded);

%The polynomials have been constructed; the optimization problem can now be
%solved.
varNames=[xList;rVars(rangeEqAdded(1:numRangeMeasAdded))];

zCart=real(solvePolySysWithExtProg(xPolys,varNames,algorithm,opts,scratchFolderPath,execPath));

%If no solutions were found, then the algorithm failed.
if(isempty(zCart))
    tLocEst=[];
    return; 
end

%Given the solutions, extraneous solutions arising only due to the squaring
%of some of the equations can be discarded. These include solutions that
%1) Have negative range estimates (for when using frequency ratio
%   measurements and unknown ranges become new variables/equations).
%2) When converted to TDOA or Doppler emitter the values have the opposite
%   signs of the original measurements.
%Additionally, complex solutions (for which we have discarded the complex
%part) are present. We do not completely discard complex solutions as one
%might want such solutions whose real parts are very close to valid
%solutions. "very close" is determined by RelTol and AbsTol.

%First, we remove any estimates with negative range estimates.
if(numRangeMeasAdded>0)
    sel=all(zCart((numDim+1):end,:)>=0,1);
    
    zCart=zCart(1:numDim,sel);
end

%Now, we remove solutions that are not close enough to the measurement
%values when transformed into the measurmeent domain.
numSol=size(zCart,2);
sel=true(numSol,1);
for curMeas=1:numDim
    switch(measTypes(curMeas))
        case 0%TDOA
            for curSol=1:numSol
                lRx1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                lRx2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
                
                TDOAComp=(norm(zCart(:,curSol)-lRx1)-norm(zCart(:,curSol)-lRx2))/c;
                absDiff=abs(TDOAComp-z(curMeas));
                
                if(~(absDiff<AbsTol || absDiff<RelTol*abs(z(curMeas))))
                    sel(curSol)=false;
                end
            end
        case 1%Range
            lRx1=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
            lRx2=sensorStates(1:numDim,sensorIdxLists(2,curMeas));
            
            for curSol=1:numSol
                rComp=norm(zCart(:,curSol)-lRx1)+norm(zCart(:,curSol)-lRx2);
                absDiff=abs(rComp-z(curMeas));
                
                if(~(absDiff<AbsTol || absDiff<RelTol*abs(z(curMeas))))
                    sel(curSol)=false;
                end
            end
        case 2%Emitter Doppler
            for curSol=1:numSol
                l=sensorStates(1:numDim,sensorIdxLists(1,curMeas));
                lDot=sensorStates((numDim+1):(2*numDim),sensorIdxLists(1,curMeas));
                
                diff=zCart(:,curSol)-l;
                rrComp=-lDot'*diff/norm(diff);

                absDiff=abs(rrComp-z(curMeas));

                if(~(absDiff<AbsTol || absDiff<RelTol*abs(z(curMeas))))
                    sel(curSol)=false;
                    break;
                end
            end
        otherwise%No tests for other measurement types.
    end
end

%Get rid of invalid solutions.
tLocEst=zCart(:,sel);

%Undo the scaling and add back in the offset that was present due to
%centering the coordinate system around the sensors.
tLocEst=tLocEst*scale;
tLocEst=bsxfun(@plus,tLocEst,centerOfRegion);
end

function xPoly=TDOAMeas2DEq(TDOA,lRx1,lRx2,xList,c)
    u=(lRx1+lRx2)/2;
    v=(lRx1-lRx2)/2;
    v1=v(1);
    v2=v(2);
    delta=c*TDOA/2;

    cTilde=delta^4-delta^2*norm(v)^2+(u'*v)^2-delta^2*norm(u)^2;
    lTilde=2*delta^2*u-2*(u'*v)*v;
    lt1=lTilde(1);
    lt2=lTilde(2);

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=lt2;
    coeffs(4)=2*v1*v2;
    coeffs(5)=(v1-delta)*(v1+delta);
    coeffs(6)=(v2-delta)*(v2+delta);

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));
    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{2},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{1},'^2+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2'];
end

function xPoly=TDOAMeas3DEq(TDOA,lRx1,lRx2,xList,c)
    u=(lRx1+lRx2)/2;
    v=(lRx1-lRx2)/2;

    v1=v(1);
    v2=v(2);
    v3=v(3);
    delta=c*TDOA/2;

    cTilde=delta^4-delta^2*norm(v)^2+(u'*v)^2-delta^2*norm(u)^2;
    lTilde=2*delta^2*u-2*(u'*v)*v;
    lt1=lTilde(1);
    lt2=lTilde(2);
    lt3=lTilde(3);

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=lt2;
    coeffs(4)=2*v1*v2;
    coeffs(5)=lt3;
    coeffs(6)=2*v1*v3;
    coeffs(7)=2*v2*v3;
    coeffs(8)=(v1-delta)*(v1+delta);
    coeffs(9)=(v2-delta)*(v2+delta);
    coeffs(10)=(v3-delta)*(v3+delta);

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{2},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{3},'+(',...
           num2str(coeffs(6),16),')*',xList{1},'*',xList{3},'+(',...
           num2str(coeffs(7),16),')*',xList{2},'*',xList{3},'+(',...
           num2str(coeffs(8),16),')*',xList{1},'^2+(',...
           num2str(coeffs(9),16),')*',xList{2},'^2+(',...
           num2str(coeffs(10),16),')*',xList{3},'^2'];
end

function xPoly=rangeMeas2DEq(r,lTx,lRx,xList)
    u=(lTx+lRx)/2;
    v=(lTx-lRx)/2;

    v1=v(1);
    v2=v(2);
    delta=r/2;

    cTilde=delta^4-delta^2*norm(v)^2+(u'*v)^2-delta^2*norm(u)^2;
    lTilde=2*delta^2*u-2*(u'*v)*v;
    lt1=lTilde(1);
    lt2=lTilde(2);

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=lt2;
    coeffs(4)=2*v1*v2;
    coeffs(5)=(v1-delta)*(v1+delta);
    coeffs(6)=(v2-delta)*(v2+delta);

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{2},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{1},'^2+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2'];
end

function xPoly=rangeMeas3DEq(r,lTx,lRx,xList)
    u=(lTx+lRx)/2;
    v=(lTx-lRx)/2;

    v1=v(1);
    v2=v(2);
    v3=v(3);
    delta=r/2;

    cTilde=delta^4-delta^2*norm(v)^2+(u'*v)^2-delta^2*norm(u)^2;
    lTilde=2*delta^2*u-2*(u'*v)*v;
    lt1=lTilde(1);
    lt2=lTilde(2);
    lt3=lTilde(3);

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=lt2;
    coeffs(4)=2*v1*v2;
    coeffs(5)=lt3;
    coeffs(6)=2*v1*v3;
    coeffs(7)=2*v2*v3;
    coeffs(8)=(v1-delta)*(v1+delta);
    coeffs(9)=(v2-delta)*(v2+delta);
    coeffs(10)=(v3-delta)*(v3+delta);

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{2},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{3},'+(',...
           num2str(coeffs(6),16),')*',xList{1},'*',xList{3},'+(',...
           num2str(coeffs(7),16),')*',xList{2},'*',xList{3},'+(',...
           num2str(coeffs(8),16),')*',xList{1},'^2+(',...
           num2str(coeffs(9),16),')*',xList{2},'^2+(',...
           num2str(coeffs(10),16),')*',xList{3},'^2'];
end

function xPoly=emitterDopplerMeas2DEq(rDot,l,lDot,xList)
    lDot1=lDot(1);
    lDot2=lDot(2);

    lTilde=2*lDot*(l'*lDot)-2*rDot^2*l;
    lt1=lTilde(1);
    lt2=lTilde(2);
    cTilde=rDot^2*norm(l)^2-(l'*lDot)^2;

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=rDot^2-lDot1^2;
    coeffs(4)=lt2;
    coeffs(5)=-2*lDot1*lDot2;
    coeffs(6)=rDot^2-lDot2^2;

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{1},'^2+(',...
           num2str(coeffs(4),16),')*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2'];
end

function xPoly=emitterDopplerMeas3DEq(rDot,l,lDot,xList)
    lDot1=lDot(1);
    lDot2=lDot(2);
    lDot3=lDot(3);

    lTilde=2*lDot*(l'*lDot)-2*rDot^2*l;
    lt1=lTilde(1);
    lt2=lTilde(2);
    lt3=lTilde(3);
    cTilde=rDot^2*norm(l)^2-(l'*lDot)^2;

    coeffs(1)=cTilde;
    coeffs(2)=lt1;
    coeffs(3)=rDot^2-lDot1^2;
    coeffs(4)=lt2;
    coeffs(5)=-2*lDot1*lDot2;
    coeffs(6)=rDot^2-lDot2^2;
    coeffs(7)=lt3;
    coeffs(8)=-2*lDot1*lDot3;
    coeffs(9)=-2*lDot2*lDot3;
    coeffs(10)=rDot^2-lDot3^2;

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',xList{1},'+(',...
           num2str(coeffs(3),16),')*',xList{1},'^2+(',...
           num2str(coeffs(4),16),')*',xList{2},'+(',...
           num2str(coeffs(5),16),')*',xList{1},'*',xList{2},'+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2+(',...
           num2str(coeffs(7),16),')*',xList{3},'+(',...
           num2str(coeffs(8),16),')*',xList{1},'*',xList{3},'+(',...
           num2str(coeffs(9),16),')*',xList{2},'*',xList{3},'+(',...
           num2str(coeffs(10),16),')*',xList{3},'^2'];
end

function rPoly=create2DRangeEq(lRx,xList,rVar)
    ljx=lRx(1);
    ljy=lRx(2);

    %The coefficients for an equation involving the r terms
    coeffs(1)=-ljx^2-ljy^2;%Constant term
    coeffs(2)=1;%rj^2
    coeffs(3)=2*ljx;%xa
    coeffs(4)=-1;%xa^2
    coeffs(5)=2*ljy;%xb
    coeffs(6)=-1;%xb^2

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    rPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',rVar,'^2+(',...
           num2str(coeffs(3),16),')*',xList{1},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'^2+(',...
           num2str(coeffs(5),16),')*',xList{2},'+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2'];
end

function rPoly=create3DRangeEq(lRx,xList,rVar)
    ljx=lRx(1);
    ljy=lRx(2);
    ljz=lRx(3);

    %The coefficients for an equation involving the r terms
    coeffs(1)=-ljx^2-ljy^2-ljz^2;%Constant term
    coeffs(2)=1;%rj^2
    coeffs(3)=2*ljx;%xa
    coeffs(4)=-1;%xa^2
    coeffs(5)=2*ljy;%xb
    coeffs(6)=-1;%xb^2
    coeffs(7)=2*ljz;%xc
    coeffs(8)=-1;%xc^2;

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    rPoly=[num2str(coeffs(1),16),'+(',...
           num2str(coeffs(2),16),')*',rVar,'^2+(',...
           num2str(coeffs(3),16),')*',xList{1},'+(',...
           num2str(coeffs(4),16),')*',xList{1},'^2+(',...
           num2str(coeffs(5),16),')*',xList{2},'+(',...
           num2str(coeffs(6),16),')*',xList{2},'^2+('...
           num2str(coeffs(7),16),')*',xList{3},'+(',...
           num2str(coeffs(8),16),')*',xList{3},'^2'];
end

function xPoly=freqRatMeas2DEq(fRat,l,lDot,xList,rList,c)
    lRx1=l{1};
    lRx2=l{2};
    lRxDot1=lDot{1};
    lRxDot2=lDot{2};

    l1x=lRx1(1);
    l1y=lRx1(2);
    l1Dotx=lRxDot1(1);
    l1Doty=lRxDot1(2);
    ljx=lRx2(1);
    ljy=lRx2(2);
    ljDotx=lRxDot2(1);
    ljDoty=lRxDot2(2);

    %The coefficients for the measurement equation.                
    coeffs(1)=-fRat*(ljDotx*ljx+ljDoty*ljy);%r1 term
    coeffs(2)=(l1Dotx*l1x+l1Doty*l1y);%rj
    coeffs(3)=c*(fRat-1);%r1*rj
    coeffs(4)=fRat*ljDotx;%r1*xa
    coeffs(5)=-l1Dotx;%rj*xa
    coeffs(6)=fRat*ljDoty;%r1*xb
    coeffs(7)=-l1Doty;%rj*xb

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=['(',num2str(coeffs(1),16),')*',rList{1},'+(',...
               num2str(coeffs(2),16),')*',rList{2},'+(',...
               num2str(coeffs(3),16),')*',rList{1},'*',rList{2},'+(',...
               num2str(coeffs(4),16),')*',rList{1},'*',xList{1},'+(',...
               num2str(coeffs(5),16),')*',rList{2},'*',xList{1},'+(',...
               num2str(coeffs(6),16),')*',rList{1},'*',xList{2},'+(',...
               num2str(coeffs(7),16),')*',rList{2},'*',xList{2}];
end

function xPoly=freqRatMeas3DEq(fRat,l,lDot,xList,rList,c)
    lRx1=l{1};
    lRx2=l{2};
    lRxDot1=lDot{1};
    lRxDot2=lDot{2};

    l1x=lRx1(1);
    l1y=lRx1(2);
    l1z=lRx1(3);
    l1Dotx=lRxDot1(1);
    l1Doty=lRxDot1(2);
    l1Dotz=lRxDot1(3);
    ljx=lRx2(1);
    ljy=lRx2(2);
    ljz=lRx2(3);
    ljDotx=lRxDot2(1);
    ljDoty=lRxDot2(2);
    ljDotz=lRxDot2(3);

    %The coefficients for the measurement equation.                
    coeffs(1)=-fRat*(ljDotx*ljx+ljDoty*ljy+ljDotz*ljz);%r1 term
    coeffs(2)=(l1Dotx*l1x+l1Doty*l1y+l1Dotz*l1z);%rj
    coeffs(3)=c*(fRat-1);%r1*rj
    coeffs(4)=fRat*ljDotx;%r1*xa
    coeffs(5)=-l1Dotx;%rj*xa
    coeffs(6)=fRat*ljDoty;%r1*xb
    coeffs(7)=-l1Doty;%rj*xb
    coeffs(8)=fRat*ljDotz;%r1*xc
    coeffs(9)=-l1Dotz;%rj*xc

    %Scale the equation.
    coeffs=coeffs/max(abs(coeffs));

    xPoly=['(',num2str(coeffs(1),16),')*',rList{1},'+(',...
           num2str(coeffs(2),16),')*',rList{2},'+(',...
           num2str(coeffs(3),16),')*',rList{1},'*',rList{2},'+(',...
           num2str(coeffs(4),16),')*',rList{1},'*',xList{1},'+(',...
           num2str(coeffs(5),16),')*',rList{2},'*',xList{1},'+(',...
           num2str(coeffs(6),16),')*',rList{1},'*',xList{2},'+(',...
           num2str(coeffs(7),16),')*',rList{2},'*',xList{2},'+(',...
           num2str(coeffs(8),16),')*',rList{1},'*',xList{3},'+(',...
           num2str(coeffs(9),16),')*',rList{2},'*',xList{3}];
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
