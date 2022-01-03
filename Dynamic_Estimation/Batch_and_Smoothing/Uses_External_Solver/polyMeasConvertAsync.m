function xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,c,HMats,RMats,algorithm,opts,AbsTol,RelTol,scratchFolderPath,execPath)
%%POLYMEASCONVERTASYNCH Given measurements consisting of time-delay of
%           arrival and/ or bistatic range, determine a state vector for a
%           target. Unlike the function polyMeasConvert, this function can
%           handle non-simultaneous measurements (and thus estimate
%           velocity or other components).
%
%INPUTS: z An xDimX1 array of the stacked measurements. Note that the
%          combined dimensionality of the measurements must equal the
%          dimensionality of the state being estimated (redundancies
%          cannot be used). There cannot be more than 52 variables.
%          However, the solvers will be too slow to use long before the
%          number of variables is that high.
%  measTypes An xDimX1 or 1XxDim list specifying which type each of the
%          numDim measurements is. Possible values for each of the xDim
%          entries are
%          0 Time delay of arrival (TDOA) measurement.
%          1 Bistatic range (collocate the transmitter and receiver for
%            monostatic).
%          2 A direction of arrival (DOA) measurement in the form of
%            direction cosines in three dimensions (meaning a u-v
%            measurement). When this measurement type is used, it is
%            assumed that the position components of the state are 3D.
%   FhMats A posDimXxDimXxDim set of xDim combined matrices H*F, one for
%          each measurement. That propagate the target state from the time
%          where the estimate is desired to the time of the measurement
%          (the F part) and then extract the position components (the H
%          part. For measurements at the time of the state, the F part
%          would be the identity matrix, and  if position components come
%          before velocity, in 3D, H=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0];
%          In this case, Fh=H for that time.
% sensorIdxLists A 2XxDim matrix of indices of which sensor locations
%          in the collection sensorLocs are associated with each
%          measurement. sensorIdxLists(:,i) concerns the ith measurement.
%          For TDOA measurements sensorIdxLists(1,i) is the reference
%          sensor from which the measured delay is taken with respect to
%          sensorIdxLists(2,i). For range measurements, the two indices
%          indicate which sensors are the transmitter and receiver (order
%          does not matter).
% sensorStates An xDimXnumSensors list of the Cartesian sensor locations
%          associated with the measurements. These sensors are selected via
%          the sensorIdxLists input.
%        c The speed of signal propagation. This is needed when using TDOA
%          measurements and is otherwise ignored. if this parameter is
%          omitted or an empty matrix is passed, then c=299792458 (the
%          speed of light in a vacuum in meters per second) is used.
% HMats,RMats These two cell arrays are only used if DOA measurements are
%          provided. DOAs measurements are assumed to be u-v measurement in
%          the local coordinate system of the recever. The matrix RMats{i}
%          is such that RMats{i}*v rotates vector
%          v=FhMats(:,:,i)*x-sensorStates(:,sensorIdxLists(1,i)
%          (vector from the sensor to the target) into the local coordinate
%          system of the receiver. Matrix HMats{i} then selects the two
%          components representing u and v in the local coordinate system.
%          Such matrices only need to be given for DOA measurements. The
%          third component in the local coordinate system fo the receiver
%          is assumed to be positive (the target is in front of the
%          receiver) and solutions with negative values are discarded.
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
%          FinalTol is set to 1e-12 by default.
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
%          respectively 1e-11 and 1e-8. If Bertini has been compiled with
%          support for the message passing interface and the mpirun command
%          exists on
%          this computer, then one can include the parameter MPIRunProcs set
%          to an integer above 0 and that many processes will be used to
%          parallelize Bertini. Setting the number too high will result in
%          slower performance. This parameter can be omitted or an empty
%          matrix passed if default values should not be changed.
% scratchFolderPath An optional parameter specifying the folder in which
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
%OUTPUTS: xEst An xEstXnumEst matrix of the solutions to the problem.
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
%Bertini can solve all of the examples below. The other solvers are not
%necessarily as reliable.
%
%EXAMPLE 1:
%A TDOA-only example in 2D.
%We have three stationary sensors forming TDOA measurements at four times.
% sensorStates=zeros(2,2);
% sensorStates(:,1)=[9e3;39e3];
% sensorStates(:,2)=[65e3;10e3];
% sensorStates(:,3)=[10e3;65e3];
% measTypes=zeros(4,1);
% 
% xTrue=[27e3;42e3;300;0];%The true target location and velocity
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(4,4);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,4,1);%State transition matrix from time 0 to 1.
% F2=FPolyKal(2*T,4,1);%State transition matrix to time 2.
% F3=FPolyKal(3*T,4,1);%State transition matrix to time 3.
% H=[1, 0, 0, 0;
%    0, 1, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(2,4,4);%There are four "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F2;
% FhMats(:,:,4)=H*F3;
% c=299792458;%m/s The speed of light.
% 
% %For the demonstration, we will use noise-free TDOA measurements.
% z=zeros(4,1);
% sensorIdxLists=[1,1,1,1;
%                 2,2,3,3];
% for i=1:2
%     z(i)=(norm(FhMats(:,:,i)*xTrue-sensorStates(:,1))-norm(FhMats(:,:,i)*xTrue-sensorStates(:,2)))/c;
% end
% for i=3:4
%     z(i)=(norm(FhMats(:,:,i)*xTrue-sensorStates(:,1))-norm(FhMats(:,:,i)*xTrue-sensorStates(:,3)))/c;
% end
% 
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 2:
%Here we have TDOA-only estimation in 3D. We have three sensors observing
%the target at four times.
% sensorStates=zeros(3,3);
% sensorStates(:,1)=[9e3;9e3;100];
% sensorStates(:,2)=[65e3;-10e3;1e3];
% sensorStates(:,3)=[64e3;-30e3;8e3];
% measTypes=zeros(6,1);
% 
% xTrue=[27e3;42e3;5e3;300;0;0];%The true target location and velocity
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(6,6);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,6,1);%State transition matrix from time 0 to 1.
% F3=FPolyKal(3*T,6,1);%State transition matrix to time 3.
% F4=FPolyKal(4*T,6,1);%State transition matrix to time 4.
% H=[1, 0, 0, 0, 0, 0;
%    0, 1, 0, 0, 0, 0
%    0, 0, 1, 0, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(3,6,6);%There are four "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F1;
% FhMats(:,:,4)=H*F3;
% FhMats(:,:,5)=H*F4;
% FhMats(:,:,6)=H*F4;
% c=299792458/1e3;%m/s The speed of light.
% 
% %For the demonstration, we will use noise-free TDOA measurements.
% z=zeros(6,1);
% sensorIdxLists=[1,1,1,2,1,1;
%                 2,2,3,3,2,3];
% for i=1:6
%    z(i)=(norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i)))-norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(2,i))))/c;
% end
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,c)
%
%EXAMPLE 3:
%A range-only example in 2D.
%We have three stationary sensors forming bistatic range measurements at
%four times.
% sensorStates=zeros(2,2);
% sensorStates(:,1)=[9e3;39e3];
% sensorStates(:,2)=[65e3;10e3];
% sensorStates(:,3)=[10e3;65e3];
% measTypes=ones(4,1);
% 
% xTrue=[27e3;42e3;300;0];%The true target location and velocity
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(4,4);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,4,1);%State transition matrix from time 0 to 1.
% F2=FPolyKal(2*T,4,1);%State transition matrix to time 2.
% F3=FPolyKal(3*T,4,1);%State transition matrix to time 3.
% H=[1, 0, 0, 0;
%    0, 1, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(2,4,4);%There are four "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F2;
% FhMats(:,:,4)=H*F3;
% 
% %For the demonstration, we will use noise-free range measurements.
% z=zeros(4,1);
% sensorIdxLists=[1,1,1,1;
%                 2,2,3,3];%Sensor 1 is always the transmitter.
% for i=1:2
%     z(i)=norm(FhMats(:,:,i)*xTrue-sensorStates(:,1))+norm(FhMats(:,:,i)*xTrue-sensorStates(:,2));
% end
% for i=3:4
%     z(i)=norm(FhMats(:,:,i)*xTrue-sensorStates(:,1))+norm(FhMats(:,:,i)*xTrue-sensorStates(:,3));
% end
% 
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates)
%
%EXAMPLE 4:
%Here we have range-only estimation in 3D. We have three sensors observing
%the target at four times.
% sensorStates=zeros(3,3);
% sensorStates(:,1)=[9e3;9e3;100];
% sensorStates(:,2)=[65e3;-10e3;1e3];
% sensorStates(:,3)=[64e3;-30e3;8e3];
% measTypes=ones(6,1);
% 
% xTrue=[27e3;42e3;5e3;300;0;0];%The true target location and velocity
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(6,6);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,6,1);%State transition matrix from time 0 to 1.
% F3=FPolyKal(3*T,6,1);%State transition matrix to time 3.
% F4=FPolyKal(4*T,6,1);%State transition matrix to time 4.
% H=[1, 0, 0, 0, 0, 0;
%    0, 1, 0, 0, 0, 0
%    0, 0, 1, 0, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(3,6,6);%There are six "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F1;
% FhMats(:,:,4)=H*F3;
% FhMats(:,:,5)=H*F4;
% FhMats(:,:,6)=H*F4;
% 
% %For the demonstration, we will use noise-free range measurements.
% z=zeros(6,1);
% sensorIdxLists=[1,1,1,2,1,1;
%                 2,2,3,3,2,3];
% for i=1:6
%     z(i)=norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i)))+norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(2,i)));
% end
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates)
%
%EXAMPLE 5:
%Here DOA-only measurements from three sensors at two times are used to
%estimate the state of a moving target with an assumed linear dynamic
%model. Since DOA measurements need a rotation matrix to define what is in
%"front" of the sensor, a somewhat detailed example near Hawaii is used.
% %Three DOA measurements.
% measTypes=2*ones(3,1);
% %Latitude/longitudes and ellipsoidal altitude of the sensors
% sensorLocs=[[20.765382;-155.978592;0],...
%             [20.231756;-155.761268;0],...
%             [20.126059;-155.555275;0]]*(pi/180);
% 
% %Latitude,longitude of the target.
% targetLatLon=[20.75;-155.5]*(pi/180);
% targetLatLonAlt=[targetLatLon;8000];%8km up in the air.
% targetLocCart=ellips2Cart(targetLatLonAlt);
% 
% %Say that the target is going North at 300m/s.
% u=getENUAxes(targetLatLonAlt);
% xTrue=[targetLocCart;300*u(:,2)];
% sensorStates=ellips2Cart(sensorLocs);
% 
% %Rotation matrices.
% %All sensors point 15 degrees up from the local horizonal.
% stdEl=15*(pi/180);
% East=90*(pi/180);
% North=0*pi/180;
% M1=findRFTransParam(sensorLocs(:,1),East,stdEl);
% M2=findRFTransParam(sensorLocs(:,2),North,stdEl);
% M3=findRFTransParam(sensorLocs(:,3),North,stdEl);
% HM=[1,0,0;
%     0,1,0];%Extract only the u-v components.
% 
% HMats=cell(3,1);
% HMats{1}=HM;
% HMats{2}=HM;
% HMats{3}=HM;
% 
% RMats=cell(3,1);
% RMats{1}=M1;
% RMats{2}=M2;
% RMats{3}=M3;
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(6,6);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,6,1);%State transition matrix from time 0 to 1.
% F2=FPolyKal(2*T,6,1);%State transition matrix to time 2
% H=[1, 0, 0, 0, 0, 0;
%     0, 1, 0, 0, 0, 0
%     0, 0, 1, 0, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(3,6,3);%There are three "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F2;
% 
% sensorIdxLists=[1,2,3;
%                 0,0,0];
% 
% %Create the DOA (u-v) measurements;
% z=zeros(6,1);
% zStart=1;
% for i=1:3
%     diff=FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i));
%     z(zStart:(zStart+1))=(1/norm(diff))*HMats{i}*RMats{i}*diff;
%     zStart=zStart+2;
% end
% 
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,[],HMats,RMats,0)
%
%EXAMPLE 6:
%Here we have monostatic range and DOA from two sensors at two times,
%observing a target with an assumed linear dynamic model. Since DOA
%measurements need a rotation matrix to define what is in "front" of the
%sensor, a somewhat detailed example near Hawaii is used. Note that it is
%not hard to find an explicit solution to this problem that is much faster
%than the solution computed here.
% %Two DOA measurements and two bistatic range measurements.
% measTypes=[2;2;1;1];
% %Latitude/longitudes and ellipsoidal altitude of the sensors
% sensorLocs=[[20.765382;-155.978592;0],...
%             [20.231756;-155.761268;0]]*(pi/180);
% 
% %Latitude,longitude of the target.
% targetLatLon=[20.75;-155.5]*(pi/180);
% targetLatLonAlt=[targetLatLon;8000];%8km up in the air.
% targetLocCart=ellips2Cart(targetLatLonAlt);
% 
% %Say that the target is going North at 300m/s.
% u=getENUAxes(targetLatLonAlt);
% xTrue=[targetLocCart;300*u(:,2)];
% sensorStates=ellips2Cart(sensorLocs);
% 
% %Rotation matrices.
% %All sensors point 15 degrees up from the local horizontal.
% stdEl=15*(pi/180);
% East=90*(pi/180);
% North=0*pi/180;
% M1=findRFTransParam(sensorLocs(:,1),East,stdEl);
% M2=findRFTransParam(sensorLocs(:,2),North,stdEl);
% HM=[1,0,0;
%     0,1,0];%Extract only the u-v components.
% 
% HMats=cell(4,1);
% HMats{1}=HM;
% HMats{2}=HM;
% HMats{3}=[];
% HMats{4}=[];
% 
% RMats=cell(4,1);
% RMats{1}=M1;
% RMats{2}=M2;
% RMats{3}=[];
% RMats{4}=[];
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(6,6);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,6,1);%State transition matrix from time 0 to 1.
% H=[1, 0, 0, 0, 0, 0;
%     0, 1, 0, 0, 0, 0
%     0, 0, 1, 0, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(3,6,4);%There are three "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F0;
% FhMats(:,:,4)=H*F1;
% 
% sensorIdxLists=[1,2,1,2;
%                 0,0,1,2];
% 
% %Create the DOA (u-v) measurements;
% z=zeros(6,1);
% zStart=1;
% for i=1:2
%     diff=FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i));
%     z(zStart:(zStart+1))=(1/norm(diff))*HMats{i}*RMats{i}*diff;
%     zStart=zStart+2;
% end
% 
% %Now, create the monostatic range measurements.
% for i=3:4
%     diff=FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i));
%     z(i+2)=2*norm(diff);
% end
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,[],HMats,RMats)
%
%EXAMPLE 7:
%Here DOA from two sensors at two times, and TDOA at two times using three
%sensors. The reference sensor does not make a DOA observation.
% %Two DOA measurements and two TDOA measurements.
% measTypes=[2;2;0;0];
% %Latitude/longitudes and ellipsoidal altitude of the sensors
% sensorLocs=[[20.765382;-155.978592;0],...
%             [20.231756;-155.761268;0],...
%             [20.126059;-155.555275;0]]*(pi/180);
% 
% %Latitude,longitude of the target.
% targetLatLon=[20.75;-155.5]*(pi/180);
% targetLatLonAlt=[targetLatLon;8000];%8km up in the air.
% targetLocCart=ellips2Cart(targetLatLonAlt);
% 
% %Say that the target is going North at 300m/s.
% u=getENUAxes(targetLatLonAlt);
% xTrue=[targetLocCart;300*u(:,2)];
% sensorStates=ellips2Cart(sensorLocs);
% 
% %Rotation matrices.
% %All sensors point 15 degrees up from the local horizonal.
% stdEl=15*(pi/180);
% East=90*(pi/180);
% North=0*pi/180;
% M1=findRFTransParam(sensorLocs(:,1),East,stdEl);
% M2=findRFTransParam(sensorLocs(:,2),North,stdEl);
% HM=[1,0,0;
%     0,1,0];%Extract only the u-v components.
% 
% HMats=cell(4,1);
% HMats{1}=HM;
% HMats{2}=HM;
% HMats{3}=[];
% HMats{4}=[];
% 
% RMats=cell(4,1);
% RMats{1}=M1;
% RMats{2}=M2;
% RMats{3}=[];
% RMats{4}=[];
% 
% %Now we construct the state transition matrices. A constant velocity
% %model is used.
% T=1;
% F0=eye(6,6);%State transition matrix at time 0 is just the identity
%             %matrix.
% F1=FPolyKal(T,6,1);%State transition matrix from time 0 to 1.
% H=[1, 0, 0, 0, 0, 0;
%     0, 1, 0, 0, 0, 0
%     0, 0, 1, 0, 0, 0];%Matrix to extract position components of the state.
% FhMats=zeros(3,6,4);%There are three "times".
% FhMats(:,:,1)=H*F0;
% FhMats(:,:,2)=H*F1;
% FhMats(:,:,3)=H*F0;
% FhMats(:,:,4)=H*F1;
% 
% sensorIdxLists=[1,2,3,3;
%                 0,0,1,2];
% c=299792458/1e3;%m/s The speed of light.
% 
% %Create the DOA (u-v) measurements;
% z=zeros(4,1);
% zStart=1;
% for i=1:2
%     diff=FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i));
%     z(zStart:(zStart+1))=(1/norm(diff))*HMats{i}*RMats{i}*diff;
%     zStart=zStart+2;
% end
% 
% %Now, create the TDOA measurements.
% for i=3:4
%     z(2+i)=(norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(1,i)))-norm(FhMats(:,:,i)*xTrue-sensorStates(:,sensorIdxLists(2,i))))/c;
% end
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,c,HMats,RMats,0)
%
%EXAMPLE 8:
%Here, we show that using the correct matrices, one can also do static
%localization as well.
% measTypes=[0;0;0];%All three TDOA
% sensorStates=zeros(3,4);
% sensorStates(:,1)=[9e3;9e3;100];
% sensorStates(:,2)=[65e3;-10e3;1e3];
% sensorStates(:,3)=[64e3;-30e3;8e3];
% sensorStates(:,4)=[-12e3;6e3;12];
% tTrue=[27e3;42e3;5e3];%The true target location.
% c=299792458;%m/s The speed of light.
% 
% FhMats=zeros(3,3,3);%There is only one "time".
% FhMats(:,:,1)=eye(3);
% FhMats(:,:,2)=eye(3);
% FhMats(:,:,3)=eye(3);
% FhMats(:,:,4)=eye(3);
% 
% z=zeros(3,1);
% sensorIdxLists=zeros(2,3);
% sensorIdxLists(:,1)=[1;2];%First TDOA pair.
% z(1)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,2)))/c;
% sensorIdxLists(:,2)=[1;3];%Second TDOA pair.
% z(2)=(norm(tTrue-sensorStates(:,1))-norm(tTrue-sensorStates(:,3)))/c;
% sensorIdxLists(:,3)=[4;3];%Third TDOA pair.
% z(3)=(norm(tTrue-sensorStates(:,4))-norm(tTrue-sensorStates(:,3)))/c;
% xEst=polyMeasConvertAsync(z,measTypes,FhMats,sensorIdxLists,sensorStates,c)
%
%REFERENCES:
%[1] D. F. Crouse, "General Multivariate Polynomial Target Localization and
%    Initial Estimation," Journal of Advances in Information Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(c))
   c=299792458;%m/s The speed of light in a vacuum.
end

if(nargin<9||isempty(algorithm))
    algorithm=0;
end

if(nargin<10)
    opts=[];
end

if(nargin<11||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<12||isempty(RelTol))
    RelTol=1e-8;
end 

if(nargin<13)
    scratchFolderPath=[];
end

if(nargin<14)
    execPath=[];
end

if(algorithm==0&&isempty(opts))%If using Bertini, set default tolerances
    %Default options when using Bertini
    opts=struct('SecurityMaxNorm',1e9,'EndpointFiniteThreshold',1e9,'PathTruncationThreshold',1e9,'FinalTol',1e-12);
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

    if(~isfield(opts,'FinalTol'))
        opts.FinalTol=1e-12;
    end
end

%We subtract out the locations of the sensors and then put them back in at
%the end. This should help with numerical sensitivity issues when
%localizing things on the ground but using ECEF coordinates. The offset is
%removed in the end. We then also scale everything with regard to the most
%distant sensor.
centerOfRegion=mean(sensorStates,2);
sensorStates=bsxfun(@minus,sensorStates,centerOfRegion);

scale=max(sqrt(sum(sensorStates.*sensorStates,1)));
sensorStates=sensorStates/scale;
numMeas=length(measTypes);
curZIdx=1;
for curMeas=1:numMeas
   if(measTypes(curMeas)==1)
       z(curZIdx)=z(curZIdx)/scale;
       curZIdx=curZIdx+1;
   elseif(measTypes(curMeas)==2)
       curZIdx=curZIdx+2;
   else
       curZIdx=curZIdx+1;
   end
end
c=c/scale;

xDim=size(FhMats,2);

%Labels for variables. The maximum number of letters ofr labels is why the
%number of variables is limited to 52.
abcList='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';

%Get the variable names for the components of the state
xList=cell(1,xDim);
for curDim=1:xDim
    xList{curDim}=['x',integer2CustomCharString(curDim,abcList)];
end

xPolys=cell(xDim,1);

curPoly=1;
zStart=1;
for curMeas=1:numMeas
    F=FhMats(:,:,curMeas);
    
    l1=sensorStates(:,sensorIdxLists(1,curMeas));

    curType=measTypes(curMeas);
    
    if(curType==0)
        zCur=z(zStart)*c;
        zStart=zStart+1;
    elseif(curType==1)
        curType=0;
        zCur=z(zStart);
        zStart=zStart+1;
    end
    
    switch(curType)
        case 0%TDOA or bistatic range measurement
            l2=sensorStates(:,sensorIdxLists(2,curMeas));
            u=(1/2)*(l1+l2);
            v=(1/2)*(l1-l2);
            delta=zCur/2;
            
            termMat=zeros(1+xDim,1+xDim+xDim*(xDim+1)/2);
            
            curTerm=1;
            
            %The constant term.
            constTerm=-delta^2*(dot(u,u)+dot(v,v))+delta^4+(u'*v)^2;
            if(constTerm~=0);
                termMat(1,curTerm)=constTerm;
                termMat(2:end,curTerm)=zeros(xDim,1);
                curTerm=curTerm+1;
            end

            %Add in the linear terms for 2*delta^2*u'*F*t-2*(t'*F'*v)*(u'*v)
            a=2*delta^2*u'*F-2*(u'*v)*v'*F;
            for curDim=1:xDim
                if(a(curDim)~=0)
                    termMat(1,curTerm)=a(curDim);
                    termMat(curDim+1,curTerm)=1;
                    curTerm=curTerm+1;
                end
            end
            
            %Add in the quadratic terms for (t'*F'*v)^2-delta^2*t'*F'*F*t
            A=F'*(v*v')*F-delta^2*(F'*F);
            sel=zeros(xDim,1);
            for curDim1=1:xDim
                %The diagonal term first.
                coeff=A(curDim1,curDim1);
                sel(curDim1)=2;
                if(coeff~=0)
                    termMat(1,curTerm)=coeff;
                    termMat(2:end,curTerm)=sel;
                    curTerm=curTerm+1;
                end
                sel(curDim1)=1;
                
                for curDim2=(curDim1+1):xDim 
                    %The off-diagonal terms.
                    coeff=2*A(curDim1,curDim2);
                    sel(curDim2)=1;
                    if(coeff~=0)
                        termMat(1,curTerm)=coeff;
                        termMat(2:end,curTerm)=sel;
                        curTerm=curTerm+1;
                    end
                    sel(curDim2)=0;
                end
                sel(curDim1)=0;
            end
            
            %Shink to fit.
            termMat=termMat(:,1:(curTerm-1));
            
            %Scale the term matrix.
            termMat(1,:)=termMat(1,:)/max(abs(termMat(1,:)));
            
            xPolys{curPoly}=termMat;
            curPoly=curPoly+1;
        case 2%A u-v DOA measurement
            zCur=z(zStart:(zStart+1));
            zStart=zStart+2;
            Hr=HMats{curMeas}*RMats{curMeas};
            u=zCur;
            
            %Two equations are produced by a u-v DOA measurement.
            termMat1=zeros(1+xDim,1+xDim+xDim*(xDim+1)/2);
            termMat2=zeros(1+xDim,1+xDim+xDim*(xDim+1)/2);
            
            curTerm1=1;
            curTerm2=1;
            
            %The constant term vector.
            constTerm=(u.*u)*dot(l1,l1)-(Hr*l1).^2;
            if(constTerm(1)~=0)
               termMat1(1,curTerm1)=constTerm(1);
               termMat1(2:end,curTerm1)=0;
               curTerm1=curTerm1+1;
            end
            
            if(constTerm(2)~=0)
               termMat2(1,curTerm2)=constTerm(2);
               termMat2(2:end,curTerm2)=0;
               curTerm2=curTerm2+1;
            end
            
            %The coefficient of the linear term 2*(Hr*F*t).*(Hr*l1)
            linTerm1=bsxfun(@times,2*(Hr*F),Hr*l1);
            
            %The coefficient of the linear term -2*(u.*u)*t'*F'*l1
            linTerm2=bsxfun(@times,-2*(u.*u),l1'*F);
            
            linTerm=linTerm1+linTerm2;
            
            for curDim=1:xDim
                if(linTerm(1,curDim)~=0)
                    termMat1(1,curTerm1)=linTerm(1,curDim);
                    termMat1(1+curDim,curTerm1)=1;
                    curTerm1=curTerm1+1;
                end
                
                if(linTerm(2,curDim)~=0)
                    termMat2(1,curTerm2)=linTerm(2,curDim);
                    termMat2(1+curDim,curTerm2)=1;
                    curTerm2=curTerm2+1;
                end
            end

            %We will add in the effects of the quadratic term
            %-(Hr*F*t).*(Hr*F*t)
            A=Hr*F;
            %and there is also a quadratic term (u.*u)*(t'*F'*F*t). We will
            %use two parts for this one.
            quadTerm21=(F'*F);
            uuTerm=(u.*u);
            
            for curDim1=1:xDim
                coeff1=quadTerm21(curDim1,curDim1)*uuTerm(1)-A(1,curDim1)^2;
                
                if(coeff1~=0)
                    termMat1(1,curTerm1)=coeff1;
                    termMat1(1+curDim1,curTerm1)=2;
                    curTerm1=curTerm1+1;
                end
                
                coeff2=quadTerm21(curDim1,curDim1)*uuTerm(2)-A(2,curDim1)^2;
                
                if(coeff2~=0)
                    termMat2(1,curTerm2)=coeff2;
                    termMat2(1+curDim1,curTerm2)=2;
                    curTerm2=curTerm2+1;
                end
                
                for curDim2=(curDim1+1):xDim
                    coeff1=2*quadTerm21(curDim1,curDim2)*uuTerm(1)-2*A(1,curDim1)*A(1,curDim2);
                    
                    if(coeff1~=0)
                        termMat1(1,curTerm1)=coeff1;
                        termMat1(1+curDim1,curTerm1)=1;
                        termMat1(1+curDim2,curTerm1)=1;
                        curTerm1=curTerm1+1;
                    end
                    
                    coeff2=2*quadTerm21(curDim1,curDim2)*uuTerm(2)-2*A(2,curDim1)*A(2,curDim2);
                    
                    if(coeff2~=0)
                        termMat2(1,curTerm2)=coeff2;
                        termMat2(1+curDim1,curTerm2)=1;
                        termMat2(1+curDim2,curTerm2)=1;
                        curTerm2=curTerm2+1;
                    end
                end
            end
            
            %Shrink to fit.
            termMat1=termMat1(:,1:(curTerm1-1));
            termMat2=termMat2(:,1:(curTerm2-1));
            
            %Scale the term matrices.
            termMat1(1,:)=termMat1(1,:)/max(abs(termMat1(1,:)));
            termMat2(1,:)=termMat2(1,:)/max(abs(termMat2(1,:)));
            
            xPolys{curPoly}=termMat1;
            xPolys{curPoly+1}=termMat2;
            curPoly=curPoly+2;
        otherwise
            error('Unknown measurement type specified')
    end
end


for curEq=1:xDim
    xPolys{curEq}=terms2String(xPolys{curEq},xList);
end
    
xEst=real(solvePolySysWithExtProg(xPolys,xList,algorithm,opts,scratchFolderPath,execPath));

posDim=size(sensorStates,1);

%Next, we check for the consistency of the converted measurements to remove
%false solutions.
numSol=size(xEst,2);
sel=true(numSol,1);
zStart=1;
for curMeas=1:numMeas
    switch(measTypes(curMeas))
        case 0%TDOA
            lRx1=sensorStates(:,sensorIdxLists(1,curMeas));
            lRx2=sensorStates(:,sensorIdxLists(2,curMeas));
            zCur=z(zStart);
            zStart=zStart+1;
            
            for curSol=1:numSol
                TDOAComp=(norm(FhMats(:,:,curMeas)*xEst(:,curSol)-lRx1)-norm(FhMats(:,:,curMeas)*xEst(:,curSol)-lRx2))/c;
                absDiff=abs(TDOAComp-zCur);
                
                if(~(absDiff<AbsTol || absDiff<RelTol*abs(zCur)))
                    sel(curSol)=false;
                end
            end
        case 1%Range
            lRx1=sensorStates(:,sensorIdxLists(1,curMeas));
            lRx2=sensorStates(:,sensorIdxLists(2,curMeas));
            zCur=z(zStart);
            zStart=zStart+1;
            
            for curSol=1:numSol
                rComp=norm(FhMats(:,:,curMeas)*xEst(:,curSol)-lRx1)+norm(FhMats(:,:,curMeas)*xEst(:,curSol)-lRx2);
                absDiff=abs(rComp-zCur);
                
                if(~(absDiff<AbsTol || absDiff<RelTol*abs(zCur)))
                    sel(curSol)=false;
                end
            end
        case 2%DOA
            lRx1=sensorStates(:,sensorIdxLists(1,curMeas));
            F=FhMats(:,:,curMeas);
            H=HMats{curMeas};
            R=RMats{curMeas};
            zCur=z(zStart:(zStart+1));
            zStart=zStart+2;
            
            for curSol=1:numSol
                diff=F*xEst(:,curSol)-lRx1;
                uComp=(1/norm(diff))*R*diff;
                
                %If the DOA is behind the sensor
                if(uComp(3)<=0)
                    sel(curSol)=false;
                    continue;
                end
                
                absDiff=abs(H*uComp-zCur);
                if(~all(absDiff<AbsTol))
                    sel(curSol)=false;
                end
            end
        otherwise%No tests for other measurement types.
    end
end

%Get rid of invalid solutions.
xEst=xEst(:,sel);

if(~isempty(xEst))
    %Add back in the offset that was present due to centering the coordinate
    %system around the sensors.
    xEst=xEst*scale;
    xEst(1:posDim,:)=bsxfun(@plus,xEst(1:posDim,:),centerOfRegion);
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
