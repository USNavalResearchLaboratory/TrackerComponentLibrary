function PV=readJPLEphem(TDB1,TDB2,objectNumber,centerNumber,units,ephemerisPath)
%%READJPLEPHEM Read information, such as the location and/or velocity of a
%       planet, from one of the National Aeronautics and Space
%       Administration (NASA) Jet Propulsion Laboratory's (JPL's)
%       Development Ephemerides (DE). These are typically referred to using
%       DE followed by a number and sometimes a letter, such as DE200 or
%       DE430t. This function only reads ephemerides given in binary format
%       for Linux (but will work under Mac OS X and Windows). The
%       coordinate system in which the locations are given varies depending
%       on the ephemeris chosen. The DE430, DE430t and DE431 ephemerides
%       are given with respect to the International Celestial Reference
%       Frame 2 (ICRF2). That is, the orientation of the axes is with
%       respect to the International Celestial Reference System (ICRS).
% 
%INPUTS: TDB1, TDB2 Two parts of the Julian date at which the ephemeris
%           value is desired, given in barycentric dynamical time. The full
%           date is the sum of both terms. The date is broken into two
%           parts to provide more bits of precision. It does not matter how
%           the date is split.
% objectNumber, centerNumber Numbers indicating which object or point's
%           position/velocity is desired and with respect to which
%           other object or point it is to be taken. For values that are
%           not positions or velocities (e.g. TDB-TT), only objectNumber is
%           used; centerNumber can be an empty matrix. All options in terms
%           of solar system bodies return [x;y;z;xDot;yDot;zDot] position
%           and velocity. The available values vary from one ephemeris set
%           to another. Some possible values are:
%           1  Mercury
%           2  Venus 
%           3  Earth (geocenter)
%           4  Mars (system barycenter)
%           5  Jupiter (system barycenter)
%           6  Saturn (system barycenter)
%           7  Uranus (system barycenter)
%           8  Neptune (system barycenter)
%           9  Pluto (system barycenter)
%           10 Moon
%           11 Sun
%           12 Solar System Barycenter
%           13 Earth-Moon Barycenter
%           14 1980 IAU nutation angles (longitude, obliquity, longitue
%                                       rate and obliquity rate)
%           15 Lunar libration (Euler) angles (phi, theta, psi, dphi/dt,
%                                       dtheta/dt, and dpsi/dt)
%           16 Lunar angular velocity (omega_x, omega_y, omega_y and
%                                       domega_x/dt, domega_y/dt,
%                                       domega_z/dt)
%           17 TT-TDB (at geocenter) The value of TT-TDB and its rate.
%     units An optional parameter selecting the units of the returned
%           values. Possible values are
%           0 (The default if omitted or an empty matrix is passed) For
%             position, use meters; for time use seconds.
%           1 For position use astronomical units (AU); for time use Julian
%             days (TDB). The length of 1AU in kilometers can be queried
%             from the ephemerides using the function readJPLEphemHeader.
% ephemerisPath This is a character string of the path to the ephemeris
%             file. The formatting is the same as used by the fopen
%             function. If this parameter is omitted or an empty matrix is
%             passed, it is assumed that the DE430t ephemeris is given in a
%             data folder that is in the same folder as this function and
%             the file is named linux_p1550p2650.430t.
%
%OUTPUTS: PV If The relative location and velocity between two solar system
%            objects is requested, then this is a 6X1 vector in the format
%            [x;y;z;xDot;yDot;zDot]. If other types of information are
%            requested, then the formatting is as described above. The
%            units depend on the units input.
%
%Ephemerides that can be read by this function can be downloaded from
%ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux
%Only the binary file is needed. For example, when considering the DE430
%ephemerides, the file that should be used is called linux_p1550p2650.430.
%Accompanying files that are not needed include testpo.430 and
%header.430_572. The file header.430_572 is an ASCII test version of the
%header and testpo.430 is an ASCII text file containing parameters related
%to testing the algorithm.
%
%This file is loosely based on the Fortran routine called testeph1.f that
%was created by the California Institute of Technology (CIT) under a U.S.
%government contract with NASA and is available for download in the folder
%ftp://ssd.jpl.nasa.gov/pub/eph/planets/
%
%EXAMPLE:
%Find the position and velocity of the Earth with respect to the Sun with
%axes oriented according to the ICRS.
% objectNumber=3;
% centerNumber=11;
% 
% %A TDB time that happens to be the reference epoch for the DE430t
% %ephemeris.
% TDB1=2440400;
% TDB2=0.5;
% %Get position and velocity in units of meters and meters per second.
% units=0;
% posVel=readJPLEphem(TDB1,TDB2,objectNumber,centerNumber,units)
% %Get position and velocity in units of astronomical units and astronomical
% %units per TDB Julian day.
% units=1;
% posVel=readJPLEphem(TDB1,TDB2,objectNumber,centerNumber,units)
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If no ephemeris is given, assume that the DE430t is in a folder named data
%that is located in the same folder as this file.
if(nargin<6||isempty(ephemerisPath)) 
%Get the path to this file.
    ScriptPath=mfilename('fullpath');
    ScriptFolder=fileparts(ScriptPath);
    ephemerisPath=[ScriptFolder,'/data/linux_p1550p2650.430t'];
end

if(nargin<5||isempty(units))
   units=0;%Use kilometers and seconds 
end

%%%%%
%%%First, validate the inputs
%%%%%

if(objectNumber<=0||objectNumber>17)
  error('objectNumber must be between 1 and 17.')  
end

%centerNumber is not used if nutations, librations, libration rates, or
%TT-TDB are requested. Those are all centerNumber>=14.
if(objectNumber<14&&(centerNumber<=0||centerNumber>13))
    error('centerNumber must be between 1 and 13..')
end

%If the coordinate system coincides with the object whose location is
%desired.
if(objectNumber<=13&&objectNumber==centerNumber)
    PV=zeros(6,1);
    return;
end

%%%%%
%%%Read the header information of the ephemerides.
%%%%%
[~,~,kSize,ss,ipt,au_km,EarthMoonRatio]=readJPLEphemHeader(ephemerisPath);
%The number of Chebyshev coefficients per record.
numCoeffs=kSize/2;

startJulEphemDate=ss(1);
endJulEphemDate=ss(2);
recordLengthDays=ss(3);%(Given in Julian days)

%This number of seconds per day is used in the ephemerides.
secondsPerDay=86400;
%The number of seconds per record.
secondSpan=secondsPerDay*recordLengthDays;

if(units)
    %The scale for position so that results are obtained in astronomical
    %units.
    xScale=1/au_km;
    %The scale for the velocity so that results are obtained in terms of
    %kilometers per second.
    vScale=xScale*secondsPerDay;
else
    %The scale for position so that results are obtained in meters.
    xScale=1000;
    %The scale for the velocity so that results are obtained in terms of
    %meters per second.
    vScale=1000;
end

%%%%%
%%%Determine which data record is required.
%%%%%
sumT=TDB1+TDB2;

if(sumT<startJulEphemDate)
    error('The input time is before the earliest time in the ephemerides.')
end
    
if(sumT>endJulEphemDate)
    error('The input time is after the latest time in the ephemerides.')
end
    
recordIdx=fix((sumT-startJulEphemDate)/recordLengthDays)+3;
if(sumT==endJulEphemDate)
    recordIdx=recordIdx-1;
end

fid=fopen(ephemerisPath,'r');
fseek(fid,4*kSize*(recordIdx-1),'bof');
data=fread(fid,numCoeffs,'double');
fclose(fid);

TF1=double(fix(TDB1));
TF2=TDB1-TF1;
TF=(TF1-data(1))/recordLengthDays;
TF=TF+((TF2+TDB2)/recordLengthDays);

%If we are very close to the end of an interval, so that it got pushed to
%the next interval. An example of where this type of problem occurs is when
%using the linux_p1550p2650.430t ephemerides and setting the time to
%TDB1=2440400+32; TDB2=0.5-1000*eps(32.5);
if(TF<0)
    if(recordIdx==3)
        error('The input time is before the earliest time in the ephemerides.')
    else
        %We should be at the end of the previous record.
        recordIdx=recordIdx-1;
        fid=fopen(ephemerisPath,'r');
        fseek(fid,4*kSize*(recordIdx-1),'bof');
        data=fread(fid,numCoeffs,'double');
        fclose(fid);

        TF1=double(fix(TDB1));
        TF2=TDB1-TF1;
        TF=(TF1-data(1))/recordLengthDays;
        TF=TF+((TF2+TDB2)/recordLengthDays);
    end
end

%%%%
%%Consider lookups for nutations, librations, or TDB-TT, which do not
%%use the centerNumber input. 
%%%%

%Nutation
if(objectNumber==14)
    if(ipt(1,12)>0&&ipt(2,12)*ipt(3,12)~=0)
        NCF=ipt(2,12);
        NCM=2;
        NSC=ipt(3,12);
        numData=NCF*NCM*NSC;
        dataCur=reshape(data(ipt(1,12):(ipt(1,12)+numData-1)),[NCF,NCM,NSC]);
        
        PV=INTCHB(dataCur,TF,secondSpan);

        if(units)
            %Convert radian/second to radian/day
            PV(3:4)=secondsPerDay*PV(3:4);
        end
        return
    else
        warning('Requested nutations not found.')
        PV=[];
        return
    end
end

%Lunar mantle Euler angles
if(objectNumber==15)
    if(ipt(1,13)>0&&ipt(2,13)*ipt(3,13)~=0)
        NCF=ipt(2,13);
        NCM=3;
        NSC=ipt(3,13);
        numData=NCF*NCM*NSC;
        dataCur=reshape(data(ipt(1,13):(ipt(1,13)+numData-1)),[NCF,NCM,NSC]);

        PV=INTCHB(dataCur,TF,secondSpan);
        
        if(units)
            %Convert radians/(second*day) to radians/day^2
            PV(4:6)=secondsPerDay*PV(4:6);
        else
            %Convert radians/day to radians/second
            PV(1:3)=PV(1:3)/secondsPerDay;
            %Convert radians/(second*day) to radians/second^2
            PV(4:6)=PV(4:6)/secondsPerDay;
        end
        return
    else
        warning('Requested lunar mantle Euler angles not found.')
        PV=[];
        return
    end
end

%Lunar mantle angular velocity
if(objectNumber==16)
    if(ipt(1,14)>0&&ipt(2,14)*ipt(3,14)~=0)
        NCF=ipt(2,14);
        NCM=3;
        NSC=ipt(3,14);
        numData=NCF*NCM*NSC;
        dataCur=reshape(data(ipt(1,14):(ipt(1,14)+numData-1)),[NCF,NCM,NSC]);
        
        PV=INTCHB(dataCur,TF,secondSpan);
        
        if(units)
            %Convert radian/second to radian/day
            PV(4)=secondsPerDay*PV(4);
            PV(5)=secondsPerDay*PV(5);
            PV(6)=secondsPerDay*PV(6);
        end
        return
    else
        warning('Requested lunar mantle angular velocity not found.')
        PV=[];
        return
    end
end

%TT-TDB
if(objectNumber==17)
    if(ipt(1,15)>0&&ipt(2,15)*ipt(3,15)~=0)
        NCF=ipt(2,15);
        NCM=1;
        NSC=ipt(3,15);
        numData=NCF*NCM*NSC;
        dataCur=reshape(data(ipt(1,15):(ipt(1,15)+numData-1)),[NCF,NCM,NSC]);
        PV=INTCHB(dataCur,TF,secondSpan);

        if(units)
            %Convert seconds/second to seconds/day
            PV(2)=secondsPerDay*PV(2);
        end
        return
    else
        warning('Requested TT-TDB not found.')
        PV=[];
        return
    end
end

%%%%
%%Consider lookups for ceklestial bodies, which use the centerNumber input.
%%%%

%If the object is the Moon and the center is the Earth
if(objectNumber==10&&centerNumber==3)
    NCF=ipt(2,10);
    NCM=3;
    NSC=ipt(3,10);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,10):(ipt(1,10)+numData-1)),[NCF,NCM,NSC]);

    PV=INTCHB(dataCur,TF,secondSpan);
    PV(1:3)=PV(1:3)*xScale;
    PV(4:6)=PV(4:6)*vScale;
    return
end

%If the object is the Earth and the center is the Moon
if(objectNumber==3&&centerNumber==10)
    NCF=ipt(2,10);
    NCM=3;
    NSC=ipt(3,10);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,10):(ipt(1,10)+numData-1)),[NCF,NCM,NSC]);
    PV=INTCHB(dataCur,TF,secondSpan);
    
    PV(1:3)=-PV(1:3)*xScale;
    PV(4:6)=-PV(4:6)*vScale;
    return
end

PV1=zeros(6,1);
PV2=zeros(6,1);

%If the object is the Moon or the Earth OR the center is the Moon or the
%Earth.
if(objectNumber==3||objectNumber==10||centerNumber==3||centerNumber==10)
    NCF=ipt(2,10);
    NCM=3;
    NSC=ipt(3,10);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,10):(ipt(1,10)+numData-1)),[NCF,NCM,NSC]);
    PVM=INTCHB(dataCur,TF,secondSpan);

    NCF=ipt(2,3);
    NCM=3;
    NSC=ipt(3,3);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,3):(ipt(1,3)+numData-1)),[NCF,NCM,NSC]);
    PVB=INTCHB(dataCur,TF,secondSpan);
    
    FACTE=1/(1+EarthMoonRatio);
    FACTM=FACTE-1;
    
    %If it is the object that is the Moon or the Earth
    if(objectNumber==3||objectNumber==10) 
        if(centerNumber<12)
            NCF=ipt(2,centerNumber);
            NCM=3;
            NSC=ipt(3,centerNumber);
            numData=NCF*NCM*NSC;
            dataCur=reshape(data(ipt(1,centerNumber):(ipt(1,centerNumber)+numData-1)),[NCF,NCM,NSC]);

            PV2=INTCHB(dataCur,TF,secondSpan);
        elseif(centerNumber==13)
            PV2=PVB;
        end

        if(objectNumber==3)
            PV=PVB-FACTE*PVM-PV2;
        else
            PV=PVB-FACTM*PVM-PV2;
        end
    else 
        %If here, then it is the center that is the Moon or the Earth.
        if(objectNumber<12)
            NCF=ipt(2,objectNumber);
            NCM=3;
            NSC=ipt(3,objectNumber);
            numData=NCF*NCM*NSC;
            dataCur=reshape(data(ipt(1,objectNumber):(ipt(1,objectNumber)+numData-1)),[NCF,NCM,NSC]);

            PV1=INTCHB(dataCur,TF,secondSpan);
        elseif(objectNumber==13)
            PV1=PVB;
        end

      if(centerNumber==3)
          PV=PV1-(PVB-FACTE*PVM);
      else
          PV=PV1-(PVB-FACTM*PVM);
      end
    end
    PV(1:3)=PV(1:3)*xScale;
    PV(4:6)=PV(4:6)*vScale;
    return
end
      
%If here, neither the object nor the center is the Moon or the Earth,
%though it could be the Earth-Moon barycenter.

if(objectNumber<12)
    NCF=ipt(2,objectNumber);
    NCM=3;
    NSC=ipt(3,objectNumber);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,objectNumber):(ipt(1,objectNumber)+numData-1)),[NCF,NCM,NSC]);
    
    PV1=INTCHB(dataCur,TF,secondSpan);
elseif(objectNumber==13)
    NCF=ipt(2,3);
    NCM=3;
    NSC=ipt(3,3);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,3):(ipt(1,3)+numData-1)),[NCF,NCM,NSC]);

    PV1=INTCHB(dataCur,TF,secondSpan);
end

if(centerNumber<12)
    NCF=ipt(2,centerNumber);
    NCM=3;
    NSC=ipt(3,centerNumber);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,centerNumber):(ipt(1,centerNumber)+numData-1)),[NCF,NCM,NSC]);
    
    PV2=INTCHB(dataCur,TF,secondSpan);
elseif(centerNumber==13)
    NCF=ipt(2,3);
    NCM=3;
    NSC=ipt(3,3);
    numData=NCF*NCM*NSC;
    dataCur=reshape(data(ipt(1,3):(ipt(1,3)+numData-1)),[NCF,NCM,NSC]);

    PV2=INTCHB(dataCur,TF,secondSpan);
end

PV=PV1-PV2;
PV(1:3)=PV(1:3)*xScale;
PV(4:6)=PV(4:6)*vScale;
end

function PV=INTCHB(BUF,T,LINT)
%%INTCHB This is the INTCHB function which takes the buffer of Chebyshev
%        coefficients for the ephemeris data and a time and returns the
%        interpolated position and velocity. This is analogous to the
%        INTCHB function in testeph1.f. However, to make the function more
%        understandable, the functions ChebyshevPolySynth and
%        ChebyshevPolyDerivCoeffs are used to manipulate the Chebyshev
%        polynomials rather than directly synthesizing them in here. The
%        function always returns a derivative value (e.g. position and
%        velocity).
%
%INPUTS: BUF The NCF X NCM X NSC set of Chebychev coefficients. NCF is the
%            number of coefficients per component. NCM is the number of
%            components per set of components (for synthesizing position
%            and velocity, this is the dimensionality of space, e.g. 3D).
%            NSC is the number of sets of coefficients within the
%            interpolation interval.
%          T The fractional time within the interval where interpolation
%            should be performed 0<=T<=1.
%       LINT The length of the interval in time units (e.g. seconds). This
%            is needed for compouting the velocity.
%
%OUTPUTS: PV A (2*NCM)X1 vector of the NCMX1 desired value and its NCMX1
%            derivative.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %NCM is the number of components per set of coefficients.
    NCM=size(BUF,2);
    %NSC is the number of sets of coefficients within the interval.
    NSC=size(BUF,3);
       
    %Allocate space
    PV=zeros(2*NCM,1);
    
    %Compute set number within interval (L), and scaled Chebyshev time
    %within that set (TC).
    NS=NSC;
    TEMP=T*NS;
    TC=2*(TEMP-fix(TEMP))-1;
    L=fix(TEMP)+1;

    if(L>NS) 
        L=L-1;
        TC=TC+2;
    end

    BMA=2*NS/LINT;
    for i = 1:NCM
        %The value (often position).
        PV(i)=ChebyshevPolySynth(TC,BUF(:,i,L));
        
        a=ChebyshevPolyDerivCoeffs(BUF(:,i,L));
        %The derivative of the value (often velocity).
        PV(i+NCM)=BMA*ChebyshevPolySynth(TC,a);
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
