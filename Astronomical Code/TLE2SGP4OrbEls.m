function [SGP4Elements,TTEpoch1,TTEpoch2,otherInfo,checksumIsBad]=TLE2SGP4OrbEls(TLELine1,TLELine2)
%%TLE2SGP4ORBELS Convert a North American Aerospace Defense Command (NORAD)
%                two-line element (TLE) satellite ephemerides into a format
%                that can be used with the propagateOrbitSGP4 function for
%                determining the position and velocity of the satellite at
%                different times using the Simplified General
%                Perturbations 4 (SGP4) model/ the Simplified Deep Space
%                Perturbations 4 (SDP4) model for distant satellites.
%                TLEs are how the Air Force Space Command (AFSPC) has
%                publicly released information on satellite orbits for
%                years.
%
%INPUTS: TLELine1, TLELine2 The character strings for each of the two lines
%                           that make up the TLE ephemeride, including all
%                           spaces as the positions of the characters are
%                           important. With checksums, these should both be
%                           69 characters long, without checksums, 68
%                           characters long.
%
%OUTPUTS: SGP4Elements The 7X1 vector of orbital elements that can be given
%                      to the propagateOrbitSGP4 function along with the
%                      epoch time to predict the location of the satellite
%                      at various times. They are similar, but not
%                      identical to, Keplerian orbital elements of the same
%                      name. If either of the TLE strings is too short,
%                      this is an empty matrix. The elements of
%                      SGP4Elements are
%                      1) eccentricity
%                      2) inclination  (radians)
%                      3) argument of perigee (radians)
%                      4) right ascension of the ascending node (radians)
%                      5) mean anomaly (radians)
%                      6) mean motion (radians per second [TT])
%                      7) BSTAR drag term. This pseudo ballistic
%                        coefficient has units of inverse Earth radii.
%                        Normally, a ballistic coefficient BC is mass per
%                        unit area divided by a drag coefficient. The BSTAR
%                        drag term is a scaled version of the normal
%                        ballistic coefficient. That is,
%                        BC=Re*rho_0/(2*BSTAR) where Re is the radius of
%                        the Earth in meters, and rho=2.461e-5 kg/m^2.
%   TTEpoch1, TTEpoch2 The epoch time of the SGP4Elements given in
%                      terrestrial time (TT) as a two-part Julian date.
%                      This will be incorrect from 2057 onward, since the
%                      TLE data only provides the final two digits of the
%                      year starting in 1957 and the century must be
%                      guessed. If either of the TLE strings is too short,
%                      these are empty matrices.
%            otherInfo A structure holding the other information present in
%                      the TLE. If either of the TLE strings is too short,
%                      this is an empty matrix. The fields of this
%                      structure are
%                      line1Number The number starting the first TLE line.
%                                  This should be 1.
%                      line2Number The number starting the second TLE line.
%                                  This should be 2.
%                      satelliteNumber1 The NORAD catalog number of the
%                                  satellite as given by the first line of
%                                  the TLE.
%                      satelliteNumber2 The NORAD catalog number of the
%                                  satellite as given by the second line of
%                                  the TLE. This should be the same as
%                                  satelliteNumber1
%                      classification A character indicating the
%                                  classification level of the satellite.
%                                  'U'=unclassified.
%                      IDLaunchYear The last two digits of the launch year
%                                  of the satellite.
%                      IDLaunchNumber The number of the launch in in the
%                                  launch year that put up the satellite.
%                      IDPieceOfLaunch A string indicating which piece of
%                                  the launch the satellite is.
%                      firstDerivMeanMotion The first derivative of the
%                                  mean motion in radians per second
%                                  squared. This is useful with the older
%                                  orbital propagators.
%                      secondDerivMeanMotion The second derivative of the
%                                  mean motion in radians per second cubed.
%                      ephemerisType Originally, this was a number from
%                                  1-5 indicating the type of orbital
%                                  propagator used to generate the
%                                  ephemeris with 1=SGP, 2=SGP4, 3=SDP4,
%                                  4=SGP8, and 5=SDP8. However, it is
%                                  always set to zero in TLEs that are
%                                  publicly released and is thus
%                                  meaningless.
%                      elementNumber A number indicating which TLE set this
%                                  TLE is from. The number is supposed to
%                                  be incremented everytime a new set comes
%                                  out. However, in practice, it is not.
%                      revolutionNumberAtEpoch The number of orbits the
%                                  satellite has made since its launch.
%        checksumIsBad A 2X1 boolean vector indicating whether the checksum
%                      value for each of the strings is bad.
%                      checksumIsBad(1) is true if the checksum of TLELine1
%                      is invalid and checksumIsBad(2) is true if the
%                      checksum of TLELine2 is invalid. If a checksums is
%                      omitted then the value will be true. If iether of
%                      the strings is too short, then this will be an empty
%                      matrix.
%
%Note that parts of the input strings that cannot be read will be filled
%with NaN values.
%
%Information of the format of the TLE sets can be found at [1]. The
%documentation discussing how the element sets are used is in [2] and [3].
%The decimal point parsing part of this function is loosely based on David
%Vallado's public-domain code for the SGP4 propagator that was downloaded
%from http://www.centerforspace.com/downloads/
%
%A sample TLE that accompanies Vallado's code is:
%TLELine1='1 11801U          80230.29629788  .01431103  00000-0  14311-1      13'
%TLELine2='2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13'
%
%REFERENCES:
%[1] T. C. for Space Standards and Innovation. (2012, 27 Sep.) NORAD two-
%    line element sets. [Online].
%    Available: http://www.celestrak.com/NORAD/elements/
%[2] F. R. Hoots and R. L. Roehrich, "Spacetrack report no. 3: Models for
%    propagation of NORAD element sets," Department of Defense, Tech. Rep.,
%    31 Dec. 1988. [Online].
%    Available: http://www.amsat.org/amsat/ftp/docs/spacetrk.pdf
%[3] D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, "Revisiting
%    spacetrack report # 3: Rev 2," in Proceedings of the AIAA/AAS
%    Astrodynamics Specialist Conference and Exhibit, Keystone, CO, 21-24
%    Aug. 2006. [Online].
%    Available: http://celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.pdf
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Multiplication coefficient to convert from degrees to radians.
deg2Rad=pi/180;

%If either of the strings is too short.
if(length(TLELine1)<68||length(TLELine2)<68)
    SGP4Elements=[];
    TTEpoch1=[];
    TTEpoch2=[];
    otherInfo=[];
    checksumIsBad=[];
    return;
end

%If the strings have no checksums at the end, then just mark the checksums
%as bad.
if(length(TLELine1)<69)
    checksumIsBad(1)=true;
else
    checksumIsBad(1)=~validateTLEChecksum(TLELine1);
end
    
if(length(TLELine2)<69)
    checksumIsBad(2)=true;
else
    checksumIsBad(2)=~validateTLEChecksum(TLELine2);
end

%Set the implied decimal points to help deal with bad values when doing a
%formatted read.
for j=11:16
    if(TLELine1(j)==' ')
        TLELine1(j)='_';
    end
end

if(TLELine1(45)~=' ')
    TLELine1(44)=TLELine1(45);
end
TLELine1(45)='.';

if(TLELine1(8)==' ')
    TLELine1(8)='U';
end

if(TLELine1(10)==' ')
    TLELine1(10)='.';
end

for j=46:50
    if (TLELine1(j)==' ')
        TLELine1(j)='0';
    end
end

if(TLELine1(52)==' ')
    TLELine1(52)='0';
end

if(TLELine1(54)~=' ')
    TLELine1(53)=TLELine1(54);
end
TLELine1(54)='.';
TLELine2(26)='.';

for j=27:33
    if(TLELine2(j)==' ')
        TLELine2(j)='0';
    end
end

if(TLELine1(63)==' ')
    TLELine1(63)='0';
end

if((length(TLELine1)<68)||(TLELine1(68)==' '))
    TLELine1(68)='0';
end

%Parse the stuff that most folks do not care about.
otherInfo.line1Number = str2double(TLELine1(1));
otherInfo.line2Number = str2double(TLELine2(1));
otherInfo.satelliteNumber1 = str2double(TLELine1(3:7));
otherInfo.satelliteNumber2 = str2double(TLELine2(3:7));

otherInfo.classification = TLELine1(8);
otherInfo.IDLaunchYear = str2double(TLELine1(10:11));
otherInfo.IDLaunchNumber= str2double(TLELine1(12:14));
otherInfo.IDPieceOfLaunch= TLELine1(15:17);

%Convert revolutions per day squared into radians per second squared,
%with 86400 seconds per TT Julian day.
otherInfo.firstDerivMeanMotion=str2double(TLELine1(34:43))*(2*pi/86400^2);
%Convert revolutions per day cubed into radians per second cubed,
%with 86400 seconds per TT Julian day.
otherInfo.secondDerivMeanMotion=str2double(TLELine1(44:50))*10.0^str2double(TLELine1(51:52))*(2*pi/86400^3); 
otherInfo.ephemerisType = str2double(TLELine1(63));
otherInfo.elementNumber = str2double(TLELine1(65:68));
otherInfo.revolutionNumberAtEpoch=str2double(TLELine2(64:68));

%Parse the stuff that most folks care about
inclination = str2double(TLELine2(8:16))*deg2Rad;
RAAscendingNode = str2double(TLELine2(17:25))*deg2Rad;
eccentricity = str2double(TLELine2(26:33));
argumentOfPerigee = str2double(TLELine2(34:42))*deg2Rad;
meanAnomaly = str2double(TLELine2(43:51))*deg2Rad;
%Convert from revolutions per day to radians per second, with 86400
%seconds per TT Julian day.
meanMotion=str2double(TLELine2(52:63))*(2*pi/86400); 
BSTARDrag= str2double(TLELine1(53:59))*10.0^str2double(TLELine1(60:61));

%Get the time of the orbital elements.
epochYear = str2double(TLELine1(19:20));
epochDays = str2double(TLELine1(21:32));

%Assume that the years only run from 1957->2057 and find the actual
%epoch year.
if (epochYear < 57)
    year=epochYear+2000;
else
    year=epochYear+1900;
end
%Convert into a two-part Julian date in TT.
[month,dayOfMonth]=dayOfYear2MonthDay(year,fix(epochDays));
dayFrac=epochDays-fix(epochDays);
[hour,minute,second]=fracDayOfMonth2HourMinSec(year,month,dayOfMonth,dayFrac);
[TTEpoch1,TTEpoch2]=Cal2TT(year,month,dayOfMonth,hour,minute,second);

SGP4Elements=zeros(7,1);
SGP4Elements(1)=eccentricity;
SGP4Elements(2)=inclination;
SGP4Elements(3)=argumentOfPerigee;
SGP4Elements(4)=RAAscendingNode;
SGP4Elements(5)=meanAnomaly;
SGP4Elements(6)=meanMotion;
SGP4Elements(7)=BSTARDrag;
end

function isValid=validateTLEChecksum(theTLELine)
%The validation comes from the explanation of how the checksum is formed
%from http://www.celestrak.com/NORAD/elements/ It is the last digit of the
%sum of all of the digits in the string plus an extra 1 for each '-'. 
    theSum=0;
    for curCharIdx=1:68
        curChar=theTLELine(curCharIdx);
        if(curChar=='-')
            theSum=theSum+1;
        elseif(curChar>='0'&&curChar<='9')
            theSum=theSum+str2double(curChar);
        end
    end

    %Get the last digit of the sum.
    sumStr=num2str(theSum);
    isValid=(sumStr(end)==theTLELine(69));
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
