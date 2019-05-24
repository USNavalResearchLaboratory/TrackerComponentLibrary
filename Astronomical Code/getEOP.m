function [xpyp,dXdY,deltaUTCUT1,deltaTTUT1,LOD]=getEOP(JulUTC1,JulUTC2,refreshFromSource,replaceEOPtxt)
%%GETEOP Get Earth orientation parameters from a table for a particular
%        date, using interpolation for points between or outside of the
%        range of tabulated values. The data can be read from
%        ./data/EOP.txt or downloaded from the internet. Once the data have
%        been loaded, they do not have to be reloaded for future function
%        calls. The Earth orientation parameters are necessary for
%        high-precision astronomical coordinate conversions. When the date
%        provided is outside of the range of tabulated values, then linear
%        interpolation is provided over deltaUTCUT1 and deltaTTUT1, and
%        piecewise cubic Hermite interpolation is provided over LOD, but
%        xpyp and dXdY are just set to zero.
%
%INPUTS: JulUTC1,JulUTC2 NX1 or 1XN vectors of two-part pseudo-Julian dates
%                    given in UTC. The units of the date are days. The full
%                    date is the sum of both terms. The date is broken into
%                    two parts to provide more bits of precision. It does
%                    not matter how the date is split.
%  refreshFromSource An optional parameter. If this parameter is given, a
%                    database of Earth orientation parameters is loaded
%                    from the selected source, replacing any previously
%                    loaded database. As the external sources are updated
%                    no more than once per day, this parameter should NOT
%                    be provided often, so as to minimize internet traffic.
%                    If this parameter is anything but -1, a warning is
%                    given. If available, the more accurate Bulletin B
%                    values are used. Otherwise, Bulletin A values are
%                    used. Possible values of this parameter are:
%                    -1) (The default) Load the data from the local file
%                       ./data/EOP.txt. No internet connection is required.
%                    0) Get the finals2000A.daily file from the US Naval
%                       Observatory. This has EOP parameters for the last
%                       90 days and prediction for the next 90 days.
%                    1) Get the finals2000A.daily file from the IERS data
%                       center.
%                    2) Get the finals.data file from the US Naval
%                       Observatory. This has EOP parameters since 1992 and
%                       predictions for the next year.
%                    3) Get the finals.data file from the IERS data center.
%      replaceEOPtxt An optional boolean argument. If this parameter is
%                    true, then the file ./EOP.txt will be replaced with
%                    the downloaded data file when the refreshFromSource
%                    parameter is provided and is not -1. replacement only
%                    occurs after the downloaded file is scanned, so if
%                    scanning fails, then EOP.txt will remain unchanged.
%
%OUTPUTS: xpyp A 2XN vector with each column being of the form xpyp=[xp;yp]
%              for the corresponding row in JulUTC1, JulUTC2. These are the
%              polar motion coordinates in radians including the effects of
%              tides and librations. For dates outside of the tabulated
%              values, xpyp=[0;0] will be used, since these values can be
%              difficult to predict forward.
%         dXdY Each column is dXdY=[dX;dY], the celestial pole offsets with
%              respect to the IAU 2006/2000A precession/nutation model in
%              radians, for each date supplied. For dates outside of the
%              tabulated values, dXdY=[0;0] will be used, because these
%              values are difficult to predict and are generally small.
%  deltaUTCUT1 An NX1 vector holding the difference between coordinated
%              universal time (UTC) time and UT1 in seconds for each date
%              provided.
%   deltaTTUT1 An NX1 vector holding the difference between terrestrial
%              time and UT1 in seconds for each date provided. This uses
%              the function cumLeapSec for computing the leap seconds.
%              Thus, the tables within the SOFA library for leapseconds
%              that is used by cumLeapSec must be up to date.
%          LOD An NX1 vector holding the difference between the length of
%              the day using terrestrial time, international atomic time,
%              or UTC without leap seconds and the length of the day in
%              UT1. This is an instantaneous parameter (a derivative). The
%              units are seconds.
%
%Once loaded, the data stays in memory until this function is cleared.
%
%x and y, sometimes called px, py or PMx, PMy are the polar motion
%coordinates and do not include tidal or libration effects. dX and dY are
%the celestial pole offsets with respect to the IAU 2006/2000A precession/
%nutation model to account for free-core nutation, which is poorly modeled.
%LOD is the offset in the length of the day in UT1 (in seconds in a uniform
%scale of TT, TAI, or UTC without leap seconds.
%
%The tidal and libration effects are added using the subroutines
%PMUT1_OCEANS and PM_GRAVI that are originally from the function interp.f
%that are part of the IERS 2010 Conventions and are mentioned in Section
%5.1. The original Fortran routine is available at
%ftp://hpiers.obspm.fr/eop-pc/models/interp.f
%As documented in the "License for IERS Code" file in the
%3rd_Party_Libraries folder, these routines are in the public domain. As
%described in the PMUT1_OCEANS function below, some coefficients from the
%interp.f file have been modified to agree with that is in the IERS 2010
%conventions.
%
%The citation for the IERS 2010 Conventions is [1].
%
%This function is most commonly called as
%[xpyp,dXdY,deltaUTCUT1,deltaTTUT1,LOD]=getEOP(JulUTC1,JulUTC2)
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    persistent dataFile
    
    %Make the inputs into column-vectors.
    JulUTC1=JulUTC1(:);
    JulUTC2=JulUTC2(:);
    
    if(nargin<3)
        refreshFromSource=-1;
    end
    
    if(refreshFromSource~=-1)
       warning('Downloading data from server. There is no need to download again until this function is cleared.')
    end
    
    if(nargin<4)
        replaceEOPtxt=false;
    end
    
    if(isempty(dataFile)||nargin>2)
        switch(refreshFromSource)
            case -1
                %Load the database of Earth Orientation Parameters from the
                %local file.
                ScriptPath=mfilename('fullpath');
                ScriptFolder = fileparts(ScriptPath);
                
                [dataFile,rawText]=processEOPData([ScriptFolder,'/data/EOP.txt']);
            otherwise
                [dataFile,rawText]=processEOPData(refreshFromSource);
        end
    end
    %The coefficient to convert arcseconds to radians.
    as2Rad=(1/60)*(1/60)*pi/180;
    
    JulTable=dataFile(:,4);%Modified Julian Date UTC Julian dates for the tabulated data.
    xpTable=dataFile(:,5);
    ypTable=dataFile(:,6);
    UT1UTCTable=dataFile(:,7);
    LODTable=dataFile(:,8);
    dXTable=dataFile(:,9);
    dYTable=dataFile(:,10);
    
    %Convert the provided two-part date to a modified Julian date in one
    %part for searching and interpolation.
    JulDes=(JulUTC1-2400000.5)+JulUTC2;

    %Interpolate to the date in JulDes and add in tidal and libration effects.
    [xpInt,ypInt,UT1UTCInt,dXInt,dYInt,LODInt]=interpEOP(JulTable,xpTable,ypTable,dXTable,dYTable,UT1UTCTable,LODTable,JulDes);

    %Convert the units of the parameters to return
    xpyp=[xpInt';ypInt']*as2Rad;
    deltaUTCUT1=-UT1UTCInt;
    dXdY=[dXInt';dYInt']*as2Rad;
    LOD=LODInt;

    [year,month,day,dayFrac]=UTC2Cal(JulUTC1,JulUTC2,true);
    
    leapSeconds=cumLeapSec(year,month,day,dayFrac);
    %The 32.184 is the offset of the zero mark of TT versus UTC and UT1.
    deltaTTUT1=deltaUTCUT1+32.184+leapSeconds;
    
    %Overwrite the EOP.txt file, if requested.
    if(replaceEOPtxt==true&&refreshFromSource~=-1)
        ScriptPath=mfilename('fullpath');
        ScriptFolder = fileparts(ScriptPath);
        fileID = fopen([ScriptFolder,'/data/EOP.txt'],'w');
    
        if(fileID==-1)
           error('Could not open EOP.txt to overwrite it.');
        end
        
        fprintf(fileID,'%s',rawText);
        fclose(fileID);
    end
end

function [xpInt,ypInt,UT1UTCInt,dXInt,dYInt,LODInt]=interpEOP(JulTable,xpTable,ypTable,dXTable,dYTable,UT1UTCTable,LODTable,JulDes)
%%INTERPEOP Given tables of Earth orientation parameters, interpolate the
%           values at a particular date that may or may not be withing the
%           tables. This function uses Matlab's interp1 function to
%           perform the interpolation and calls subroutines translated from
%           interp.f in the IERS 2010 conventions to adjust the values for
%           tidal and librational effects.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C

	numDates=size(JulDes,1);

    xpInt=interp1(JulTable,xpTable,JulDes,'pchip',NaN);
    ypInt=interp1(JulTable,ypTable,JulDes,'pchip',NaN);
    UT1UTCInt=interp1(JulTable,UT1UTCTable,JulDes,'pchip',NaN);
    LODInt=interp1(JulTable,LODTable,JulDes,'pchip');
    
    dXInt=interp1(JulTable,dXTable,JulDes,'pchip',NaN);
    dYInt=interp1(JulTable,dYTable,JulDes,'pchip',NaN);

    for curDate=1:numDates
        %Correct for the oceanic effect      
        [cor_x,cor_y,cor_ut1,cor_lod]=PMUT1_OCEANS(JulDes(curDate));

        xpInt(curDate)=xpInt(curDate)+cor_x;
        ypInt(curDate)=ypInt(curDate)+cor_y;
        UT1UTCInt(curDate)=UT1UTCInt(curDate) + cor_ut1;
        LODInt(curDate)=LODInt(curDate)+cor_lod;

        %Correct for the lunisolar effect 
        [cor_x,cor_y]=PM_GRAVI(JulDes(curDate));

        xpInt(curDate)=xpInt(curDate) + cor_x;
        ypInt(curDate)=ypInt(curDate) + cor_y;
    end
    
    %Set values for dates outside of the tabulated range to zero, except
    %for UT1.
    sel=~isfinite(xpInt);
    xpInt(sel)=0;
    ypInt(sel)=0;
    dXInt(sel)=0;
    dYInt(sel)=0;
    
    %Outside of the valid range, use simple linear interpolation for UT1
    UT1UTCInt(sel)=interp1(JulTable,UT1UTCTable,JulDes(sel),'linear','extrap');
end

function [cor_x,cor_y,cor_ut1,cor_lod]=PMUT1_OCEANS(rjd)
%This function is taken from interp.f, which is mentioned in the 
%IERS 2010 conventions, and has been converted to Matlab.
%
%The origin of the lunisolar arguments are from [1] and is what was in
%interp.f. However, the iers conventions modified some of the terms.
%For example, the constant -6962890.2665 in the fundamental argument Omega
%in [1] has been changed to -6962890.5431. Additionally, a coefficient of 
%-0.001037 in interp.f disagrees with the value in [1] and in [2]. Thus,
%the code below has those coefficients modified compared to the original
%ones in interp.f. Comments for those parameters have been added and a
%comment in French has been replaced by one in English.
%
%A link to the original interp.f function is on the site:
%http://hpiers.obspm.fr/eop-pc/models/models_ru.html
%
%REFERENCES:
%[1] J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touzé, G. Francou,
%    and J. Laskar, "Numerical expressions for precession formulae and mean
%    elements for the moon and the planets," Astronomy and Astrophysics,
%    vol. 282, no. 2, pp. 663?683, 1994.
%[2] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%
%----------------------------------------------------------------
%     SUBROUTINE PMUT1_OCEANS (rjd,cor_x,cor_y,cor_ut1,cor_lod)
%
%    This subroutine provides, in time domain, the diurnal/subdiurnal
%    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
%    listed in the program above, have been extracted from the procedure   
%    ortho_eop.f coed by Eanes in 1997.
%    
%    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
%
%    These corrections should be added to "average"
%    EOP values to get estimates of the instantaneous values.
%
%     PARAMETERS ARE :
%     rjd      - epoch of interest given in mjd
%     cor_x    - tidal correction in x (sec. of arc)
%     cor_y    - tidal correction in y (sec. of arc)
%     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
%     cor_lod  - tidal correction in length of day (sec. of time)
%
%     coded by Ch. Bizouard (2002), initially coded by McCarthy and 
%     D.Gambis(1997) for the 8 prominent tidal waves.  
      
      ARG=zeros(6,1); %Array of the tidal arguments   
      DARG=zeros(6,1); %Array of their time derivative 
            
      halfpi = pi/2;
      secrad=2*halfpi/(180*3600);	

%  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
%  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
     data=[1,-1, 0,-2,-2,-2,  -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078;
     1,-2, 0,-2, 0,-1,   0.06,   0.64,  -0.64,   0.06,  0.195, -0.059;
     1,-2, 0,-2, 0,-2,   0.30,   3.42,  -3.42,   0.30,  1.034, -0.314;
     1, 0, 0,-2,-2,-1,   0.08,   0.78,  -0.78,   0.08,  0.224, -0.073;
     1, 0, 0,-2,-2,-2,   0.46,   4.15,  -4.15,   0.45,  1.187, -0.387;
     1,-1, 0,-2, 0,-1,   1.19,   4.96,  -4.96,   1.19,  0.966, -0.474;
     1,-1, 0,-2, 0,-2,   6.24,  26.31, -26.31,   6.23,  5.118, -2.499;
     1, 1, 0,-2,-2,-1,   0.24,   0.94,  -0.94,   0.24,  0.172, -0.090;
     1, 1, 0,-2,-2,-2,   1.28,   4.99,  -4.99,   1.28,  0.911, -0.475;
     1, 0, 0,-2, 0, 0,  -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070;
     1, 0, 0,-2, 0,-1,   9.22,  25.06, -25.06,   9.22,  3.025, -2.280;
     1, 0, 0,-2, 0,-2,  48.82, 132.91,-132.90,  48.82, 16.020,-12.069;
     1,-2, 0, 0, 0, 0,  -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078;
     1, 0, 0, 0,-2, 0,  -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154;
     1,-1, 0,-2, 2,-2,  -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074;
     1, 1, 0,-2, 0,-1,  -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050;
     1, 1, 0,-2, 0,-2,  -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271;
     1,-1, 0, 0, 0, 0,  -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751;
     1,-1, 0, 0, 0,-1,  -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151;
     1, 1, 0, 0,-2, 0,  -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137;
     1, 0,-1,-2, 2,-2,   1.54,   3.03,  -3.03,   1.54,  0.315, -0.189;
     1, 0, 0,-2, 2,-1,  -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035;
     1, 0, 0,-2, 2,-2,  26.13,  51.25, -51.25,  26.13,  5.512, -3.095;
     1, 0, 1,-2, 2,-2,  -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025;
     1, 0,-1, 0, 0, 0,  -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070;
     1, 0, 0, 0, 0, 1,   1.54,   3.00,  -3.00,   1.54,  0.348, -0.171;
     1, 0, 0, 0, 0, 0, -77.48,-151.74, 151.74, -77.48,-17.620,  8.548;
     1, 0, 0, 0, 0,-1, -10.52, -20.56,  20.56, -10.52, -2.392,  1.159;
     1, 0, 0, 0, 0,-2,   0.23,   0.44,  -0.44,   0.23,  0.052, -0.025;
     1, 0, 1, 0, 0, 0,  -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065;
     1, 0, 0, 2,-2, 2,  -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111;
     1,-1, 0, 0, 2, 0,  -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043;
     1, 1, 0, 0, 0, 0,  -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187;
     1, 1, 0, 0, 0,-1,  -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037;
     1, 0, 0, 0, 2, 0,  -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005;
     1, 2, 0, 0, 0, 0,  -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005;
     1, 0, 0, 2, 0, 2,  -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037;
     1, 0, 0, 2, 0, 1,  -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023;
     1, 0, 0, 2, 0, 0,  -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005;
     1, 1, 0, 2, 0, 2,  -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024;
     1, 1, 0, 2, 0, 1,  -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015;
     2,-3, 0,-2, 0,-2,  -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011;
     2,-1, 0,-2,-2,-2,  -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032;
     2,-2, 0,-2, 0,-2,  -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177;
     2, 0, 0,-2,-2,-2,  -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222;
     2, 0, 1,-2,-2,-2,  -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015;
     2,-1,-1,-2, 0,-2,   0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013;
     2,-1, 0,-2, 0,-1,   2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058;
     2,-1, 0,-2, 0,-2, -56.87, -12.93,  11.15,  32.88, -3.795, -1.556;
     2,-1, 1,-2, 0,-2,  -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015;
     2, 1, 0,-2,-2,-2, -11.01,  -2.40,   1.89,   6.41, -0.698, -0.298;
     2, 1, 1,-2,-2,-2,  -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014;
     2,-2, 0,-2, 2,-2,   0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022;
     2, 0,-1,-2, 0,-2,   1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025;
     2, 0, 0,-2, 0,-1,  12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266;
     2, 0, 0,-2, 0,-2,-330.15, -26.96,  37.58, 195.92,-16.195, -7.140;
     2, 0, 1,-2, 0,-2,  -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021;
     2,-1, 0,-2, 2,-2,   2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034;
     2, 1, 0,-2, 0,-2,   9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117;
     2,-1, 0, 0, 0, 0,  -2.35,   0.37,   0.47,   1.41, -0.106, -0.029;
     2,-1, 0, 0, 0,-1,  -1.04,   0.17,   0.21,   0.62, -0.047, -0.013;
     2, 0,-1,-2, 2,-2,  -8.51,   3.50,   3.29,   5.11, -0.437, -0.019;
     2, 0, 0,-2, 2,-2,-144.13,  63.56,  59.23,  86.56, -7.547, -0.159;
     2, 0, 1,-2, 2,-2,   1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000;
     2, 0, 0, 0, 0, 1,   0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001;
     2, 0, 0, 0, 0, 0, -38.48,  19.14,  17.72,  23.11, -2.104,  0.041;
     2, 0, 0, 0, 0,-1, -11.44,   5.75,   5.32,   6.87, -0.627,  0.015;
     2, 0, 0, 0, 0,-2,  -1.24,   0.63,   0.58,   0.75, -0.068,  0.002;
     2, 1, 0, 0, 0, 0,  -1.77,   1.79,   1.71,   1.04, -0.146,  0.037;
     2, 1, 0, 0, 0,-1,  -0.77,   0.78,   0.75,   0.45, -0.064,  0.017;
     2, 0, 0, 2, 0, 2,  -0.33,   0.62,   0.65,   0.19, -0.049,  0.018];
 
     nlines=size(data,1);

     narg=data(:,1:6);
     XSIN=data(:,7);
     XCOS=data(:,8);
     YSIN=data(:,9);
     YCOS=data(:,10);
     UTSIN=data(:,11);
     UTCOS=data(:,12);
 
      T = (rjd - 51544.5)/36525; %julian century

% Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
% and their time derivatives.

      ARG(1) = (67310.54841 +(876600*3600 + 8640184.812866)*T +0.093104*T^2 -6.2e-6*T^3)*15.0 + 648000.0;
      ARG(1)= mod(ARG(1),1296000)*secrad; 
   
      DARG(1) = (876600*3600 + 8640184.812866 + 2 * 0.093104 * T - 3 * 6.2-6*T^2)*15;
      DARG(1) = DARG(1)* secrad / 36525.0;   % rad/day

      %(l) The mean anomaly of the Moon.
      ARG(2) = -0.00024470*T^4 + 0.051635*T^3 + 31.8792*T^2+ 1717915923.2178*T + 485868.249036;
      ARG(2) = mod(ARG(2),1296000)*secrad;
      
      DARG(2) = -4*0.00024470*T^3 + 3*0.051635*T^2 + 2*31.8792*T + 1717915923.2178;
      DARG(2) = DARG(2)* secrad / 36525.0;   % rad/day

      %(l') The mean anomaly of the Sun. Note that the sign of the 0.000136
      %coefficient in interp.f has been flipped so that it agrees with [1]
      %and [2].
      ARG(3) = -0.00001149*T^4 + 0.000136*T^3 -  0.5532*T^2+ 129596581.0481*T + 1287104.79305;
      ARG(3) = mod(ARG(3),1296000)*secrad;

      DARG(3) = -4*0.00001149*T^3 + 3*0.000136*T^2 -  2*0.5532*T + 129596581.0481;
      DARG(3) = DARG(3)* secrad / 36525.0;   % rad/day
      
      %(F) This is the mean longitude of the Moon minus Omega, the mean
      %longitude of the ascending node of the Moon.
      ARG(4) = 0.00000417*T^4 - 0.001037*T^3 - 12.7512*T^2+ 1739527262.8478*T + 335779.526232;
      ARG(4) = mod(ARG(4),1296000)*secrad;

      DARG(4) = 4*0.00000417*T^3 - 3*0.001037*T^2 - 2*12.7512*T + 1739527262.8478;
      DARG(4) = DARG(4)* secrad / 36525.0;   % rad/day
    
      %(D) Mean elongation of the Moon from the Sun.
      ARG(5) = -0.00003169*T^4 + 0.006593*T^3 - 6.3706*T^2+ 1602961601.2090*T + 1072260.70369;
      ARG(5) = mod(ARG(5),1296000)*secrad;

      DARG(5) = -4*0.00003169*T^3 + 3*0.006593*T^2- 2*6.3706*T + 1602961601.2090;
      DARG(5) = DARG(5)* secrad / 36525.0;   % rad/day

      %(Omega) Mean longitude of the ascending node of the Moon. Note that
      %the 6962890.2665 coefficient in [1] has been modified from the
      %interp.f file to agree with the IERS 2010 conventions in [2].
      ARG(6) = -0.00005939*T^4 + 0.007702*T^3+ 7.4722*T^2- 6962890.5431*T + 450160.398036;
      ARG(6) = mod(ARG(6),1296000)*secrad;

      DARG(6) = -4*0.00005939*T^3 + 3*0.007702*T^2+ 2*7.4722*T - 6962890.5431;
      DARG(6) = DARG(6)* secrad / 36525.0;   % rad/day

% CORRECTIONS

	cor_x  = 0;
	cor_y  = 0;
	cor_ut1= 0;
	cor_lod= 0;

 	for j=1:nlines
        ag  = 0;
        dag = 0;
        for i=1:6
            ag  = ag  + narg(j,i)*ARG(i);
            dag = dag + narg(j,i)*DARG(i);
        end
        ag=mod(ag,4*halfpi);

        cor_x  = cor_x   + XCOS(j) *cos(ag) + XSIN(j) * sin(ag);
        cor_y  = cor_y   + YCOS(j) *cos(ag) + YSIN(j) * sin(ag);
        cor_ut1= cor_ut1 + UTCOS(j)*cos(ag) + UTSIN(j)* sin(ag);
        cor_lod= cor_lod -(-UTCOS(j) * sin(ag) + UTSIN(j) * cos(ag) ) * dag;

   end
  
   cor_x   = cor_x * 1e-6;%arcsecond (")
   cor_y   = cor_y * 1e-6;%arcsecond (")
   cor_ut1 = cor_ut1 * 1e-6;%second (s)
   cor_lod = cor_lod * 1e-6;%second (s)
end


function [cor_x,cor_y]=PM_GRAVI (rjd)
%This function is taken from interp.f, which is mentioned in the 
%IERS 2010 conventions, and has been converted to Matlab. The same changes
%to the constants in PMUT1_OCEANS have been made to the constants used in
%this function.
%
%As discussed in the comments to PMUT1_OCEANS, constants relating to
%celestial body positions have been modified from interp.f, which is based
%on [1] to match the values in [2].
%
%REFERENCES:
%[1] J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touzé, G. Francou,
%    and J. Laskar, "Numerical expressions for precession formulae and mean
%    elements for the moon and the planets," Astronomy and Astrophysics,
%    vol. 282, no. 2, pp. 663?683, 1994.
%[2] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%
%----------------------------------------------------------------
%     SUBROUTINE PM_GRAVI (rjd,cor_x,cor_y)
%
%    This subroutine provides, in time domain, the diurnal
%    lunisolar effet on polar motion (")
%    
%    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
%
%    These corrections should be added to "average"
%    EOP values to get estimates of the instantaneous values.
%
%     PARAMETERS ARE :
%     rjd      - epoch of interest given in mjd
%     cor_x    - tidal correction in x (sec. of arc)
%     cor_y    - tidal correction in y (sec. of arc)
%
%     coded by Ch. Bizouard (2002)
      
      halfpi = pi/2;
      secrad=2*halfpi/(180*3600);

%  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
%  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
     data=[1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44;
      1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31;
      1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44;
      1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14;
      1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36;
      1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84;
      1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76;
      1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27;
      1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93;
      1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76];
     
     nlines=size(data,1);
     narg=data(:,1:6);
     XSIN=data(:,7);
     XCOS=data(:,8);
     YSIN=data(:,9);
     YCOS=data(:,10);

      T = (rjd - 51544.5)/36525.0;%  julian century

% Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
% and their time derivatives.

      ARG(1) = (67310.54841 +(876600*3600 + 8640184.812866)*T +0.093104*T^2 -6.2e-6*T^3)*15.0 + 648000.0;
      ARG(1)=mod(ARG(1),1296000)*secrad;
   
      %(l) The mean anomaly of the Moon.
      ARG(2) = -0.00024470*T^4 + 0.051635*T^3 + 31.8792*T^2+ 1717915923.2178*T + 485868.249036;
      ARG(2) = mod(ARG(2),1296000)*secrad;
      
      %(l') The mean anomaly of the Sun. Note that the sign of the 0.000136
      %coefficient in interp.f has been flipped so that it agrees with [1]
      %and [2].
      ARG(3) = -0.00001149*T^4 + 0.000136*T^3 -  0.5532*T^2+ 129596581.0481*T + 1287104.79305;
      ARG(3) = mod(ARG(3),1296000)*secrad;

      %(F) This is the mean longitude of the Moon minus Omega, the mean
      %longitude of the ascending node of the Moon.
      ARG(4) = 0.00000417*T^4 - 0.001037*T^3 - 12.7512*T^2+ 1739527262.8478*T + 335779.526232;
      ARG(4) = mod(ARG(4),1296000)*secrad;

      %(D) Mean elongation of the Moon from the Sun.
      ARG(5) = -0.00003169*T^4 + 0.006593*T^3 - 6.3706*T^2 + 1602961601.2090*T + 1072260.70369;
      ARG(5) = mod(ARG(5),1296000)*secrad;

      %(Omega) Mean longitude of the ascending node of the Moon. Note that
      %the 6962890.2665 coefficient in [1] has been modified from the
      %interp.f file to agree with the IERS 2010 conventions in [2].
      ARG(6) = -0.00005939*T^4 + 0.007702*T^3+ 7.4722*T^2- 6962890.5431*T + 450160.398036;
      ARG(6) = mod(ARG(6),1296000)*secrad;


% CORRECTIONS

	cor_x  = 0;
	cor_y  = 0;

 	for j=1:nlines
        ag  = 0;
        for i=1:6
            ag  = ag  + narg(j,i)*ARG(i);
        end

        ag=mod(ag,4*halfpi);

        cor_x =cor_x+XCOS(j)*cos(ag)+XSIN(j)*sin(ag);
        cor_y =cor_y+YCOS(j)*cos(ag)+YSIN(j)*sin(ag);
   end
  
      cor_x = cor_x * 1e-6;   % arcsecond (")
      cor_y = cor_y * 1e-6;   % arcsecond (")
end


function [dataRet,rawText]=processEOPData(source)
%%PROCESSEOPDATA Process a text file containing Earth orientation parameter
%                data. Either a path can be provided to the text file, or
%                the data can be downloaded from an online source (assuming
%                an available internet connection). The sources provide
%                parameter data for the IERS 2000 models of precession and
%                nutation.
%
%INPUTS: source A number indicating the type and source of the data. This
%               can either be a character string containing the path to a
%               text file to load or it can be an integer indicating
%               whence the data should be obtained. If a text file is
%               provided, it must have the same formatting as one of the
%               files that can be downloaded. Possible integer values
%               indicating online sources are
%               0) (The default) Get the finals2000A.daily file from the
%                  US Naval Observatory. This has EOP parameters for the
%                  last 90 days and prediction for the next 90 days.
%               1) Get the finals2000A.daily file from the IERS data
%                  center.
%               2) Get the finals.data file from the US Naval Observatory.
%                  This has EOP parameters since 1992 and predictions for
%                  the next year.
%               3) Get the finals.data file from the IERS data center.
%
%OUTPUTS: dataRet A matrix where each row is a set of EOP. It has the
%                 format:
%[Year(UTC), month(UTC), day(UTC), Modified Julian Day (UTC),
%x (arcseconds), y (arcseconds), UT1-UTC (seconds), LOD (milliseconds),
%dX (arcseconds), dY (arcseconds)]
%         rawText The raw ASCII text of the downloaded data file. In the
%                 raw file, dX and dY are given in milliarcseconds.
%
%The values of dX and dY that are those where the free-core nutation has
%not been removed. The Err parameters are the errors in the other
%parameters having the same units as the parameters. It is not documented
%whether the errors are one, two or more standard deviations.
%
%Bulletin B values, which are more accurate, are used if available.
%Otherwise, Bulletin A values are used. The difference between Bulletin A
%and B values is given at
%http://hpiers.obspm.fr/iers/bul/bulb/explanatory.html
%
%Some values are not always available, especially as predictions. When a
%value is not available, it is set to zero and the error value associated
%with it is set to infinity. Polar motion coordinates and UT1-UTC should
%always be available.
%
%Documentation for the format of the raw files are given at
%http://maia.usno.navy.mil/ser7/readme.finals2000A
%http://datacenter.iers.org/eop/-/somos/5Rgv/getMeta/10/finals2000A.data
%http://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%http://hpiers.obspm.fr/eoppc/bul/bulb/explanatory.html
%
%When interpolating to a given date, tidal effects need to be added back to
%the values as described in the IERS conventions in [1].
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<1)
   source=0; 
end

%If the data should be loaded from a file.
if(ischar(source))
    fileID = fopen(source);
    
    if(fileID==-1)
       error('Could not open the file at the given path');
    end
    
    rawText = fread(fileID,'*char')';
    fclose(fileID);
else
    options=weboptions('ContentType','text');
    switch(source)
        case 0
            rawText=webread('http://maia.usno.navy.mil/ser7/finals.daily',options);
        case 1
            rawText=webread('https://datacenter.iers.org/eop/-/somos/5Rgv/latest/13',options);
        case 2
            rawText=webread('http://maia.usno.navy.mil/ser7/finals.data',options);
        case 3
            rawText=webread('https://datacenter.iers.org/eop/-/somos/5Rgv/latest/10',options);
        otherwise
            error('Invalid data source given')
    end
end

%Find the indices of all of the carriage returns in the file. These will
%determine when each line ends and the number of them determine the number
%of entries in the file.
returnList=strfind(rawText,10);%10 is the value of the newline character
                                %used in the file. It is not followed or
                                %preceded by a carriage return (13).
numEntry=length(returnList);
dataRet=zeros(numEntry,10);

baseIdx=0;
for curEntry=1:numEntry
    idx=baseIdx+(1:2);
    dataRet(curEntry,1)=sscanf(rawText(idx),'%f');%year
    
    idx=baseIdx+(3:4);
    dataRet(curEntry,2)=sscanf(rawText(idx),'%f');%month
    
    idx=baseIdx+(5:6);
    dataRet(curEntry,3)=sscanf(rawText(idx),'%f');%day
    
    idx=baseIdx+(8:15);%The fractional modified Julian date
    dataRet(curEntry,4)=sscanf(rawText(idx),'%f');
    
    %If the bulletin B data is provided, then use that instead of the
    %Bulletin A data. However, in files where Bulletin B is provided, when
    %it stops being provided, it is replaced with spaces, so we have to
    %test for the length of the line and then also test to make sure that
    %there are not just spaces at the end.
    if(returnList(curEntry)>baseIdx+135)
        idx=baseIdx+(135:144);%x-coordinate of polar motion in arcseconds
        val=sscanf(rawText(idx),'%f');
    else
       val=[]; 
    end
    if(~isempty(val))
        dataRet(curEntry,5)=val;

        idx=baseIdx+(145:154);%y-coordinate of polar motion in arcseconds
        dataRet(curEntry,6)=sscanf(rawText(idx),'%f');
        
        idx=baseIdx+(155:165);
        dataRet(curEntry,7)=sscanf(rawText(idx),'%f');%UT1-UTC
        
        idx=baseIdx+(80:86);
        val=sscanf(rawText(idx),'%f');%LOD, not always filled.
        if(~isempty(val))
            dataRet(curEntry,8)=val;
        else
            %Mark as unavailable by setting it to Inf.
            dataRet(curEntry,8)=Inf;
        end

        idx=baseIdx+(166:175);
        %dX with free core nutation not removed, not always filled.
        %The 1000 converts milliarcseconds to arcseconds.
        dataRet(curEntry,9)=sscanf(rawText(idx),'%f')/1000;
        
        idx=baseIdx+(176:185);
        %dY with free core nutation not removed.
        %The 1000 converts milliarcseconds to arcseconds.
        dataRet(curEntry,10)=sscanf(rawText(idx),'%f')/1000;

        %Go to the next entry.
        baseIdx=returnList(curEntry);
        continue;
    end

    idx=baseIdx+(19:27);%x-coordinate of polar motion in arcseconds
    
    %In finals.data, some of the final rows are just dates with no data. If
    %this is not filled, then we have reached the end of the data.
    val=sscanf(rawText(idx),'%f');
    if(~isempty(val))
        dataRet(curEntry,5)=val;
    else
        dataRet=dataRet(1:curEntry,:);
        break;
    end

    idx=baseIdx+(38:46);%y-coordinate of polar motion in arcseconds
    dataRet(curEntry,6)=sscanf(rawText(idx),'%f');

    idx=baseIdx+(59:68);
    dataRet(curEntry,7)=sscanf(rawText(idx),'%f');%UT1-UTC

    if(returnList(curEntry)>baseIdx+80)
        idx=baseIdx+(80:86);
        val=sscanf(rawText(idx),'%f');%LOD, not always filled.
        if(~isempty(val))
            dataRet(curEntry,8)=val;
        end
        
        if(returnList(curEntry)>baseIdx+98)
            idx=baseIdx+(98:106);
            %dX with free core nutation not removed, not always filled.
            %The 1000 converts milliarcseconds to arcseconds.
            val=sscanf(rawText(idx),'%f')/1000;
            if(~isempty(val))
                dataRet(curEntry,9)=val;

                idx=baseIdx+(117:125);
                %dY with free core nutation not removed.
                %The 1000 converts milliarcseconds to arcseconds.
                dataRet(curEntry,10)=sscanf(rawText(idx),'%f')/1000;
            end
        else%Mark the data as unavailable by setting the error to Inf.
            dataRet(curEntry,15)=Inf;
            dataRet(curEntry,16)=Inf;
        end
    else%Mark the data as unavailable by setting the error to Inf.
        dataRet(curEntry,14)=Inf;
        dataRet(curEntry,15)=Inf;
        dataRet(curEntry,16)=Inf;
    end
    
    baseIdx=returnList(curEntry);
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
