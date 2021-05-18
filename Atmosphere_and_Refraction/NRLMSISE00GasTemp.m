function [gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A)
%%NRLMSISE00GASTEMP  Obtain the gas constituents and temperature
%                    parameters of the atmosphere at a particular location
%                    and time using the NRLMSISE-00 atmospheric model for
%                    dry air.
%
%INPUTS: dayOfYear The integer day of the year in the Gregorian calendar
%                  in universal coordinated time (UTC). Counting starts at
%                  1. the resolution of the model is not sufficient for it
%                  to matter whether 365 or 366 is given at the day if it
%                  isn't/is a leap year.
%   secondOfTheDay The second of the day. This starts at zero. The
%                  resolution of the model is not high enough for leap
%                  seconds to matter, so values above 86400.0 are just
%                  clipped to 86400.0.
%        latLonAlt The location under consideration given in WGS-84
%                  ellipsoidal coordinates of latitude and longitude in
%                  radians and ellipsoidal height in meters.
%               Ap An optional parameter specifying the average daily
%                  geomagnetic index. If omitted, the value 4.0, which is
%                  suitable for models below 90km is used. The index at a
%                  particular time can be obtained from a service of the
%                  International Service of Geomagnetic Indices at
%                  http://www-app3.gfz-potsdam.de/kp_index/
%                  A scalar value of AP indicates standard mode. If a 7X1
%                  vector is passed, then it is assumed that  storm mode
%                  is desired. In storm mode, Ap(1) is the daily Ap index.
%                  Ap(2) through Ap(5) are 3 hours Ap indices respectively
%                  for the current time and 3, 6, and 9 hours before the
%                  current time. Ap(6) is the average of 8 3-hour Ap
%                  indices from 12 to 33 hours prior to the current time
%                  and Ap(7) is the average of 8 3-hour Ap indices from 36
%                  to 59 hours prior to the current time.
%             F107 An optional parameter specifying the 10.7cm solar radio
%                  flux for the previous day. If omitted, a value of 150
%                  is used, which should be suitable for models below
%                  90km. Values of the index can be obtained from
%                  http://www.swpc.noaa.gov/ftpdir/indices/quar_DSD.txt
%            F107A An optional parameter specifying the 81 day average of
%                  the 10.7cm radio solar flux. If omitted, a value of 150
%                  is used, which should be suitable for models below
%                  90km.
%
%OUTPUTS: gasTable A cell array of constituent elements of the atmosphere
%                  and their number densities. gasTable{i,1} is a string
%                  listing the name of the ith constituent element of the
%                  atmosphere. This can take the values
%                  'He' Helium
%                  'O'  Elemental Oxygen
%                  'N2' Nitrogen
%                  'O2' Oxygen
%                  'Ar' Argon
%                  'H'  Elemental Hydrogen
%                  'N'  Elemental Nitrogen
%                  'O*' Anomalous Oxygen
%                  gasTable{i,2} is the corresponding number density of
%                  the ith constituent element in units of atoms per cubic
%                  meter.
%                t Temperature variables. t(1) is the exospheric
%                  temperature and t(2) is the temperature at altitude.
%                  The temperature is in units of Kelvin.
%                d A 9X1 vector containing the same information as
%                  gasTable but without labels. The entries in d are
%                  d(1) - HE NUMBER DENSITY
%                  d(2) - O NUMBER DENSITY
%                  d(3) - N2 NUMBER DENSITY
%                  d(4) - O2 NUMBER DENSITY
%                  d(5) - AR NUMBER DENSITY                       
%                  d(6) - TOTAL MASS DENSITY  without d(9)
%                  d(7) - H NUMBER DENSITY
%                  d(8) - N NUMBER DENSITY
%                  d(9) - Anomalous oxygen NUMBER DENSITY(M-3)
%                  where all number densities are in units of atoms per
%                  cubic meter and the units of d(6) are kilograms per
%                  cubic meter. Note that d(6) is computed using a
%                  definition of the atomic mass unit and definitions of
%                  the atomic masses that are not necessarily the most up
%                  to date. 
%
%This function is essentially a wrapper for the C-code implementation of
%the  NRLMSISE-00 atmospheric model using appropriate functions to put the
%time and location in the correct format. The model is described in [1].
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt);
%or
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP);
%or
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107);
%or
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
%
%EXAMPLE 1:
%dayOfYear=172;
%secondOfTheDay=29000;
%latLonAlt=[60*(pi/180);-70*(pi/180);400e3];
%F107=150;
%F107A=150;
%AP=4;
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
%
%EXAMPLE 2
%In most cases, one might be given a Julian Date or a Calendar Date.
%Here, we convert the date appropriately.
%year=2012;
%month=11;
%day=16;
%hour=0;
%minute=0;
%second=0;
%F107=150;
%F107A=150;
%AP=4;
%latLonAlt=[60*(pi/180);-70*(pi/180);400e3];
%[Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
%[~,dayOfYear,secondOfTheDay]=UTC2DayCount(Jul1,Jul2);
%[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
%
%REFERENCES:
%[1] J. M. Picone, A. E. Hedin, D. P. Drob, and A. C. Aikin, "NRLMSISE-00
%    empirical model of the atmosphere: Statistical comparisons and
%    scientific issues," Journal of Geophysical Research: Space Physics,
%    vol. 107, no. A12, Dec. 2002.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
