function [constantNames,constantValues,kSize,ss,ipt,au_km,EarthMoonRatio,DENum,IDString]=readJPLEphemHeader(ephemerisPath)
%%READJPLEPHEMHEADER Read the header information of one of the National
%       Aeronautics and Space Administration (NASA) Jet Propulsion
%       Laboratory's (JPL's) Development Ephemerides (DE). These are
%       typically referred to using DE followed by a number and sometimes a
%       letter, such as DE200 or DE430t. This function only reads
%       ephemerides given in binary format for Linux (but will work under
%       Mac OS X and Windows).
%
%INPUTS: ephemerisPath This is a character string of the path to the
%                      ephemeris file. The formatting is the same as used
%                      by the fopen function. If this parameter is omitted
%                      or an empty matrix is passed, it is assumed that the
%                      DE430t ephemeris is given in a data folder that is
%                      in the same folder as this function and the file is
%                      named linux_p1550p2650.430t.
%
%OUTPUTS: constantNames This is a numConstantsX6 character array such that 
%                       constantNames(i,:) is the name of the ith constant
%                       that is given in the header file. The value of the
%                       ith constant is given in constantValues(i).  The
%                       available constants can vary. The meanings of some
%                       of the names that might be present are given below.
%        constantValues A numConstantsX1 vector containing the values of the
%                       constants specified in constantNames.
%                 kSize The number of values in each ephemeris record in
%                       terms of 4-byte binary words.
%                    ss A 3X1 vector such that
%                       ss(1) is the starting Julian ephemeris date of the
%                       ephemeris file.
%                       ss(2) is the ending Julian ephemeris date of the
%                       ephemeris file.
%                       ss(3) is the length of a file record in Julian
%                       days.
%                   ipt This is a 3X15 matrix whose values describe offsets
%                       in each of the data records to access the Chebyshev
%                       interpolation coefficients of the ephemeris.
%                 au_km The length of an astronomical unit as used in the
%                       ephemerides, given in kilometers.
%        EarthMoonRatio The ratio of the mass of the Earth to the mass of
%                       the Moon.
%                 DENum The number of the JPL Development Ephemerides that
%                       are being used.
%              IDString This is a character string that says the name of
%                       ephemerides used as well as the staring and
%                       stopping Julian epoch dates.
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
%The format of the header was deduced by looking at the Fortran routine
%called testeph1.f that was created by the California Institute of
%Technology (CIT) under a U.S. government contract with NASA and is
%available for download in the folder
%ftp://ssd.jpl.nasa.gov/pub/eph/planets/
%
%The meanings of some of the names that might begiven in constantNames are:
%'DENUM ' The number of the planetary ephemeris loaded. The value of this
%         should equal the DENum output of this function.
%'LENUM ' The number of the lunear ephemeris loaded. This usually equals
%         DENum.
%'CLIGHT' The speed of light in km/s.
%'AU    ' The number of kilometers per astronomical unit (AU). The value of
%         this should equal the value au_km of this function.
%'EMRAT ' The Earth-Moon mass ratio. The value should equal the
%         EarthMoonRatio value returned by this function.
%'JDEPOC' The TDB Julian date epoch at which initial conditions for
%         planets and other objects are given.'
%'CENTER' The reference center for the initial position conditions. This is
%         typically 11 for the Sun  and either 0 or 12 for the solar
%         system barycenter
%'GMB   ' The universal gravitational constant times the mass of the Earth-
%         Moon system (the Earth-Moon barycenter) in AU^3/day^2.
%'GMS   ' The universal gravitational constant times the mass of the Sun in
%         AU^3/day^2.
%'GMi   ' The universal gravitational constant times the ith planet in the
%         solar system going away from the Sun. i=1 denotes Mercury, i=3,
%         the Earth, i=8 Neptune and i=9 Pluto (regardless of whether one
%         really considers it a planet). Values typically do not go over 9.
%'Maiiii' The universal gravitational constant times the mass of the
%         asteroid given by the 4 digit number iiii.
%'Xi', 'Yi', 'Zi' The position components of the ith body in the solar
%         system at the Julian epoch in  (X,Y,Z) Cartesian coordinates
%         taken with respect to the location given by 'CENTER' given in AU.
%'XDi', 'YDi', 'ZDi' The velocity components of the ith body in the solar
%         system at the Julian epoch taken with respect to the location
%         given by 'CENTER' given in astronomical units given in AU/day.
%'XS', 'YS', 'ZS' and 'XDS', 'YDS', 'ZDS' The location of the Sun at the
%         Julian epoch taken with respect to the location given by 'CENTER'
%         given in AU and AU/day.
%'XM', 'YM', 'ZM' and 'XDM', 'YDM', 'ZDM' The location of the Moon taken
%         with respect to the Earth at the Julian epoch given in AU and
%         AU/day.
%'PHI', 'THT', 'PSI' Euler angles of orientation in radians of the lunar
%         mantle at the Julian epoch.
%'PHIC', 'THTC', 'PSIC' Euler angles of orientation in radians of the lunar
%         core at the Julian epoch.
%'RE    ' The radius of the Earth in kilometers.
%'AM    ' The radius of the Sun in kilometers.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If no ephemeris is given, assume that the DE430t is in a folder named data
%that is located in the same folder as this file.
if(nargin<1||isempty(ephemerisPath)) 
%Get the path to this file.
    ScriptPath=mfilename('fullpath');
    ScriptFolder=fileparts(ScriptPath);
    ephemerisPath=[ScriptFolder,'/data/linux_p1550p2650.430t'];
end

%The maximum number of constants given in older ephemeris files. Newer
%ephemeris files can have more. We will know what type of file we have from
%the numConst value (below). 
oldNumConstLimit=400;

fid=fopen(ephemerisPath,'r');

%ttl holds a text string with the name of the ephemerides, the start epoch
%and the final epoch as strings.
IDString=fread(fid, 252,'*char').';
constantNames=fread(fid, 6*oldNumConstLimit,'*char');
%ss(1) is the starting Julian ephemeris date of the ephemeris file.
%ss(2) is the ending Julian ephemeris date of the ephemeris file.
%ss(3) is the length in days of a file record.
ss=fread(fid,3,'double');
%The total number of names of constants in the ephemeris file.
numConst=fread(fid,1,'int');
%The number of kilometers per astronomical unit
au_km=fread(fid,1,'double');
%Earth-to-Moon mass ratio
EarthMoonRatio=fread(fid,1,'double');

ipt=zeros(3,15);
ipt(:,1:12)=fread(fid,[3,12],'int');
%This is the planetary ephemeris number. For example, 430 for the DE430
%ephemerides.
DENum=fread(fid,1,'int');
ipt(:,13)=fread(fid,3,'int');

if(numConst<=oldNumConstLimit)
    %If one is using an older ephemeris file that can contain at most
    %oldNumConstLimit constants.
    constantNames=reshape(constantNames(1:6*numConst),[6,numConst]).';
    ipt(:,14)=fread(fid,3,'int');
    ipt(:,15)=fread(fid,3,'int');
else
    %If one is using a newer format file with more than oldNumConstLimit
    %constants, then a bunch of additional constant names are given before
    %the final two columns of the ipt matrix.
    constantNames=reshape([constantNames;fread(fid,numConst*6-6*oldNumConstLimit,'*char')],[6,numConst]).';
    ipt(:,14)=fread(fid,3,'int');
    ipt(:,15)=fread(fid,3,'int');
end

%Scan ahead to where the first constant starts. We assume that the first
%constant is not zero.
lastVal=0;
while(lastVal==0)
    lastVal=fread(fid,1,'double');
end

%This is the length of one record in this file in terms of binary words.
kSize=(ftell(fid)-8)/4;

%Read in all the constants that are associated with the names in cnam.
constantValues=[lastVal;fread(fid,numConst-1,'double')];

fclose(fid);
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
