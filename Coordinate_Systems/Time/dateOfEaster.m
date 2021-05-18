function [day,month]=dateOfEaster(year)
%%DATEOFEASTER Determine the month and day of the holiday Easter under the
%              western system (not the orthodox system). This is only valid
%              after 1582, which is the year that the calendar used in much
%              of Europe switched from Julian to Gregorian. Note that
%              different countries made the switch at different times. 
%
%INPUTS: The integer year number given in the Gregorian calendar.
%
%OUTPUTS: day The day of the month, starting at 1, when Easter occurs.
%       month The month when Easter occurs. This is either 3 (March) or 4
%             (April).
%
%This function implements algorithm E of Chapter 1.3.2 of [1].
%
%EXAMPLE:
% [day,month]=dateOfEaster(2017)
%One will see that Easter was April 16th.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Step E1
%The golden number of the year in the 19-year metonic cycle.
G=mod(year,19)+1;

%Step E2
%The century number.
C=floor(year/100)+1;

%Step E3
X=floor(3*C/4)-12;
Z=floor((8*C+5)/25)-5;

%Step E4
D=floor(5*year/4)-X-10;

%Step E5
E=mod(11*G+20+Z-X,30);
if(E==25&&G>1||E==24)
    E=E+1;
end

%Step E6
N=44-E;
if(N<21)
   N=N+30;
end

%Step E7
N=N+7-mod(D+N,7);

%Step E8
if(N>31)
    day=N-31;
    month=4;
else
    day=N;
    month=3;
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
