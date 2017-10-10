function [deltaC,deltaS]=gravOceanTideOffset(TT1,TT2)
%%GRAVOCEANTIDEOFFSET Compute the offsets to add the fully normalized
%                     spherical harmonic gravitational coefficients to
%                     handle the effects of the ocean Earth tides. Solid
%                     Earth tides and other effects need to be added
%                     separately to get the full offset due to tides.
%
%%INPUTS TT1,TT2 Two parts of a Julian date given in terrestrial time (TT). 
%                The units of the date are days. The full date is the sum
%                of both terms. The date is broken into two parts to
%                provide more bits of precision. It does not matter how the
%                date is split.
%
%OUTPUTS:deltaC An array class holding the offsets to the fully normalized
%               coefficient terms that are multiplied by cosines in the
%               spherical harmonic expansion of the gravitational
%               potential. If given to a CountingClusterSet class, then
%               C(n+1,m+1) is the coefficient of degree n and order m. The
%               offsets only go up to the degree and order of the FES2004
%               coefficients.
%        deltaS An array holding the coefficient terms that are multiplied
%               by sines in the spherical harmonic expansion. The format of
%               of S is the same as that of C.
%
%This implements the FES2004 tide model that is described in [1]. The
%FES2004 coefficients are loaded from the file fes2004 Cnm-Snm.dat that
%was downloaded from
%http://maia.usno.navy.mil/conv2010/convupdt/chapter6/tidemodels/fes2004_Cnm-Snm.dat
%and placed in the ./data folder.
%
%To get the full offset in the gravitational coefficients due to tides,
%the functions gravSolidTideOffset, gravOceanTideOffset and
%gravPoleTideOffset have to be combined.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Read the FES2004 coefficient file, which should be located in a data
%folder that is in the same folder as this file.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

fileID=fopen([ScriptFolder,'/data/fes2004_Cnm-Snm.dat']);
data=textscan(fileID,'%s','CommentStyle','#','whitespace',' ','delimiter','\n');
fclose(fileID);
data=data{1};

%Put all of the elements for each row into a cell array.
numRows=length(data);
rowData=cell(numRows-4,1);
for curRow=5:numRows
    V=textscan(data{curRow},'%s');
    rowData{curRow-4}=V{1};
end

numRows=numRows-4;%The number of data rows.
%Allocate a matrix to hold the data values of the table. The first six
%values are the Doodson parameters.
dataMat=zeros(numRows,12);
for curRow=1:numRows
    %Extract the digits of the Doodson parameters. The Doodson parameters
    %have all but the first digit augmented by 5.
    VScan=textscan(rowData{curRow}{1},'%f');
    DoodsonVal=VScan{1};
    %Extract the digits. The -5 gets rid of the offset of 5 on all but the
    %first digit.
    dataMat(curRow,1)=fix(mod(fix(DoodsonVal/100),10));
    dataMat(curRow,2)=fix(mod(fix(DoodsonVal/10),10))-5;
    dataMat(curRow,3)=fix(mod(fix(DoodsonVal),10))-5;
    dataMat(curRow,4)=fix(mod(fix(DoodsonVal*10),10))-5;
    dataMat(curRow,5)=fix(mod(fix(DoodsonVal*100),10))-5;
    dataMat(curRow,6)=fix(mod(fix(DoodsonVal*1000),10))-5;
    
    for n=7:12
        %The -4 skips the Darw field in the data file.
        VScan=textscan(rowData{curRow}{n-4},'%f');
        dataMat(curRow,n)=VScan{1};
    end    
end

%The units of the coefficients are 10^(-12).
dataMat(:,9:12)=dataMat(:,9:12)*10^(-12);

%%NEXT, SYNTHESIZE THE COEFFICIENTS USING EQUATION 6.15
%The maximum order of the coefficients that will be offset.
M=max(dataMat(:,7));
totalNumCoeffs=(M+1)*(M+2)/2;
emptyData=zeros(totalNumCoeffs,1);
deltaC=CountingClusterSet(emptyData);
deltaS=CountingClusterSet(emptyData);

%The next few things before the loop are needed to compute the thetaF
%values.

%Convert the terrestrial time to Greenwich mean sidereal time in radians.
GMST=TT2GMST(TT1,TT2);

%Turn terrestrial time into barycentric dynamical time TDB.
[TDB1,TDB2]=TT2TDB(TT1,TT2);

%Get the Delaunay variables in radians
vec=DelaunayVar(TDB1,TDB2);
l=vec(1);
lp=vec(2);
F=vec(3);
D=vec(4);
Omega=vec(5);

s=F+Omega;
beta(2)=s;%s: Moon's mean longitude
beta(1)=GMST+pi-s;%tau: time angle in lunar days reckoned from lower transit
beta(3)=s-D;%h: Sun's mean longitude
beta(4)=s-l;%p: longitude of the Moon's mean perigee
beta(5)=-Omega;%N':negative longitude of the Moon's mean node
beta(6)=s-D-lp;%pl: Longitude of the Sun's mean perigee.

%Evaluate Equation 6.15
for curData=1:numRows
    thetaF=sum(beta.*dataMat(curData,1:6));
   
    n=dataMat(curData,7);
    m=dataMat(curData,8);
    
    sumVal=0;
    %First, the upper signs in Equation 6.15
    CPlusnm=dataMat(curData,9);
    SPLusnm=dataMat(curData,10);
    sumVal=sumVal+(CPlusnm-1i*SPLusnm)*exp(1i*thetaF);
    
    %Then the lower signs
    CMinusnm=dataMat(curData,11);
    SMinusnm=dataMat(curData,12);
    sumVal=sumVal+(CMinusnm+1i*SMinusnm)*exp(-1i*thetaF);
    
    deltaC(n+1,m+1)=deltaC(n+1,m+1)+real(sumVal);
    deltaS(n+1,m+1)=deltaS(n+1,m+1)-imag(sumVal);
end

deltaC=deltaC.clusterEls;
deltaS=deltaS.clusterEls;

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
