function [C,S,a,c]=getEMMCoeffs(M,year,fullyNormalize)
%%GETEMMCOEFFS Obtain spherical harmonic coefficients for the 2017 
%              version of the National Oceanic and Atmospheric
%              Administration's (NOAA's) Enchaned Magnetic Model (EMM) at a
%              particular time or at the reference epoch (2017). The model
%              is considered valid down to 10km underground.
%
%INPUTS: M The integer maximum order of the spherical harmonic coefficients
%          obtained. This is a value between 1 and 740. If this parameter
%          is omitted, the default value is 740. If one wishes to load all
%          coefficients, one can also just pass Inf.
%     year A decimal number indicating a year in the Gregorian calendar as
%          specified by universal coordinated time (UTC). For example,
%          halfway through the year 2017 would be represented as 2017.5.
%          The precision of the model is not sufficiently fine that leap
%          seconds matter. If this parameter is omitted, then the
%          reference year of the EMM model is used. In this instance,
%          2017.0.
% fullyNormalize Geomagnetic models are normally given in terms of Schmidt
%          semi-normalized Legendre functions. If fullyNormalize=true, then
%          the coefficients are converted for use with fully normalized
%          associated Legendre functions, as are commonly used with
%          spherical harmonic algorithms for gravitational models. If this
%          parameter is omitted, the default value is true.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the harmonic expansion. This can be given to a
%           CountingClusterSet so that C(n+1,m+1) is the coefficient of
%           degree n and order m. The coefficients have units of Tesla. The
%           coefficients are normalized according to the fullyNormalize
%           term.
%         S An array holding the coefficient terms that are multiplied by
%           sines in the harmonic expansion. The format of S is the same as
%           that of C.
%         a The numerator in the (a/r)^n term in the spherical harmonic sum
%           having units of meters.
%         c The constant value by which the spherical harmonic series is
%           multiplied, having units of squared meters.
%
%Details on the normalization of the coefficients is given in the comments
%to the function spherHarmonicEval.
%
%This function first checks for a .mat file with the coefficients in it.
%The .mat file is small and can be read quite quickly. However, if one does
%not exist, then it tries to read the EMM2015.COF and EMM2015SV.COF text
%files that one can obtain directly from the NOAA. Reading from the text
%files is very slow.
%
%Documentation is in [1] and [2]. Documentation as well as some code
%containing containing the coefficients is available at
%http://www.ngdc.noaa.gov/geomag/EMM/
%Being created by the U.S. government, the model is not subject to
%copyright.
%
%The data is kept in zipped files in the ./data folder. If all of the data
%is being loaded for a particular year (M=Inf), then a .mat file will be
%created in the ./data folder with the data for that year so that it can be
%loaded more quickly in the future.
%
%REFERENCES:
%[1] S. Maus, "An ellipsoidal harmonic representation of Earth's
%    lithospheric magnetic field to degree and order 720," Geochemistry,
%    Geophysics, Geosystems, vol. 11, no. 6, Jun. 2010.
%[2] Information on Schmidt semi-normalized Legendre functions is given in
%    D. E. Winch, D. J. Ivers, J. P. R. Turner, and R. J. Stening,
%    "Geomagnetism and Schmidt quasi-normalization," Geophysical Journal
%    International, vol. 160, no. 2, pp. 487-504, Feb. 2005.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(M))
    M=Inf;%Load all of the coefficients.
end

yearRef=2015.0;
if(nargin<2||isempty(year))
    year=yearRef;
end

if(nargin<3||isempty(fullyNormalize))
    fullyNormalize=true;
end

saveMatFile=(M==Inf);

%The EMM2015 magnetic coefficient data file, should be located in a data 
%folder that is in the same folder as this file. This find the path to this
%file.
ScriptPath=mfilename('fullpath');
ScriptFolder=fileparts(ScriptPath);

%We have to load the version of the EMM that predicts for the year
%given. To know which one to load, we will see which years are in the
%zip archive. The years are given in the filenames.
filePathAndName=fileNamesInZipArchive([ScriptFolder,'/data/EMM_Coefficients.zip']);
secVarFilePathAndName=fileNamesInZipArchive([ScriptFolder,'/data/EMM_Secular_Variations.zip']);

numFiles=length(filePathAndName);
EMMYears=zeros(numFiles,1);
for curFile=1:numFiles
    [~,fileName] = fileparts(filePathAndName{curFile});
    %The filenames all begin with EMM, so skip the first three
    %characters to get the years.
    EMMYears(curFile)=sscanf(fileName(4:end),'%f',1);
end

%We have to find the EMM model year that is equal to the integer part
%of year. If there is none, then we have to find the closest year.
yearIdx2Load=find(EMMYears==fix(year),1);
if(isempty(yearIdx2Load)||yearIdx2Load~=year)
    [~,yearIdx2Load]=min(abs(EMMYears-fix(year)));
end
yearRef=EMMYears(yearIdx2Load);

%First, see if a .mat file with all of the data for the appropriate year
%exists. If so, then use that to load the coefficients from which
%interpolation must be performed.
matFile=[ScriptFolder,'/data/EMM',int2str(yearRef),'.mat'];

if(exist(matFile,'file'))
    load(matFile,'CCoeffs','SCoeffs','C1Coeffs','S1Coeffs');

    %Keep only as many coefficients as the maximum order provided.
    totalNumCoeffs=(M+1)*(M+2)/2;
    if(length(CCoeffs)<totalNumCoeffs)
        totalNumCoeffs=length(CCoeffs);
    end
    C=CountingClusterSet(CCoeffs(1:totalNumCoeffs));
    S=CountingClusterSet(SCoeffs(1:totalNumCoeffs));
    M=C.numClust-1;

    %If the drift coefficients need to be reduced in size.
    totalNumDriftCoeffs=min(length(C1Coeffs),totalNumCoeffs);
    if(length(C1Coeffs)<totalNumCoeffs)
        totalNumDriftCoeffs=length(C1Coeffs);
    end
    
    C1=CountingClusterSet(C1Coeffs(1:totalNumDriftCoeffs));
    S1=CountingClusterSet(S1Coeffs(1:totalNumDriftCoeffs));
else
    %Otherwise, just read the data from the text files.
    
    %We now know which coefficients and secular variation terms to load.
    temp=readZipArchive([ScriptFolder,'/data/EMM_Coefficients.zip'],filePathAndName{yearIdx2Load});
    coefficientFile=char(temp{2});
    temp=readZipArchive([ScriptFolder,'/data/EMM_Secular_Variations.zip'],secVarFilePathAndName{yearIdx2Load});
    secVarTerms=char(temp{2});
    
    %The first line of the coefficient file is just a description, so find
    %the first occurence of a return character and only take the file after
    %that character
    startIdx=1;
    %char(10) is the newline character used in the file.
    while(strcmp(coefficientFile(startIdx),newline)~=1)
        startIdx=startIdx+1;
    end
    startIdx=startIdx+1;
    %Load all of the coefficients.
    dataCoeffs=sscanf(coefficientFile(startIdx:end),'%d %d %f %f',Inf);
    numCoeffEntries=length(dataCoeffs)/4;
    dataCoeffs=reshape(dataCoeffs,4,numCoeffEntries);
    
    %The number of coefficients.
    M=min(M,dataCoeffs(1,end));
    totalNumCoeffs=(M+1)*(M+2)/2;
    
    %Load all of the drift coefficients.
    driftCoeffs=sscanf(secVarTerms,'%d %d %f %f',Inf);
    numDriftEntries=length(driftCoeffs)/4;
    driftCoeffs=reshape(driftCoeffs,4,numDriftEntries);
    MDrift=driftCoeffs(1,end);
    totalNumDriftCoeffs=(MDrift+1)*(MDrift+2)/2;
    
    %Allocate space for the coefficients.
    emptyData=zeros(totalNumCoeffs,1);
    C=CountingClusterSet(emptyData);
    S=CountingClusterSet(emptyData);
    %These will hold the time-varying terms.
    emptyData=zeros(totalNumDriftCoeffs,1);
    C1=CountingClusterSet(emptyData);
    S1=CountingClusterSet(emptyData);
    
    %Put the data into the CountingClusterSet class instances. 
    for curRow=1:numCoeffEntries
        n=dataCoeffs(1,curRow);%The degree of the coefficient.
        m=dataCoeffs(2,curRow);%The order of the coefficient.
        C(n+1,m+1)=dataCoeffs(3,curRow);
        S(n+1,m+1)=dataCoeffs(4,curRow);
    end
    
    for curRow=1:numDriftEntries
        n=driftCoeffs(1,curRow);%The degree of the coefficient.
        m=driftCoeffs(2,curRow);%The order of the coefficient.
        C1(n+1,m+1)=driftCoeffs(3,curRow);
        S1(n+1,m+1)=driftCoeffs(4,curRow);
    end

    %Change the units fron Nanotesla to Tesla.
    C(:)=10^(-9)*C(:);
    S(:)=10^(-9)*S(:);
    C1(:)=10^(-9)*C1(:);
    S1(:)=10^(-9)*S1(:);
    
    %If the coefficients should be saved into a .mat file for more
    %efficient loading in the future.
    if(saveMatFile)
        CCoeffs=C.clusterEls;
        SCoeffs=S.clusterEls;

        C1Coeffs=C1.clusterEls;
        S1Coeffs=S1.clusterEls;

        save(matFile,'CCoeffs','SCoeffs','C1Coeffs','S1Coeffs');
    end
end

%If interpolation to other dates must be performed.
if(year~=yearRef)
    yearDiff=year-yearRef;
    %Perform linear interpolation.
    idx=1:totalNumDriftCoeffs;
    C.clusterEls(idx)=C.clusterEls(idx)+yearDiff*C1.clusterEls(idx);
    S.clusterEls(idx)=S.clusterEls(idx)+yearDiff*S1.clusterEls(idx);
end

%If the coefficients should be fully normalized.
if(fullyNormalize~=false)
     for n=0:M
        k=1/sqrt(1+2*n);
        C(n+1,:)=k*C(n+1,:);
        S(n+1,:)=k*S(n+1,:);
     end
end

%Return S and C as arrays, not as a CountingClusterSet classes.
C=C.clusterEls;
S=S.clusterEls;

%The EMM2015 model uses the same reference ellipse as the WMM2010.
a=Constants.WMM2010SphereRad;%meters
c=a^2;
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
