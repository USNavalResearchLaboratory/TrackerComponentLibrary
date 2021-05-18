function [C,S,a,c]=getIGRFCoeffs(year,fullyNormalize)
%%GETIGRFCOEFFS Obtain spherical harmonic coefficients for the
%               13th generation International Geomagnetic Reference Field
%               (IGRF) at a particular time or at the latest reference
%               epoch. Note that the coefficients must be adjusted for
%               local magnetization effects (buildings, cars, power lines,
%               ionospheric currents) in order to hope to obtain something
%               close to what a local measurement might provide.
%
%INPUTS: year A decimal number indicating a year in the Gregorian
%             calendar as specified by universal coordinated time (UTC).
%             For example, halfway through the year 2013 would be
%             represented as 2013.5. The precision of the model is not
%             sufficiently fine that leap seconds matter. If this parameter
%             is omitted or an empty matrix is passed, then the reference
%             year of the IGRF model is used. In this instance, 2020.0. For
%             predicting the magnetic field into the future, it is
%             suggested that one use the World Magnetic Model (WMM) not the
%             IGRF.
% fullyNormalize Geomagnetic models are normally given in terms of Schmidt
%             semi-normalized Legendre functions. If fullyNormalize =true,
%             then the coefficients are converted for use with fully
%             normalized associated Legendre functions, as are commonly
%             used with spherical harmonic algorithms for gravitational
%             models. If this parameter is omitted, the default value is
%             true.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the harmonic expansion. This can be given to a
%           CountingClusterSet class so that C(n+1,m+1) is the coefficient
%           of degree n and order m. The coefficients have units of Tesla.
%           The coefficients are normalized according to the fullyNormalize
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
%Documentation for the IGRF model is given in [1] and information on
%Schmidt semi-normalized Legendre functions is given in [2]. The
%coefficients originate from
%http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%and being created by the government are not subject to copyright
%protection.
%
%Note that the IGRF is not a crustal magnetic model; it is just the main
%field of the Earth. The Enhanced Magnetic Model (EMM) from the National
%Oceanic and Atmospheric Administration's (NOAA's) National Geophysical
%Data Center is a complete magnetic model.
%
%REFERENCES:
%[1] International Association of Geomagnetism and Aeronomy, Working
%    Group V-MOD, "International geomagnetic reference field: Eleventh
%    generation," Geophysical Journal International, vol. 183, no. 3, pp.
%    1216-1230, Dec. 2010.
%[2] D. E. Winch, D. J. Ivers, J. P. R. Turner, and R. J. Stening,
%    "Geomagnetism and Schmidt quasi-normalization," Geophysical Journal
%    International, vol. 160, no. 2, pp. 487-504, Feb. 2005.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Read the IGRF data file, which should be located in a data folder that is
%in the same folder as this file.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

fileID=fopen([ScriptFolder,'/data/igrf13coeffs.txt']);
data=textscan(fileID,'%s','CommentStyle','#','whitespace',' ','delimiter','\n');
fclose(fileID);
data=data{1};
%data{1} is just a bunch of labels for IGRF, DGRF and SV. IGRF and DGRF
%represent the coefficients at particular epochs. SV is the last column and
%provides coefficients for predicting the most recent model forward in
%time.
%data{2} is the row of dates for each of the models, with the last one being a
%range for the validity of the prediction coefficients.
%The rest of the M cells are all of the data. A cell can be labeled either
%g or h, which corresponds to C and S coefficients.

%Put all of the header elements into a cell array.
HeaderData=textscan(data{2},'%s','delimiter',' ');
HeaderData=HeaderData{1};
%Get rid of the empty cells that arise when there are two spaces.
sel=~cellfun('isempty',HeaderData);
HeaderData=HeaderData(sel);

%Put all of the elements for each row into a cell array.
numRows=length(data)-2;
rowData=cell(numRows,1);
for curRow=1:numRows
    V=textscan(data{curRow+2},'%s','delimiter',' ');
    temp=V{1};
    %Get rid of the empty cells that arise when there are two spaces.
    sel=~cellfun('isempty',temp);
    rowData{curRow}=temp(sel);
end

%The last row should have the highest degree and order. Thus, we shall
%extract the maximum degree and order to know the size of the ClusterSet
%classes that should be constructed to hold the coefficients.
temp=textscan(rowData{numRows}{2},'%f');
M=temp{1};

%Allocate space for two sets of coefficients. these are needed so that
%interpolation between the years can be performed, if necessary.
totalNumCoeffs=(M+1)*(M+2)/2;
emptyData=zeros(totalNumCoeffs,1);
C=CountingClusterSet(emptyData);
S=CountingClusterSet(emptyData);
C1=CountingClusterSet(emptyData);
S1=CountingClusterSet(emptyData);

%Next, the closest two years of coefficients for the data must be read in
%so that one can linearly interpolate to the proper time. If the requested
%after the final year, then the SV values must be read in for future
%interpolation. The given Julian date must then be converted to a decimal
%year. If no date is given, then the final reference date in the model is
%used and no interpolation is performed.

if(nargin<1||isempty(year))
    temp=textscan(HeaderData{end-1},'%f');
    year=temp{1};
end

if(nargin<2||isempty(fullyNormalize))
   fullyNormalize=true; 
end

%Next, the years in the table have to be converted to decimal values. The
%first three entries in the header just label the type of coefficient and
%the degree and order. Also, the last entry just marks the valid year range
%for the prediction coefficients.
numYears=length(HeaderData)-4;
yearList=zeros(numYears,1);
for curYear=1:numYears
   temp=textscan(HeaderData{curYear+3},'%f');
   yearList(curYear)=temp{1};
end

%Now, find the indices of the two closest years, unless the value given is
%exactly on a reference year, in which case, we can just use the values of
%the reference year itself or unless the year is after the last reference
%year, in which case, we have to extrapolate past the last reference year.
if(year<yearList(1))
   warning('The year is before the first year of the IGRF models. The values will be extrapolated backwards using the first two years for which data are available')
end

if(year>yearList(end))
    %If the year is after the final year, the prediction coefficients
    %in the final column of rowData must be used to predict forward the end
    %coefficients.
    idx=length(HeaderData)-1;
    putCoeffsIntoCS(rowData,C,S,idx);
    
    %The slopes for interpolation.
    putCoeffsIntoCS(rowData,C1,S1,idx+1);
    yearDiff=year-yearList(end);
    
    %Perform linear interpolation.
    for n=0:M
        for m=0:n
            C(n+1,m+1)=C(n+1,m+1)+yearDiff*C1(n+1,m+1);
            S(n+1,m+1)=S(n+1,m+1)+yearDiff*S1(n+1,m+1);
        end
    end
else
    [val,idx]=sort(abs(yearList-year));
    
    if(val(1)==0)%If it perfectly matched a year, then no extrapolation is needed.
        putCoeffsIntoCS(rowData,C,S,idx(1)+3);
    else%Otherwise, get the final two years and perform linear extrapolation between them.
        putCoeffsIntoCS(rowData,C,S,idx(1)+3);
        putCoeffsIntoCS(rowData,C1,S1,idx(2)+3);

        yearSpan=yearList(idx(2))-yearList(idx(1));
        yearDiff=year-yearList(idx(1));
        
        %Perform linear interpolation between the points and put the result
        %into S and C.
        for n=0:M
           for m=0:n
               CDiff=C1(n+1,m+1)-C(n+1,m+1);
               SDiff=S1(n+1,m+1)-S(n+1,m+1);

               C(n+1,m+1)=C(n+1,m+1)+yearDiff*CDiff/yearSpan;
               S(n+1,m+1)=S(n+1,m+1)+yearDiff*SDiff/yearSpan;
           end
        end
    end
end

%Change the units fron Nanotesla to Tesla.
C(:)=10^(-9)*C(:);
S(:)=10^(-9)*S(:);

%If the coefficients should be fully normalized.
if(fullyNormalize~=false)
     for n=0:M
        k=1/sqrt(1+2*n);
        for m=0:n
            C(n+1,m+1)=k*C(n+1,m+1);
            S(n+1,m+1)=k*S(n+1,m+1);
        end
     end
end

%Return S and C as arrays, not as a CountingClusterSet classes.
C=C.clusterEls;
S=S.clusterEls;

a=Constants.WMM2010SphereRad;%meters
c=a^2;
end

function putCoeffsIntoCS(rowData,C,S,idx)
    %Put the data from column idx for all rows of rowData into the S and C
    %matrices.
    numRows=length(rowData);
    
    for curRow=1:numRows
        %Extract the order and degree.
        temp=textscan(rowData{curRow}{2},'%f');
        n=temp{1};
        temp=textscan(rowData{curRow}{3},'%f');
        m=temp{1};
        
        temp=textscan(rowData{curRow}{1},'%c');
        if(temp{1}=='g')%If the row consists of coefficients for C
            temp=textscan(rowData{curRow}{idx},'%f');
            C(n+1,m+1)=temp{1};
        elseif(temp{1}=='h')%If the row consists of coefficients for S
            temp=textscan(rowData{curRow}{idx},'%f');
            S(n+1,m+1)=temp{1};
        else
            error('Invalid data value read.') 
        end
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
