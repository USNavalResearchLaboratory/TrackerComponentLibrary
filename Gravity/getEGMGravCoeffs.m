function [C,S,a,c,CStdDev,SStdDev]=getEGMGravCoeffs(M,isTideFree,modelType,TT1,TT2,adjType,useIERSC20)
%%GETEGMGRAVCOEFFS     Get fully normalized spherical harmonic coefficients
%                      for the National Geospatial Intelligence Agency's
%                      (NGA's) Earth Gravitation Model 2008 (EGM2008) or
%                      Earth Gravitation Model 1996 (EGM96) models up to a
%                      desired order. The coefficients can be for the
%                      zero-tide or the tide-free models. Note that
%                      time-varying components except for the drift of the
%                      coefficients over time and the pole tide (if
%                      requested) are not included and must be added
%                      separately. The coefficients can be used with the
%                      spherHarmonicEval function to provide an
%                      acceleration with respect to the WGS-84 coordinate
%                      system (the ITRS) but WITHOUT the Coriolis term or
%                      any other terms that would arise due to the fact
%                      that the coordinate system is non-inertial.
%
%INPUTS:    M       The integer maximum order of the spherical harmonic
%                   coefficients obtained. This is a value between 2 and
%                   2190 for the EGM2008 model and between 2 and 360 for
%                   the EGM96 model. If this parameter is omitted or an
%                   empty matrix is passed, the default value being the
%                   total number of coefficients in the model is used.
%       isTideFree  A boolean value indicating whether the coefficients
%                   should be for a theoretical tide-free Earth or for a
%                   zero-tide model. If this parameter is omitted, the
%                   default value is isTideFree=true;
%        modelType  An optional parameter specifying the coefficient model
%                   to load. Possible values are
%                   0 (the default if omitted) Load the EGM2008 model.
%                   1 Load the EGM96 model.
%         TT1,TT2   An optional input. Two parts of a Julian date given in
%                   terrestrial time (TT). The units of the date are days.
%                   The full date is the sum of both terms. The date is
%                   broken into two parts to provide more bits of
%                   precision. It does not matter how the date is split.
%                   The use of the date depends on the subsequent parameter
%                   adjType.
%           adjType An optional parameter specifying what adjustments
%                   (if any) should be made to the coefficients. Possible
%                   values are
%                   0) (The default if omitted) Do not make any
%                      adjustments. The TT values (if provided) will not be
%                      used.
%                   1) Adjust for the drift of the low-order coefficients
%                      over time as given by the respective models. The
%                      drifts over time occur for many reasons, such as
%                      because of the rebound of the tectonic plates after
%                      the glaciers receded at the end of the last ice age.
%                      The TT values are required. The
%                      getGravCoeffOffset4Drift function is used. The
%                      results are more accurate (but less 
%                   2) Adjust for polar motion. Polar motion
%                      is a change in the rotation axis of the Earth with
%                      respect to the crust and according to the 2010 IERS
%                      conventions, it affects the lower-order spherical
%                      harmonic terms. The getAdjustedGravCoeffs4PolarMotion
%                      function is used after the drift of the low-order
%                      coefficients is taken into account, but before the
%                      tide-free model is converted to a zero-tide model
%                      (if desired).
%                   3) Adjust for the drift of the coefficients and for
%                      polar motion.
%      useIERSC20  This input is only used if the EGM2008 model is chosen.
%                  The IERS 2010 conventions provides the drift rates for
%                  the gravitational terms, but defines a different more
%                  accurate value of C(2+1,0+1). If this is true, then the
%                  more accurate value is used. If this is false, then the
%                  value in the EGM2008 model is used. If omitted, the
%                  default is false if adjType=0 and true for any other
%                  value of adjType.
%
%OUTPUTS:   C   A ClusterSet class holding the coefficient terms that are
%               multiplied by cosines in the harmonic expansion. C(n+1,m+1)
%               is the coefficient of degree n and order m. When a maximum
%               degree of M is used, all C have values for all n from 0 to
%               M and for all m from 0 to n for each n. The coefficients
%               are unitless.
%           S   A ClusterSet class holding the coefficient terms that are
%               multiplied by sines in the harmonic expansion. The
%               format of S is the same as that of C.
%           a   The numerator in the (a/r)^n term in the spherical harmonic
%               sum, having units of meters.
%           c   The constant value by which the spherical harmonic series
%               is multiplied, having units of m^3/s^2.
%       CStdDev The standard deviation of each C coefficient.
%       SStdDev The standard deviation of each S coefficient.
%
%The EGM2008 coefficients are available from the NGA on the site:
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%where one should download the "EGM2008 Tide Free Spherical Harmonic
%Coefficients", which is the file EGM2008_to2190_TideFree.gz. The text file
%should be ungzipped and placed in the data folder that is in the same
%folder as this script. If no external program is available to unzip the
%file, the Matlab function gunzip can be used.
%
%The EGM96 coefficients are available from the NGA on the site:
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
%The coefficients should be decompressed (resulting in a text file named
%egm96 without a file extension). The decompressed file should be placed in
%the data folder that is in the same folder as this script.
%
%This function first checks for a .mat file with the coefficients in it.
%The .mat file is small and can be read quite quickly. However, if one does
%not exist, then it tries to read the EGM2008_to2190_TideFree text file (or
%the egm96 text file) that one can obtain directly from the NGA. If
%M is the maximum number of coefficients or is empty and the text file is
%read directly, then the .mat file is created so that subsequent reads are
%faster. Note that after the .mat file has ben created, the text file can
%be deleted.
%
%More on using the spherical harmonic coefficients is given in
%the comments for the function spherHarmonicEval and the format and use of
%the coefficients is also documented in
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/README_FIRST.pdf
%
%The EGM2008 model is documented in [1] and [2]. The EGM96 model is
%documented in [3].
%
%The IERS conventions are in [4].
%
%REFERENCES:
%[1] N. K. Pavlis, S. A. Holmes, S. C. Kenyon, and F. J. K., "The
%    development and evaluation of the Earth gravitational model 2008
%    (EGM2008)," Journal of Geophysical Research, vol. 117, no. B4, Apr.
%    2012.
%[2] N. K. Pavlis, S. A. Holmes, S. C. Kenyon, and F. J. K., "Correction to
%    "The development and evaluation of the Earth gravitational model 2008
%    (EGM2008)"," Journal of Geophysical Research, vol. 118, 2013.
%[3] F. G. Lemoine, S. Kenyon, J. Factor, R. G. Trimmer, N. K. Pavlis, and
%    et. al., "The development of the joint NASA GSFC and the National
%    Imagery and Mapping Agency (NIMA) geopotential model EGM96," National
%    Aeronautics and Space Administration, Goddard Space Flight Center,
%    Greenbelt, MD, Tech. Rep. NASA/TP-1998-206861, Jul. 1998.
%[4] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6)
    adjType=0;
end

if(nargin<7)
    if(adjType==0)
       useIERSC20=false;
    else
       useIERSC20=true;
    end
end

if(nargin<5)
    TT1=[];
    TT2=[];
end

if(nargin<3)
    modelType=0;%EGM2008 Model
end

if(nargin<2)
    isTideFree=true;
end

if(nargin<1||isempty(M))
    if(modelType==0)
        M=2190;%EGM2008 Model
    else
        M=360;%EGM96 Model
    end
end

if(modelType==0&&M>2190)
    error('The EGM2008 model only has coefficients up through degree 2190');
elseif(modelType==1&&M>360)
    error('The EGM96 model only has coefficients up through degree 360');
end

%The same values of a and c are used in the EGM2008 mnodel as in the EGM96
%model.
a=Constants.EGM2008SemiMajorAxis;
c=Constants.EGM2008GM;

if(modelType==0)
    fileName='/data/EGM2008_to2190_TideFree';
else
    fileName='/data/egm96';
end

totalNumCoeffs=(M+1)*(M+2)/2;

%The EGM20008 and EGM96 coefficient data files should be located in a data
%folder that is in the same folder as this file. This finds the path to
%this file. 
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

%First, see if a .mat file with all of the data exists. If so, then use
%that and ignore everything else.
if(exist([ScriptFolder,fileName,'.mat'],'file'))
    load([ScriptFolder,fileName,'.mat'],'CCoeffs','SCoeffs','CCoeffsStdDev','SCoeffsStdDev','clustSizes','offsets');
    %Create the ClusterSet classes to hold the data.
    C=ClusterSet();
    S=ClusterSet();
    CStdDev=ClusterSet();
    SStdDev=ClusterSet();
    
    C.clusterEls=CCoeffs(1:totalNumCoeffs);
    S.clusterEls=SCoeffs(1:totalNumCoeffs);
    CStdDev.clusterEls=CCoeffsStdDev(1:totalNumCoeffs);
    SStdDev.clusterEls=SCoeffsStdDev(1:totalNumCoeffs);
    
    C.clusterSizes=clustSizes(1:(M+1));
    S.clusterSizes=C.clusterSizes;
    CStdDev.clusterSizes=C.clusterSizes;
    SStdDev.clusterSizes=C.clusterSizes;
    
    C.offsetArray=offsets(1:(M+1));
    S.offsetArray=C.offsetArray;
    CStdDev.offsetArray=C.offsetArray;
    SStdDev.offsetArray=C.offsetArray;
else
    %If the .mat file does not exist, then assume that the coefficients
    %must be read from the text file provided by the NGA and, if all of the
    %coefficients were read in, then create the .mat file after reading in 
    %the coefficients so that future reads are faster

    %Allocate space for the coefficients.
    emptyData=zeros(totalNumCoeffs,1);
    clustSizes=1:(M+1);
    C=ClusterSet(emptyData,clustSizes);
    S=ClusterSet(emptyData,clustSizes);
    CStdDev=ClusterSet(emptyData,clustSizes);
    SStdDev=ClusterSet(emptyData,clustSizes);

    %Read in the data up to the specified order. The -3 deals with the fact
    %that coefficients for degrees and order (0,0), (1,0) and (1,1) have 
    %been omitted from the data file.
    fileID=fopen([ScriptFolder,fileName]);
    data=textscan(fileID,'%d %d %f %f %f %f',totalNumCoeffs-3);
    fclose(fileID);

    %Put the data into the ClusterSet class instances. The data is already
    %appropriately arranged by degree and order starting with degree 2 and
    %order 0.
    C.clusterEls(4:end)=data{3};
    S.clusterEls(4:end)=data{4};
    CStdDev.clusterEls(4:end)=data{5};
    SStdDev.clusterEls(4:end)=data{6};

    %The EGM96 and EGM2008 coefficient files omits the C_{0,0} value. It is
    %explicitly added here.
    C(0+1,0+1)=1;

    %If all of the coefficients have been read, then save the data as a
    %.mat file so that it can be read faster in the future.
    if(modelType==0&&M==2190||modelType==1&&M==360)
        CCoeffs=C.clusterEls;
        SCoeffs=S.clusterEls;
        CCoeffsStdDev=CStdDev.clusterEls;
        SCoeffsStdDev=SStdDev.clusterEls;
        clustSizes=C.clusterSizes;
        offsets=C.offsetArray;

        save([ScriptFolder,fileName,'.mat'],'CCoeffs','SCoeffs','CCoeffsStdDev','SCoeffsStdDev','clustSizes','offsets');
    end
end

%If the value for the EGM2008 C20 term from the IERS2008 model should be
%used. The C20 terms are kept separately, because it is implied in the IERS
%conventions that
if(modelType==0&&useIERSC20==true)
    %The zero-tide value is from Table 6.2 in the 2010 IERS conventions.
    C20ZeroTide=-0.48416948e-3;
    C20TideFree=-0.48416531e-3;
elseif(modelType==0)
    C20TideFree=C(2+1,0+1);
    %This value is documented on the web site providing the EGM2008
    %coefficients.
    C20ZeroTide=-0.484169317366974e-3;
else%The values for the EGM96 model
    C20TideFree=C(2+1,0+1);
    %The formula is in the readme file accompanying the EGM96 model.
    C20ZeroTide=C(2+1,0+1)-3.11080e-8*0.3/sqrt(5);
end

%If the coefficients should be updated for the drift over time.
if(~isempty(TT1)&&adjType~=0&&adjType~=2)
    %Get the offsets due to drift over time.
    [deltaC,deltaS]=getGravCoeffOffset4Drift(TT1,TT2,modelType);
    
    %Add in the offsets
    numOffset=length(deltaC.clusterEls);
    if(numOffset<=totalNumCoeffs)
        C.clusterEls(1:numOffset)=C.clusterEls(1:numOffset)+deltaC.clusterEls(:);
        S.clusterEls(1:numOffset)=S.clusterEls(1:numOffset)+deltaS.clusterEls(:);
    else
        C.clusterEls(:)=C.clusterEls(:)+deltaC.clusterEls(1:numEl);
        S.clusterEls(:)=S.clusterEls(:)+deltaS.clusterEls(1:numEl);
    end
    
    %Update the separate
    C20TideFree=C20TideFree+deltaC(2+1,0+1);
    C20ZeroTide=C20ZeroTide+deltaC(2+1,0+1);
end

%If the coefficients should be updated for the effects of polar motion.
if(adjType==2||adjType==3)
    %The IERS conventions implies that the zero-tide value of C20 should be
    %used in this computation
    C(2+1,0+1)=C20ZeroTide;
    getAdjustedGravCoeffs4PolarMotion(C,S,TT1,TT2,true);
end

%Set the C20 value for the appropriate model.
if(isTideFree==false)
    C(2+1,0+1)=C20ZeroTide;
else
    C(2+1,0+1)=C20TideFree;
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
