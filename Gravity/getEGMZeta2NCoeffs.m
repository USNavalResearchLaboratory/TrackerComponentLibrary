function [C,S]=getEGMZeta2NCoeffs(M,modelType)
%%GETEGMZETA2NCOEFFS Get fully normalized spherical harmonic coefficients
%                    for computing the parameter needed to convert height
%                    anomalies to geoid undulations when using the WGS84
%                    ellipsoid with the Earth Gravitation Model 2008
%                    (EGM2008) or the Earth Gravitation Model 1996 (EGM96)
%                    from the National Geospatial Intelligence Agency (NGA).
%
%INPUTS: M The integer maximum order of the spherical harmonic coefficients
%          obtained. This is a value between 2 and 2160 for the EGM2008
%          model and betwen 2 and 360 for the EGM96 model. If this
%          parameter is omitted or an empty matrix is passed, the default
%          value being the total number of coefficients in the model is
%          used.
% modelType An optional parameter specifying coefficient model to load.
%          Possible values are
%          0 (The default if omitted) Load the height conversion
%             coefficients for the  EGM2008 model.
%          1 Load the coefficients for the EGM96 model.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the harmonic expansion. This can be given to a
%           CountingClusterSet class so that C(n+1,m+1) is the coefficient
%           of degree n and order m. When a maximum degree of M is used,
%           all C have values for all n from 0 to M and for all m from 0 to
%           n for each n. The coefficients have units of meters.
%         S A ClusterSet class holding the coefficient terms that are
%           multiplied by sines in the harmonic expansion. The format of S
%           is the same as that of C.
%
%When using the method of [1] for determining geoid undulations from
%spherical harmonic gravitational models, the coefficients are used for
%going from the zeta term to the N term. The function getEGM2008GeoidHeight
%makes use of the coefficients to compute the geoid height.
%
%The EGM2008 coefficient zeta-to-N set can be obtained at
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%where the "Correction Model" file Zeta-to-N_to2160_egm2008.gz should be
%downloaded, ungzipped and placed in the data folder, which is in the same
%folder as this script. For the EGM96 model, the correction coefficients
%can be downloaded from
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
%where the file of interest is called corrcoef.z which should be
%decompressed and placed in the data folder.
%
%This function first checks for a .mat file with the coefficients in it.
%The .mat file is small and can be read quite quickly. However, if one doe
%not exist, then it tries to read the Zeta-to-N_to2160_egm2008 text file
%(or the corrcoef text file) that one can obtain directly from the NGA. If
%M is the maximum number of coefficients or is empty and the text file is
%read directly, then the .mat file is created so that subsequent reads are
%faster. Note that after the .mat file has ben created, the text file can
%be deleted.
%
%REFERENCES:
%[1] R. H. Rapp, "Use of potential coefficient models for geoid undulation
%    determinations using a spherical harmonic representation of the height
%    anomaly/geoid undulation difference," Journal of Geodesy, vol. 71,
%    no. 5, pp. 282-289, Apr. 1997.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    modelType=0;%EGM2008 Model
end

if(nargin<1||isempty(M))
    if(modelType==0)
        M=2160;%EGM2008 Model
    else
        M=360;%EGM96 Model
    end
end

if(modelType==0&&M>2160)
    error('The EGM2008 model only has coefficients up through degree 2190.');
elseif(modelType==1&&M>360)
    error('The EGM96 model only has coefficients up through degree 360.');
end

totalNumCoeffs=(M+1)*(M+2)/2;

%The name of the coefficient files for the EGM2008 and EGM96 models.
if(modelType==0)
    fileName='/data/Zeta-to-N_to2160_egm2008';
else
    fileName='/data/corrcoef';
end

%The EGM20008 and EGM96 coefficient data files should be located in a data 
%folder that is in the same folder as this file. This finds the path to
%this file. 
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

%First, see if a .mat file with all of the data exists. If so, then use
%that and ignore everything else.
if(exist([ScriptFolder,fileName,'.mat'],'file'))
    load([ScriptFolder,fileName,'.mat'],'C','S');

    C=C(1:totalNumCoeffs);
    S=S(1:totalNumCoeffs);
else
    %If the .mat file does not exist, then assume that the coefficients
    %must be read from the text file provided by the NGA.

    %Read in the data up to the specified order.
    fileID=fopen([ScriptFolder,fileName]);
    data=textscan(fileID,'%d %d %f %f',totalNumCoeffs);
    fclose(fileID);

    %Put the data into the ClusterSet class instances. The data is already
    %appropriately arranged by degree and order starting with degree 2 and
    %order 0.
    C=data{3};
    S=data{4};

    %If all of the coefficients have been read, then save the data as a
    %.mat file so that it can be read faster in the future.
    if(modelType==0&&M==2160||modelType==1&&M==360)

        save([ScriptFolder,fileName,'.mat'],'C','S');
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
