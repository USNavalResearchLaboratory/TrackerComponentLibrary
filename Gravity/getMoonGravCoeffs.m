function [C,S,a,c,CStdDev,SStdDev]=getMoonGravCoeffs(M,isTideFree)
%%GETMOONGRAVCOEFFS Get fully normalized spherical harmonic coefficients
%                   for NASA JPL's GL0900C gravitational model of the
%                   Moon up to a desired order. The coefficients can be for
%                   the zero-tide or the tide-free models. The coefficients
%                   are given with respect to the principle-axis lundar
%                   body fixed coordinate system consistent with NASA's
%                   DE430 ephemerides.
%
%INPUTS: M The integer maximum order of the spherical harmonic coefficients
%          obtained. This is a value between 2 and 900. If this parameter
%          is omitted or an empty matrix is passed, the default value being
%          the total number of coefficients in the model is used.
% isTideFree A boolean value indicating whether the coefficients should be
%          for a theoretical tide-free Earth or for a zero-tide model. If
%          this parameter is omitted, the default value is isTideFree=true;
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the spherical harmonic expansion. This can be given
%           to a CountingClusterSet so that C(n+1,m+1) is the coefficient
%           of degree n and order m. When a maximum degree of M is used,
%           all C have values for all n from 0 to M and for all m from 0 to
%           n for each n. The coefficients are unitless.
%         S A ClusterSet class holding the coefficient terms that are
%           multiplied by sines in the harmonic expansion. The format of S
%           is the same as that of C.
%         a The numerator in the (a/r)^n term in the spherical harmonic
%           sum, having units of meters.
%         c The constant value by which the spherical harmonic series is
%           multiplied, having units of m^3/s^2.
%   CStdDev An array holding the standard deviation of each C coefficient.
%   SStdDev An array holding the standard deviation of each S coefficient.
%
%The lunar gravitational coefficients are available on the site
%http://pds-geosciences.wustl.edu/missions/grail/default.htm
%where the file
%http://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_0900c_sha.tab
%is the actual coefficients and the file
%http://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_0900c_sha.lbl
%is documentation commenting on the coefficients. The text file
%jggrx_0900c_sha.tab should be placed in the data folder that is in the
%same folder as this script.
%
%This function first checks for a .mat file with the coefficients in it.
%The .mat file is small and can be read quite quickly. However, if one does
%not exist, it will try to read the jggrx_0900c_sha.tab file. If M is the
%maximum number of coefficients or is empty and the text file is read
%directly, then the .mat file is created so that subsequent reads are
%faster. Note that after the .mat file has ben created, the text file can
%be deleted.
%
%The GL0900C lunar gravity model is briefly described in [2] and is
%an extention of the 660-degree model of [1]. The coefficients are provided
%in a tide-free system. C20 and C22 coefficients in the zero-tide system
%are provided in [1], but that is for the 660-degree model. However, both
%models use the same second-degree Love number solution for the Moon.
%Assuming a similar solid tide model is used as in Equation 6.6 of Chapter
%6.2.1 of [3] for the Earth, if the masses of the planets remained the
%same, one would expect the offset from tide-free to zero-tide between the
%models to remain the same. However, since the 660-degree model (the
%GL0660B) used the DE421 ephemerides and the newer model uses the DE430
%ephemerides, the mass assigned the moon changed a little. Nonetheless, the
%offset is probably close to the correct value and is what is used here to
%get the zero-tide model.
%
%REFERENCES:
%[1] A. S. Konopliv, R. S. Park, D.-N. Yuan, S. W. Asmar, and et. al.,
%    "The JPL lunar gravity field to spherical harmonic degree 660 from
%    the GRAIL primary mission," Journal of Geophysical Research: Planets,
%    vol. 118, no. 7, pp. 1415-1434, Jul. 2013.
%[2] A. S. Konopliv, R. S. Park, D.-N. Yuan, S. W. Asmar, and et. al.,
%    "High-resolution lunar gravity fields from the GRAIL primary and
%    extended missions," Geophysical Research Letters, vol. 41, no. 5, pp.
%    1452-1458, 16 Mar. 2014.
%[3] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    isTideFree=2;
end

if(nargin<1||isempty(M))
    M=900;%all of the coefficients in the GL0900C model.
end

totalNumCoeffs=(M+1)*(M+2)/2;

fileName='/data/jggrx_0900c_sha.tab';
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

a=Constants.GL0900CSemiMajorAxisMoon;
c=Constants.GL0900CGMMoon;

%First, see if a .mat file with all of the data exists. If so, then use
%that and ignore everything else.
if(exist([ScriptFolder,fileName,'.mat'],'file'))
    load([ScriptFolder,fileName,'.mat'],'C','S','CStdDev','SStdDev');
    %Create the ClusterSet classes to hold the data.
    C=C(1:totalNumCoeffs);
    S=S(1:totalNumCoeffs);
    
    CStdDev=CStdDev(1:totalNumCoeffs);
    SStdDev=SStdDev(1:totalNumCoeffs);
else
    %Read in the data, which, after the first row, is just a table.
    coeffs=dlmread([ScriptFolder,fileName],',',1,0);
    %Get rid of any at the end that are not used. The -1 is because the
    %file does not contain the C_{0,0} term.
    coeffs=coeffs(1:(totalNumCoeffs-1),:);

    emptyData=zeros(totalNumCoeffs,1);
    C=emptyData;
    S=emptyData;
    CStdDev=emptyData;
    SStdDev=emptyData;
    
    C(2:end)=coeffs(:,3);
    S(2:end)=coeffs(:,4);
    CStdDev(2:end)=coeffs(:,5);
    SStdDev(2:end)=coeffs(:,6);
    
    %The model omits the C_{0,0} value. It is explicitly added here.
    C(1)=1;
    
    %If all of the coefficients were read in, then save the result into a
    %.mat file so that it can be read faster next time.
    if(M==900)
        save([ScriptFolder,fileName,'.mat'],'C','S','CStdDev','SStdDev');
    end
end
    
%If a zero-tide model is desired, add in an approximation for the permanent
%tide. 
if(isTideFree==false)
    %These are the zero-tide values from the GL0660B model.
    C20ZeroTide660B=-2.0330530e-4/sqrt(5);
    C22ZeroTide660B=2.242615e-5/sqrt(5/12);
    %These are the tide-Free values from the GL0660B model.
    C20TideFree660B=-0.9088083466222001E-04;
    C22TideFree660B=0.3467379822995000E-04;

    %We will assume the difference between the zero-tide and the tide free
    %models is approxiamtely the same for the GL0660B model and the GL0900C
    %model and will thus add in the same offset.
    %Adjust C_{2,0}
    C(4)=C(4)+(C20ZeroTide660B-C20TideFree660B);
    %Adjust C_{2,2}
    C(6)=C(6)+(C22ZeroTide660B-C22TideFree660B);
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
