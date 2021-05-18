function dataFile=getHipparcosCat(rectBound,magThresh,onlySingleStars,formatObsParam)
%%GETHIPPARCOSCAT Load the stars in the Hipparcos catalog, removing stars
%                 and parameters according to the input parameters. Only
%                 the five parameter solutions in the Astrometric catalog
%                 are loaded. The specific catalog used in I/311:
%                 Hipparcos, the New Reduction from
%                 http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/311
%                 The units are converted as given below.
%
%INPUTS: rectBound  Specifies a rectangular region (in radians) at the
%                   1991.25 epoch in which stars are returned. rectBound
%                   has the format [RAmin;RAmax;DEmin;DEmax], which
%                   represent minimum and maximum bounds on the right
%                   ascention and declination of the stars to include. An
%                   empty matrix being passed means that all stars should
%                   be included. The default value if this parameter is
%                   omitted is an empty matrix.
%        magThresh  A magnitude threshold. Stars having a Hipparcos
%                   magnitude above this value are omitted. If omitted, a
%                   threshold of 6.5 is used to get rid of stars that are
%                   probably not visible to the naked eye. Setting a
%                   magnitude of positive infinity will guarantee all stars
%                   are included.
%   onlySingleStars If true, stars with more than one component are
%                   not included. The default value is false if this
%                   parameter is omitted.
%   formatObsParam  If true, the output is formatted so that it can be
%                   directly put into the function starCat2Obs to be
%                   converted to the local coordinate system of an observer
%                   on the surface of the Earth. This involves getting rid
%                   of irrelevant fields, converting the units, and
%                   converting everything from the 1995.25 epoch of the
%                   catalog to the J2000.0 epoch. If this parameter is
%                   omitted, then the default is true.
%
%OUTPUT: dataFile   If  formatObsParam=true, then dataFile has the format
% dataFile(:,1) RArad    Right Ascension in (rad) ICRS (epoch J2000.0)
% dataFile(:,2) DErad    Declination in (rad) ICRS (epoch J2000.0)
% dataFile(:,3) Plx      Parallax (rad)
% dataFile(:,4) pmRA     Proper motion in Right Ascension (rad/yr)
%                        in the form dRA/dt and not cos(Dec)*dRA/dt
% dataFile(:,5) pmDE     Proper motion in Declination (rad/yr)
% dataFile(:,6) vRad     Radial velocity of the star in meters per second
%                        with a positive sign meaning that a star is
%                        receding (going away from) Earth. The Hipparcos
%                        catalog does not provide this information, so this
%                        column will be filled with zeros.
%                   On the other hand, if formatObsParam=false, then
%                   dataFile has the format of data in the catalog, which
%                   also has different units on some terms.
% dataFile(:,1) HIP      Hipparcos identifier
% dataFile(:,2) Sn       [0,159] Solution type new reduction (1)
% dataFile(:,3) So       [0,5] Solution type old reduction (2)
% dataFile(:,4) Nc       Number of components. A value greater than
%                        one indicates a system of more than one star.
% dataFile(:,5) RArad    Right Ascension in ICRS, Ep=1991.25 (rad)
% dataFile(:,6) DErad    Declination in ICRS, Ep=1991.25 (rad)
% dataFile(:,7) Plx      Parallax (mas)
% dataFile(:,8) pmRA     Proper motion in Right Ascension (mas/yr)
% dataFile(:,9) pmDE     Proper motion in Declination (mas/yr)
% dataFile(:,10) e_RArad Formal error on RArad (mas)
% dataFile(:,11) e_DErad Formal error on DErad (mas)
% dataFile(:,12) e_Plx   Formal error on Plx (mas)
% dataFile(:,13) e_pmRA  Formal error on pmRA (mas/yr)
% dataFile(:,14) e_pmDE  Formal error on pmDE (mas/yr)
% dataFile(:,15) Ntr     Number of field transits used
% dataFile(:,16) F2      Goodness of fit
% dataFile(:,17) F1      Percentage rejected data
% dataFile(:,18) var     Cosmic dispersion added (stochastic solution)
% dataFile(:,29) ic      Entry in one of the suppl.catalogues
% dataFile(:,20) Hpmag   Hipparcos magnitude. Higher numbers correspond to
%                        fainter stars. Stars above 6.5 probably are not
%                        visible to the naked eye.
% dataFile(:,21) e_Hpmag Error on mean Hpmag
% dataFile(:,22) sHp     Scatter of Hpmag
% dataFile(:,23) VA      [0,2] Reference to variability annex
% dataFile(:,24) B-V     Colour index in the Johnson UBV system.
% dataFile(:,25) e_B-V   Formal error on colour index
% dataFile(:,26)  V-I    V-I colour index in the Cousins' system
% dataFile(:,27:41)  UW  Upper-triangular weight matrix (G1)
%
%The documentation at
%http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/311
%similarly describes the fields.
%Some insight into the meaning of the fields can be found at
%http://www.rssd.esa.int/SA-general/Projects/Hipparcos/pstex/sect2_01.pdf
%However, that document was made for an older version of the catalog.
%
%The Johnson UBV system specifies the color of the star and is described in
%[1].
%
%REFERENCES:
%[1] Johnson, H. L. and Morgan, W. W., "Fundamental stellar photometry for
%    standards of spectral type on the revised system of the Yerkes
%    spectral atlas," Astrophysical Journal, 117, 313, May 1953.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1)
   rectBound=[]; 
end

if(nargin<2)
    magThresh=6.5;
end

if(nargin<3)
    onlySingleStars=false;
end

if(nargin<4)
   formatObsParam=true; 
end

ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

dataFile=dlmread([ScriptFolder,'/data/hip2.dat']);

%The constant to convert milliarcseconds to radians. The multiplications go
%mas->as->arcminutes->deg->rad
mas2Rad=(1/1000)*(1/60)*(1/60)*(pi/180);

%bound everything to the desired region.
if(~isempty(rectBound))
    dataFile=dataFile(~(dataFile(:,5)<rectBound(1)),:); 
    dataFile=dataFile(~(dataFile(:,5)>rectBound(2)),:);
    dataFile=dataFile(~(dataFile(:,6)<rectBound(3)),:); 
    dataFile=dataFile(~(dataFile(:,6)>rectBound(4)),:); 
end

%Get rid of star clusters, if desired.
if(onlySingleStars==true)
    dataFile=dataFile(dataFile(:,4)==1,:);
end

%Get rid of weak-magnitude stars, if desired.
if(isfinite(magThresh))
    dataFile=dataFile(dataFile(:,20)<=magThresh,:);
end

if(formatObsParam==true)
    numStars=size(dataFile,1);
    dataFile=dataFile(:,5:9);
    %The catalog does not store radial velocities, so just put in zeros.
    dataFile=[dataFile,zeros(numStars,1)];
    
    %Convert the parallax from milliarcseconds to arcseconds
    dataFile(:,3)=mas2Rad*dataFile(:,3);
    
    %Convert the proper motion in right ascention term to the format that
    %the starCat2Obs function desires. The Hipparcos catalog provides
    %values times the cosine of the declination. This divides that back out.
    dataFile(:,4)=mas2Rad*dataFile(:,4)./cos(dataFile(:,2));
    %Convert the mas/yr of the proper motion in declination into radians
    %per year.
    dataFile(:,5)=mas2Rad*dataFile(:,5);
    
    %Change the epoch to J2000.0
    dataFile=changeEpoch(dataFile);
    return;
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
