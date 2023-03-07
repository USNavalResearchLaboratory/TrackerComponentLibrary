function [geoidHeight,coeffData]=getEGMGeoidHeight(latLon,tideSys,useNGAApprox,modelType,coeffData)
%%GETEGMGEOIDHEIGHT  Given a point as latitude and longitude in WGS-84
%                  ellipsoidal coordinates, obtain the geoid height (geoid
%                  undulation) using the National Geospatial Intelligence
%                  Agency's (NGA's) Earth Gravitation Model 2008 (EGM2008)
%                  or Earth Gravitation Model 1996 (EGM96) with
%                  the associated set of coefficients for transforming the
%                  height anomaly on the reference ellipsoid into a geoid
%                  undulation.
%
%INPUTS: latLon A 2XN matrix of N points of the form [latitude;longitude]
%               given with respect to the WGS-84 reference ellipsoid where
%               the geoid height is to be evaluated. If a batch of points
%               on a grid is evaluated, the algorithm is significantly
%               faster if points having the same latitudes are grouped
%               together.
%       tideSys A number indicating the tide system in use. The values are
%               0: Conventional tide free. (The default if this parameter
%                  is omitted.)
%               1: Mean tide
%               2: Zero tide
% useNGAApprox As opposed to computing the zeroth order term of the
%              disturbing potential T and the term -(W0-U0)/gamma using the
%              gamma (magnitude of the acceleration due to gravity) that is
%              appropriate for the point in question, just use the constant
%              factor of -41cm (EGM2008) or -53cm (EGM96) that the NGA
%              uses. Also, omit a correction for differences in the value
%              of GM and the Earth's semi-major axis when computing the
%              disturbing potential T. Using this approximation makes the
%              results match the output of the National Geospatial
%              Intelligence Agency's (NGA's) code to three places after the
%              decimal point --the accuracy at which their Fortran prints
%              the output. If this parameter is omitted, the default is
%              false: the more precise method of computing the zeroth-order
%              term is used. However, one cannot assume that the estimated
%              results are more accurate.
%    modelType An optional parameter specifying coefficient model to load.
%              Possible values are
%              0 (The default if omitted) Load the EGM2008 model.
%              1 Load the EGM96 model.
%    coeffData A set of pre-loaded coefficients that can speed up the
%              computation by eliminating the need to compute them on the
%              fly. coeffData is a structure with elements C, S, CZeta and
%              SZeta where C and S have been set by the output of
%              getEGMWGS84TCoeffs and CZeta and SZeta have been set by the
%              output of getEGMZeta2NCoeffs. If this parameter is omitted,
%              the coefficients are computed during the execution of the
%              function.
%
%OUTPUTS: geoidHeight The geoid height in meters.
%           coeffData The coeffData coefficients that can be passed to
%                     another call of getEGMGeoidHeight to make it faster.
%
%Note that this function will be very slow if one hasn't called
%CompileCLibraries to compile the spherical harmonic synthesis functions.
%
%The geoid is a theoretical surface of constant gravitational potential.
%The potential used for the geoid is that implied on 
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%where the potential of the geoid is defined as the potential on the
%surface of a reference ellipsoid having semi-major axis a=6378136.58 and
%flattening factor f= 1/298.257686, presumably with the same value of the
%Earth rotation rate (omega) and the universal gravitational constant times
%the mass of the Earth (GM) as is present in the WGS-84 standard. Thus,
%using the function ellipsParam2Grav, the potential on the geoid can be
%found.
%
%The formula used is based on Equation 4 in the paper by Rapp [1], where
%the values of C_1 and C_2 in the paper are provided in terms of spherical
%harmonic coefficients by the NGA in the file Zeta-to-N_to2160_egm2008,
%which is loaded using the function getEGMZeta2NCoeffs(). The NGA did not
%document how the coefficients were determined. However, they require the,
%use of digital elevation data as well as data specifying the location of
%the coastlines. Additionally, high-precision computations require
%information on the density of the Earth's crust.
%
%However, Equation 4 in Rapp's paper can not be directly used because the
%reference ellipsoid does not use the same values of GM and a for the
%WGS-84 ellipsoid, with respect to which heights are evaluated, as for the
%EGM2008 gravitational potential coefficients and the geoid potential is
%not the same as the potential implied by the WGS-84 reference ellipsoid.
%Thus, the modified expression for computing the disturbing potential T of
%Equation 5 in [2] is used and the extra term in Brun's equation for the
%difference between the geoid and ellipsoid potentials of Equation 2 of the
%Smith paper is also used. The latter equations is also given in Chapter
%2.17 (Equation 2-347) of [3].
%
%This function computes the tide-free geoid height and then uses the
%conversion equations of [4] to convert to other systems. Note that a loss
%of precision compared to directly computing the geoid height in the other
%systems might be present. The differences between the tide models are
%complicated and are outlined in Figure 1.2 on page 17 in Section 1.2 of 
%[5].
%
%An overview of the math underlying the NGA's method for computing the
%geoid undulations is given in [6].
%
%With the EGM96 model, one can verify the routine using the test values
%from the NGA. These are:
%Latitude and longitude of points converted to radians.
% latLon= [38.628155  269.779155;  
%         -14.621217  305.021114; 
%          46.874319  102.448729; 
%         -23.617446  133.874712; 
%          38.625473  359.999500; 
%         -00.466744    0.002300]'*(pi/180);
% 
% %Use the EGM96 model.
% modelType=1;
% %Get the tide-free undulation.
% tideSys=0;
% %Use the same approximations that the NGA publishes.
% useNGAApprox=true;
% geoidHeight=getEGMGeoidHeight(latLon,tideSys,useNGAApprox,modelType);
% 
% %The output should (to all significant digits listed), match the test
% %output data given by the NGA for those points. That is:
% geoidHeight=[-31.629;
%               -2.966;
%              -43.572;
%               15.868;
%               50.065;
%               17.330];
%
%The sample values from the NGA for the EGM96 model were taken from
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
%
%Similarly, the results can be compared to the NGA's implementation using
%the test values it proves for the EGM2008 model. These are:
% latLon=[37.000000  241.000000;   
%         37.000000 -119.000000;   
%         36.000000  242.983333;   
%         90.000000    0.000000;    
%        -90.000000  359.983333;  
%        -90.000000    0.000000]'*(pi/180);
% 
% %Use the EGM2008 model.
% modelType=0;
% %Get the tide-free undulation.
% tideSys=0;
% %Use the same approximations that the NGA publishes.
% useNGAApprox=true;
% geoidHeight=getEGMGeoidHeight(latLon,tideSys,useNGAApprox,modelType);
%
% %The output should (to all significant digits listed), match the test
% %output data given by the NGA for those points. That is:
% geoidHeight=[-26.151;
%              -26.151;
%              -29.170;
%               14.899;
%              -30.150;
%              -30.150];
%
%REFERENCES:
%[1] R. H. Rapp, "Use of potential coefficient models for geoid undulation
%    determinations using a spherical harmonic representation of the height
%    anomaly/geoid undulation difference," Journal of Geodesy, vol. 71,
%    no. 5, pp. 282-289, Apr. 1997.
%[2] Smith, D. A., "There is no such thing as "The" EGM96 geoid: Subtle
%    points on the use of a global geopotential model" IGeS Bulletin No. 8,
%    International Geoid Service, Milan, Italy, p. 17-28, 1998.
%    http://www.ngs.noaa.gov/PUBS_LIB/EGM96_GEOID_PAPER/egm96_geoid_paper.html
%[3] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%[4] M. Ekman, "Impacts of geodynamic phenomena on systems for height
%    and gravity," Bulletin Géodésique, vol. 63, no. 3, pp. 281-296, 1996.
%[5] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%[6] F. G. Lemoine, S. Kenyon, J. Factor, R. G. Trimmer, N. K. Pavlis, and
%    et. al., "The development of the joint NASA GSFC and the National
%    Imagery and Mapping Agency (NIMA) geopotential model EGM96,"
%    National Aeronautics and Space Administration, Goddard Space Flight
%    Center, Greenbelt, MD, Tech. Rep. NASA/TP-1998-206861, Jul. 1998.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    tideSys=0;
end

if(nargin<3)
    useNGAApprox=false;
end

numPoints=size(latLon,2);

if(nargin<4)
    modelType=0;%EGM2008 Model
end

%The constants for the reference ellipsoid
GMEllipse=Constants.WGS84GMWithAtmosphere;
aEllipse=Constants.WGS84SemiMajorAxis;
omegaEllipse=Constants.WGS84EarthRotationRate;
fEllipse=Constants.WGS84Flattening;

%The constants for the gepotential model. The values are the same in the
%EGM96 model as in the EGM2008 model.
GMPotential=Constants.EGM2008GM;
aPotential=Constants.EGM2008SemiMajorAxis;

if(nargin<5)
    %Get the coefficients for the disturbing potential
    [C,S]=getEGMWGS84TCoeffs([],useNGAApprox,modelType);    
    coeffData.C=C;
    coeffData.S=S;
    
    %The degree and order of the correction term model.
    [C,S]=getEGMZeta2NCoeffs([],modelType);
    coeffData.CZeta=C;
    coeffData.SZeta=S;
end
    
%The point on the reference ellipsoid in Cartesian coordinates.
cartPoint=ellips2Cart([latLon;zeros(1,numPoints)],aEllipse,fEllipse);
r=sqrt(sum(cartPoint.*cartPoint,1));
spherPoint=[r;ellips2Sphere(latLon)];
%Due to precision limitations, it is possible that points having the same
%latitude will have spherical radii that differ by a value close to eps.
%This loop makes sure that points of the same latitude have the same radius
%value. That can greatly accelerate the evaluation of the spherical
%harmonic coefficients when the points are evaluated on a grid and have
%been sorted by latitude.
for curPoint=2:numPoints
    if(spherPoint(3,curPoint)==spherPoint(3,curPoint-1))
        spherPoint(1,curPoint)=spherPoint(1,curPoint-1);
    end
end

%Compute the disturbing potential at the lat-lon-point.
T=spherHarmonicEval(coeffData.C,coeffData.S,spherPoint,aPotential,GMPotential);
[U0,g]=ellipsParam2Grav(cartPoint,omegaEllipse,aEllipse,fEllipse,GMEllipse);

%Acceleration due to gravity on the reference ellipsoid under an
%ellipsoidal Earth.
gamma=sqrt(sum(g.*g,1))';%The norm of the columns of the matrix

%The height anomaly at the surface of the reference ellipsoid. We can
%compute it using the fixed -41cm (EGM2008) or -53cm (EGM96) term of the
%NGA, or we can compute it based upon the actual geopotential.
if(useNGAApprox~=false)
    if(modelType==0)
        correctionTerm=-0.41;
    else
        correctionTerm=-0.53;
    end
    
    %The second term gets rid of the contribution of the C_(0,0) term in
    %T/gamma.
    zeta=T./gamma-coeffData.C(1)*GMPotential./(gamma.*r')+correctionTerm;
else
    %This W0 is the defining gravitational potential of the geoid. It is not
    %the same as U0, the gravitational potential of the reference ellipsoid.
    W0=ellipsParam2Grav([],omegaEllipse,Constants.EGM2008GeoidSemiMajorAxis,Constants.EGM2008GeoidFlattening,GMEllipse);

    zeta=T./gamma-(W0-U0)./gamma;
end

zeta2N=spherHarmonicEval(coeffData.CZeta,coeffData.SZeta,spherPoint(2:3,:));
%The division by 100 for the EGM96 model is not documented and was
%ascertained by noticing that the corrections all appeared off by a factor
%of 100. In the NGA's Fortran function F477.F, which implements a geoid
%undulation computation, the factor of 100 division was found to be
%present.
if(modelType==0)
    geoidHeight=zeta+zeta2N;
else
    geoidHeight=zeta+zeta2N/100;
end

%Convert the appropriate tide system.
switch(tideSys)
    case 1%mean tide
        k=0.3;%A Love number.
        %The 10^(-2) converts from centimeters to meters. 
        geoidHeight=geoidHeight+(1+k)*(9.9-26.9*sin(latLon(1,:)').^2)*10^(-2);
    case 2%zero tide
        k=0.3;%A Love number.
        geoidHeight=geoidHeight+k*(9.9-29.6*sin(latLon(1,:)').^2)*10^(-2);
    otherwise%The tide-free system does not need a conversion.
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
