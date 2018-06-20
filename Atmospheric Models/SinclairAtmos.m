function [n,dndr,T,P]=SinclairAtmos(h,phl0,Rh,P0,T0,wl,ht)
%SINCLAIRATMOS   Compute atmospheric parameters for the Sinclair
%                atmospheric model. This simple model is often used for
%                computing atmospheric refraction when viewing stars. The
%                index of refraction in the model depends on the height
%                above mean sea level. In this implementation, mean sea
%                level is approximated by the surface of the WGS-84
%                reference ellipsoid. The model is not sufficiently precise
%                that geoid undulations should matter and it also shouldn't
%                matter if a geocentric latitude is substituted for the
%                geodetic latitude.
%
%INPUTS: h  The vector or matrix of heights above sea level where the
%           indices of refraction and other parameters are desired.
%      phl0 The location of an observer in WGS-84 ellipsoidal coordinates
%           [latitude;longitude;ellipsoidal height] in the troposphere
%           measuring the relative humidity, ambient air pressure and the
%           temperature. Only the latitude in radians and the ellipsoidal
%           height in meters matter; the longitude is ignored.
%       Rh  The relative humidity at the observer (between 0 and 1). If
%           this parameter is omitted or an empty matrix is passed, then
%           Constants.standardRelHumid is used.
%       P0  The atmospheric pressure at the observer in Pascals (N/m^2). If
%           this parameter is omitted or an empty matrix is passed, then
%           Constants.standardAtmosphericPressure is used.
%       T0  The air temperature at the observer in degrees Kelvin. If this
%           parameter is omitted or an empty matrix is passed, then
%           Constants.standardTemp is used.
%       wl  The wavelength at which the observation is made in units of
%           meters. If this parameter is omitted or an empty matrix is
%           passed, then a wavelength of 0.574 micrometers is used, which
%           is in the visible spectrum (a rather yellow color). This
%           parameter only matters for determining the narrowband index of
%           refraction at a given .
%       ht  The assumed height of the troposphere in meters. This parameter
%           specifies when a stratospheric model begins being used. If this
%           parameter is omitted, then the default value of 11000m is used.
%
%OUTPUTS: n A matrix providing the index of refraction at all of the
%           distances given in r.
%      dndr A matrix providing the derivative of the index of refraction
%           with respect to r evaluated at all of the values given in r.
%        T  A matrix giving the temperature in degrees Kelvin in the
%           atmospheric model at the distances given in r.
%        P  A matrix giving the barometric pressures in Pascals  in the
%           atmospheric model at the distances given in r.
%
%The Sinclair atmospheric model is described in Chapter 7.2 of [1] and in
%[2]. The original source is cited in [2] as being [3]. Reference 3 was not
%consulted in implementing this algorithm.
%
%In the references, the atmosphere is defined in terms of radial
%distances from the center of the Earth. However, the radials distances
%only appear as differences, so heights above mean sea level (MSL) can be
%substituted. MSL in this implementations is just taken to be the surface
%of the WGS-84 reference ellipsoid. The model makes a number of local
%spherical Earth approximations.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%[2] C. Y. Hohenkerk and A. T. Sinclair, "The computation of angular
%    atmospheric refraction at large zenith angles," United Kingdom
%    Hydrographic Office, HM Nautical Almanac, Tech. Rep. 63, Apr. 1985.
%    http://astro.ukho.gov.uk/data/tn/naotn63.pdf
%[3] A. T. Sinclair, "The effect of atmospheric refraction on laser ranging
%    data," United Kingdom Hydrographic Office, HM Nautical Almanac, Tech.
%    Rep. 59, 1982.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<7||isempty(ht))
       ht=11000;%Assumed top of the troposphere in meters.
    end

    if(nargin<6||isempty(wl))
       wl=0.574e-6; 
    end

    if(nargin<5||isempty(T0))
       T0=Constants.standardTemp; 
    end
    
    if(nargin<4||isempty(P0))
        P0=Constants.standardAtmosphericPressure;
    end
    
    if(nargin<3||isempty(Rh))
        Rh=Constants.standardRelHumid;
    end

    %The geodetic latitude in radians.
    phi=phl0(1);
    %The height above the reference ellipsoid in meters is used in
    %place of the height above the geoid.
    h0=phl0(3);
    phi=ellips2Sphere(phi);
    
    %%SET ATMOSPHERIC MODEL PARAMETERS.
    %Convert the pressure from Pascals to millibars.
    P0=P0*0.01;
    %Convert the wavelength from meters to micrometers.
    lambda=wl*1e6;
    %The set of constants is from Equation 7.82 in [1].
    R=1000*Constants.molarGasConstant;%(R) J/(kilomol K)
    Md=28.966;%Assumed molecular mass (im amu) of dry air.
    %The molecular mass (in amu) of water.
    Mw=2*Constants.elementAMU(1)+Constants.elementAMU(8);
    %Exponent of the temperature dependence of water vapor pressure.
    %This is used in an approximate conversion to get the partial
    %pressure of water vapor from the relative humidity.
    delta=18.36;
    alpha=0.0065;%Tropospheric lapse rate of temperature in K/m.
    
    %The following set of Equations is from Equation 7.83 in [1].
    %The partial pressure of water vapor at the observer in millibars.
    Pw0=Rh*(T0/247.1)^delta;
    gBar=9.784*(1-0.0026*cos(2*phi)-0.00000028*h0);
    A=(287.604+1.6288/lambda^2+0.0136/lambda^4)*(273.15/1013.25)*1e-6;
    C2=gBar*Md/R;
    gamma=C2/alpha;
    C5=Pw0*(1-Mw/Md)*gamma/(delta-gamma);
    C6=A*(P0+C5)/T0;
    C7=(A*C5+11.2684e-6*Pw0)/T0;
    C8=alpha*(gamma-1)*C6/T0;
    C9=alpha*(delta-1)*C7/T0;
    
    %Allocate space for the return variables.
    n=zeros(size(h));
    dndr=zeros(size(h));
    T=zeros(size(h));
    P=zeros(size(h));

    %First, deal with all of the points that are in the stratosphere.
    sel=h>ht;
    
    %First, compute the refraction and temperature at the top of the
    %troposphere using Equation 7.85 in [1].
    Tt=T0-alpha*(ht-h0);
    T(sel)=Tt;
    TRat=Tt/T0;
    nt=1+(C6*TRat^(gamma-2)-C7*TRat^(delta-2))*TRat;
    %Approximate partial vapor pressure of water at the tropopause,
    %taken from page 4 of [2].
    Pwt=Pw0*TRat^delta;
    %Air pressure at the tropopause taken from page 4 of [2].
    Pt=(P0+C5)*TRat^(gamma)-Pwt*(1-Mw/Md)*gamma/(delta-gamma);

    %The modelled refraction in the stratosphere.
    %Equation 7.86 in [1].
    n(sel)=1+(nt-1).*exp(-C2*(h(sel)-ht)./Tt);
    dndr(sel)=-(C2./Tt).*(nt-1).*exp(-C2*(h(sel)-ht)./Tt);

    %Stratospheric air pressure, taken from page 5 of [2]. The 0.01
    %term converts the pressure from millibars to Pascals.
    P(sel)=Pt.*exp(-C2*(h(sel)-ht)./Tt)/0.01;
  
    %Next, deal with the points in the troposphere.
    sel=~sel;
    %The modelled refraction in the troposphere.
    %Equation 7.85 in [1].
    T(sel)=T0-alpha*(h(sel)-h0);
    TRat=T(sel)/T0;
    n(sel)=1+(C6*TRat.^(gamma-2)-C7*TRat.^(delta-2)).*TRat;
    dndr(sel)=-C8*TRat.^(gamma-2)+C9*TRat.^(delta-2);

    %Approximate partial vapor pressure of water at altitude, taken
    %from page 4 of [2].
    Pw=Pw0*TRat.^delta;
    %Air pressure at altitude taken from page 4 of [2]. The 0.01 term
    %converts the pressure from millibars to Pascals.
    P(sel)=((P0+C5).*TRat.^(gamma)-Pw*(1-Mw/Md)*gamma/(delta-gamma))/0.01;
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
