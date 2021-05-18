function [zTrue,deltaZ]=removeAstroRefrac(algorithm,plhObs,z0,Rh,P,T,wl)
%%REMOVEASTROREFRAC Remove the effects of refraction in an observed zenith
%                   angle of an object outside the atmosphere using low-
%                   precision atmospheric models for an observer near the
%                   surface of the Earth viewing an object using a
%                   narrowband sensor.
%
%INPUTS: algorithm This specified the algorithm used. The possible values
%                  are:
%                  0: Use a numerical integration method from Chapter 7.2
%                     of [1] and from [2]. This algorithm uses the inputs
%                     plhObs,z0, Rh, P, T, and wl and makes use of the
%                     Sinclair atmospheric model. This algorithm is the
%                     most precise of all of the methods. The acceptable
%                     zenith distances are capped at 100 degrees, which
%                     should be underground. Values above 100 degrees lead
%                     to zTrue and deltaZ being empty matrices.
%                  1: Use the simple formula of [3] for the refraction
%                     experienced by an observer at sea level viewing light
%                     with a wavelength of 0.574 micrometers (yellow). This
%                     algorithm uses the inputs z0, Rh, P, T. This
%                     algorithm is only valid for zenith distances below 70
%                     degrees. For larger values, zTrue and deltaZ are
%                     returned as empty matrices
%                  2: Use the algorithm from the International Astronomical
%                     Union's (IAU) standards of fundamental astronomy
%                     library. This algorithm uses the inputs z0,Rh,P,T,wl
%                     for an observer at sea level. This algorithm will
%                     produce results for all positive zenith distances,
%                     though the results might not be very good for large
%                     values.
%           plhObs The WGS-84 ellipsoidal latitude, longitude and height
%                  of the observer in meters. Only algorithm 0 uses this
%                  parameter, and only the latitude and height above the
%                  reference ellipsoid are used (the longitude does not
%                  matter). The height is treated as an approximate height
%                  above mean sea level. Given the precision of the
%                  algorithm, the fact that geoid undulations are ignored
%                  probably does not matter.For algorithms 1 and 2, an
%                  empty matrix can be passed.
%               z0 A vector or matrix of observed positive zenith
%                  "distances" (in radians) of the refraction-corrupted
%                  celestial object being observed. A zenith distance is
%                  the observed angle of an object down from the
%                  gravitational vertical (the zenith) at the observer's
%                  location. This value is inaccurate as one approaches
%                  pi/2 radians (the horizon). z0 must be greater than
%                  zero. Given the precision of the algorithms that are
%                  available, an angle with respect to the vertical defined
%                  by the WGS-84 reference ellipsoid could probably be
%                  substituted for an angle with respect to the true
%                  gravitational vertical.
%               Rh The relative humidity at the observer (between 0 and 1).
%                  If this parameter is omitted or an empty matrix is
%                  passed, then Constants.standardRelHumid is used.
%                P The atmospheric pressure at the observer in Pascals
%                  (N/m^2). If this parameter is omitted or an empty matrix
%                  is passed, then Constants.standardAtmosphericPressure is
%                  used.
%                T The air temperature at the observer in degrees Kelvin.
%                  If this parameter is omitted or an empty matrix is
%                  passed, then Constants.standardTemp is used.
%               wl The wavelength at which the observation is made in units
%                  of meters. If this parameter is omitted or an empty
%                  matrix is passed, then a wavelength of 0.574 micrometers
%                  is used, which is in the visible spectrum (a rather
%                  yellow color).
%
%OUTPUTS: zTrue An NX1 vector of the z0 values in radians with atmospheric
%               refraction removed.
%        deltaZ The refraction correction that was applied. zTrue=z0+deltaZ.
%
%The references for Algorithm 0 are [1] and [2]. Unlike the method given in
%those references, the distance from the center of the Earth to the
%reference ellipsoid at the latitude of the observer is used for the radius
%of the spherical Earth in the model as opposed to a constant for the
%entire Earth, in the model.
%
%Algorithm 1 is described in [2] but also in much more detail in [3]. The
%model requires the partial vapor pressure of water. The low-precision
%conversion from a relative humidity to a partical vapor pressure as is
%used in the algorithm in [2] was used here. The loss in accuracy in using
%a low-precision humidity model should be significantly less than the loss
%of precision in the approximations inherent to the Saastamoinen refraction
%model.
%
%Algorithm 2 is directly taken from the IAU's standards of fundamental
%astronomy library. Specifically, whereas the
%documentation describes the value that must be added to get rid of
%refractive biases as
%DeltaZ=A*tan(z0)+B*tan(z0)^3
%where A and B are obtained using the iauRefco function from the IAU's
%library, the implentation used in the function iauAtioq uses a slightly
%modified formula so that reasonable values are obtained for zenith angles
%at or above 90 degrees. Thus, the formula for the correction is copied
%from software version 10 of the function iauAtioq in the IAU's library.
%The documentation for the library lists the accuracy compared to a
%more sophisticated ray tracing method as being  60 milliarcseconds (mas)
%at the worst case in the visible region to 300 mas in the radio region
%(wavelengths greater than 1mm). In both cases, the maximum zenith distance
%(zetaObs) considered was 75 degrees.
%
%Algorithm 2 will provide results without errors and without producing NaNs
%for all input values. The other techniques will produce errors.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%[2] C. Y. Hohenkerk and A. T. Sinclair, "The computation of angular
%    atmospheric refraction at large zenith angles," United Kingdom
%    Hydrographic Office, HM Nautical Almanac, Tech. Rep. 63, Apr. 1985.
%    http://astro.ukho.gov.uk/data/tn/naotn63.pdf
%[3] J. Saastamoinen, "Introduction to the practical computation of
%    astronomical refraction," Bulletin G�od�sique, vol. 106, no. 1,
%    pp. 383-397, Dec. 1972.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<7||isempty(wl))
       wl=0.574e-6; 
    end

    if(nargin<6||isempty(T))
       T=Constants.standardTemp; 
    end
    
    if(nargin<5||isempty(P))
        P=Constants.standardAtmosphericPressure;
    end
    
    if(nargin<4||isempty(Rh))
        Rh=Constants.standardRelHumid;
    end
    
    if(sum(sum(z0<0))>0)
        error('The observed zenith distance must be positive.');
    end
    
switch(algorithm)
    case 0%Use the basic numerical ray tracing algorithm.
        %The geodetic latitude in radians.
        phi=plhObs(1);
        %The height above the reference ellipsoid in meters is used in
        %place of the height above the geoid.
        lambda=plhObs(2);
        h0=plhObs(3);
        %Use the radius of the Earth at this latitude and longitude rather
        %than the standard one provided in the book.
        re=norm(ellips2Cart([phi;lambda;0]));%Units of meters.
        r0=re+h0;%The initial distance from the center of the Earth.
        ht=11000;%Assumed top of the troposphere in meters.
        rt=ht+re;%Radial distance to the top of the troposphere.
        hs=80000;%Height at which refraction is negligible.
        rs=hs+re;%Radial distance to top of the stratosphere.

        %If all zenith distances were allowed, then the algorithm might not
        %terminate.
        if(sum(sum(z0>100*pi/180))>0)
            zTrue=[];
            deltaZ=[];
            return;
            %error('The observed zenith distance should not be too far under the horizon (>100 degrees)');
        end

        if(h0>ht)
           error('The algorithm is not meant for observers above the height of the troposphere (11000m)'); 
        end

        %The index of refraction at the observer, who is assumed to be
        %located in the troposphere.
        n0=SinclairAtmos(h0,plhObs,Rh,P,T,wl,ht);
        %The index of refraction at the tropopause (the top of the
        %troposphere).
        nt=SinclairAtmos(ht,plhObs,Rh,P,T,wl,ht);
        %The index of refraction at the height where it becomes
        %negligible.
        ns=SinclairAtmos(hs,plhObs,Rh,P,T,wl,ht);

        numZ=length(z0(:));
        zTrue=zeros(size(z0));%Allocate space for the return variable.
        deltaZ=zeros(size(z0));%Allocate space for the return variable.

        for curZ=1:numZ
            if(z0(curZ)<1e-20)
            %This deals with numerical precision limitations when using the
            %integral function with small angles. The model has precisely
            %zero refraction when z0=0. If this condition were not here,
            %then the integral function would have problems.
                zTrue(curZ)=z0;
                deltaZ(curZ)=0;
            else
                %The zenith distance of the tropopause. Equation 7.88 in [1].
                zt=asin(n0*r0*sin(z0(curZ))/(nt*rt));
                %The zenith distance of the top of the stratopshere.
                %Equation 7.88 in [1].
                zs=asin(n0*r0*sin(z0(curZ))/(ns*rs));

                %%%IN THE TROPOSPHERE
                refracFunc=@(h)SinclairAtmos(h,plhObs,Rh,P,T,wl,ht);
                intFunc=@(z)integrand(z,z0(curZ),refracFunc);
                xit=-integral(intFunc,zt,z0(curZ));

                %%%IN THE STRATOSPHERE
                refracFunc=@(h)SinclairAtmos(h,plhObs,Rh,P,T,wl,ht);
                intFunc=@(z)integrand(z,z0(curZ),refracFunc);
                xis=-integral(intFunc,zs,zt);

                deltaZ(curZ)=xit+xis;
                zTrue(curZ)=z0(curZ)+deltaZ(curZ);
            end
        end
        return;
    case 1%Use the Saastamoinen formula.
        %Convert the pressure from Pascals to millibars.
        P=P*0.01;
        %Exponent of the temperature dependence of water vapor pressure.
        %This is used in an approximate conversion to get the partial
        %pressure of water vapor from the relative humidity.
        delta=18.36;
        Pw0=Rh*(T/247.1)^delta;
        Q=(P-0.156*Pw0)/T;
        
        if(sum(sum(z0>70*pi/180))>0)
            zTrue=[];
            deltaZ=[];
            return;
            %error('The algorithm is not meant for zenith distances above 70 degrees.');
        end
        
        arcSec2Rad=(1/60)*(1/60)*(pi/180);
        deltaZ=arcSec2Rad*(16.271*Q*tan(z0).*(1+0.0000394*Q*tan(z0).^2)-0.0000749*P*(tan(z0)+tan(z0).^3));
        zTrue=z0+deltaZ;
        return;
    case 2%The IAU formula.....
    	CELMIN = 1e-6;
        SELMIN = 0.05;
        
        %Get the refraction parameters.
        [A,B]=simpAstroRefParam(Rh,P,T,wl);
        
        %Use the same bounding on the r and z values as in the iauAtioq
        %function.
        r=max(sin(z0),CELMIN);
        z=max(cos(z0),SELMIN);
        tanZ=r./z;
        
        w=B*tanZ.*tanZ;
        
        %The formula used in the iauAtioq function.
        %It is labeled the A*tan(z)+B*tan^3(z) model, with Newton-Raphson
        %correction.
        deltaZ=(A+w).*tanZ./(1+(A+3*w)./(z.*z));
        zTrue=z0+deltaZ;
    otherwise
        error('An invalid value for the algorithm was provided');
end
    function val=integrand(z,z0,refracFunc)
    %This function is the argument of the integral for algorithm 0.
    
            %First, find the r value associated with the z value.
            rVal=zenithDist2r(refracFunc,n0,r0,z0,z);
            
            %Next, find the index of refraction and its derivative with
            %respect to altitude.
            [n,dndr]=refracFunc(rVal-re);
            
            %Equation 7.87 in [1].
            val=rVal.*dndr./(n+rVal.*dndr);
    end

    function r=zenithDist2r(nRefrac,n0,r0,z0,z)
    %This function finds a zenith distance corresponding to a given radial
    %distance from the center fo the Earth. It is needed for algorithm 0
    %and implements Equation 7.84 of [1].
    
        r=r0;
        %Assume convergence in six iterations.
        numIter=0;
        while numIter<6
            [n,dndr]=nRefrac(r-re);

            r=r-(n.*r-n0*r0*sin(z0)./sin(z))./(n+r.*dndr);
            numIter=numIter+1;
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
