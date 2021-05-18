function [V,gradV,HessianV]=spherHarmonicEval(C,S,point,a,c,systemType,spherDerivs,scalFactor,algorithm)
%%SPHERHARMONICEVAL Evaluate a real or complex potential (e.g.
%                   gravitational or magnetic field) and/ or the gradient
%                   and Hessian of a potential when the potential is
%                   expressed in terms of real or complex spherical
%                   harmonic coefficients at a particular point given in
%                   spherical coordinates. Alternatively, this function can
%                   be used to evaluate the type of spherical harmonic
%                   series used to express terrain heights when given just
%                   spherical azimuth and elevation (but not range).
%
%INPUTS: C A length (M+2)*(M+1)/2 array holding the real or complex
%          coefficient terms that are multiplied by cosines in the harmonic
%          expansion. The coefficients must be fully normalized using the
%          type of full normalization that is used in the EGM2008 model.
%          Their normalization type can be changed using the
%          changeSpherHarmonicNorm function If given to a
%          CountingClusterSet class, C(n+1,m+1) is the coefficient of
%          degree n and order m. When a maximum degree of M is used, all C
%          must have values for all n from 0 to M and for all m from 0 to n
%          for each n. If coefficients are not present for certain degrees,
%          then insert a 0. It is assumed that M>=3. For magnetic field
%          models, the units of C should generally be Tesla. For
%          gravitational models, C is generally unitless with the units of
%          the function output being determined entirely by the parameters
%          a and c.
%        S A length (M+2)*(M+1)/2 array holding the coefficient terms that
%          are multiplied by sines in the harmonic expansion. The
%          requirements on S are the same as those on C.
%    point The 3XN set of N real points at which the potential and/or
%          gradient should be evaluated given in SPHERICAL coordinates
%          consisting of [r;azimuth;elevation]; When evaluating points on a
%          grid, the algorithm will be fastest if the points are provided
%          presorted by range and then by azimuth. This reduces the amount
%          of recomputation of certain values. Alternatively, if C and S
%          are for evaluating terrain heights, then points are 2XN having
%          the format [azimuth;elevation] and it is best if the points are
%          sorted by azimuth. Azimuth is measured counterclockwise from the
%          x-axis in the x-y plane. Elevation is measured up from the x-y
%          plane (towards the z-axis). Angles must be given in radians.
%        a The real or complex numerator in the (a/r)^n term in the
%          spherical harmonic sum. Normally, this is some type of a
%          reference radius. For example, when using most gravitational
%          models, a is the semi-major axis of the reference ellipsoid. If
%          this parameter is omitted, it is assumed that one is using the
%          spherical harmonics with something like the National Geospatial
%          Intelligence Agency's (NGA's) EGM96 or EGM2008 models, in which
%          case a=Constants.EGM2008SemiMajorAxis is used unless point is
%          2D, in which case c=1 is used.
%        c The real or complex constant value by which the spherical
%          harmonic series is multiplied. For example, for gravitational
%          potentials, the value is usually GM where G is the universal
%          gravitational constant and  M is the mass of the Earth. When
%          using the International Geomagnetic Reference Field (IGRF), the
%          constant is a^2, where a is the same as the numerator in the a/r
%          term. If this parameter is omitted, it is assumed that one is
%          using the spherical harmonics with something like the NGA's
%          EGM96 or EGM2008 models, in which case c=Constants.EGM2008GM is
%          used unless point is 2D, in which case c=1 is used.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
% spherDerivs This parameter specifies whether the gradient and Hessian
%          terms should be returned in spherical coordinates rather than
%          Cartesian coordinates. Possible values are true and false. The
%          default if this parameter is omitted or an empty matrix is
%          passed is false. The spherical coordinate system used by the
%          gradient and Hessian will match that specified by the systemType
%          input.
% scalFactor An optional real scale factor used in computing the normalized
%          associated Legendre polynomials. Generally, the default value
%          (if the scalFactor parameter is omitted) of 10^(-280) that is
%          suggested in [1] is sufficient. When very high-order models are
%          used, this scale factor prevents overflows.
% algorithm An optional parameter that selects which algorithm to use.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) If only
%            the potential is desired, always use Legendre's algorithm.
%            Otherwise, for points above 88 degrees spherical latitude, use
%            Pines' algorithm and for all other points use Legendre's
%            algorithm.
%          1 Only use Legendre's algorithm regardless of the location of
%            the points. Note that the gradient is singular at the poles,
%            so numerical problems will arise.
%          2 Only use Pines' algorithm.
%
%OUTPUTS: V The NX1 set of scalar potentials as obtained from the spherical
%           harmonic series. When dealing with gravitational models, the
%           the SI units are m^2/s^2, when dealing with geomagnetic models,
%           the SI units are T*m. This is real is all of the inputs are
%           real.
%     gradV The 3XN set of gradients of the potential in Cartesian
%           coordinates as obtained using the spherical harmonic series.
%           The derivatives are in the order [dV/dx;dV/dy;dV/dz] for
%           Cartesian values and in the order [dV/dr;dV/dAz;dV/dEl] for
%           spherical values. If V is the gravitational potential of the
%           Earth, then gradV is the acceleration due to the Earth's
%           gravity. If the point in question is on the rotating Earth,
%           then the centrifugal force and the Coriolis force must be taken
%           into account to get the acceleration experienced by the
%           observer. If V is a magnetic potential, then gradV is the
%           negative of the magnetic flux density vector B.
%     HessV The 3X3XN set of Hessian matrices of the potential.  If
%           Cartesian derivatives are used, then the ordering is
%           [d2/(dxdx),d2/(dxdy),d2/(dxdz);
%            d2/(dydx),d2/(dydy),d2/(dydz);
%            d2/(dzdx),d2/(dzdy),d2/(dxdx)]; If spherical derivatives are
%            used, then the ordering is the same with (x,y,z) replaced by
%            (r,Az,El).
%
%This function implements Legendre's algorithm from [2], but evaluates the
%sums using Horner's method as in and the Legendre function ratios as
%described in [1] and [6]. This function also implements Pines' algorithm
%is from [2]. Note that [2] contains a number of errors, many of which are
%corrected in [4]. A discussion on the uses of spherical harmonics is given
%in [3].
%
%When non-normalized coefficients are used, the spherical harmonic series
%for the potential is assumed to be of the form
%V=(c/r)*sum_{n=0}^M\sum_{m=0}^n(a/r)^n*(C(n+1,m+1)*cos(m*lambda)+S(n+1,m+1)*sin(m*lambda))*P^m_n(sin(theta))
%where r, lambda, and theta are the range, azimuth and elevation in a
%spherical coordinate system and P^m_n(x) is the notation for an associated
%Legendre function of degree n and order m evaluated at the point x. In the
%above sum, the C and S terms are as in the inputs, except for the shown
%dual indexation to be used, they are assumed to have been passed to the 
%changeSpherHarmonicNorm function. Note that theta is essentially a
%latitude, not a co-latitude, as is commonly used in an alternate
%definition of spherical coordinates. When computing terrain heights,
%which are specified only by spherical latitude and longitude, the (c/r)
%term as well as the (a/r) term disappear. When fully normalized
%coefficients as in the EGM2008 model are used, which is how this function
%is implemented, then P^m_n(x) is replaced with fully normalized associated
%Legendre functions \bar{P}_{nm}(x). Note that magnetic field models often
%use Schmidt semi-normalized coefficients, in which case P^m_n(x) is
%replaced with S^m_n(x), the Schmidt semi-normalized Legendre functions.
%The function changeSpherHarmonicNorm can be used to change the
%normalization of spherical harmonic coefficients so that they are suitable
%for use in this function.
%
%Note that models will often separate out the n=0 term, since
%sin(0*lambda)=0, cos(0*lambda)=1 and P^0_0(x)=1 and \bar{P}_{00}(x)=1 for
%any x. In such an instance, a C(0+1,0+1) might need to be explicitly
%added. For example, when using the coefficients for the EGM96 and EGM2008
%gravitational models provided by the National Geospatial Intelligence
%Agency (NGA), one has to insert C(0+1,0+1)=1.
%
%To understand the type of spherical harmonics that need to be inputted,
%one must disambiguate the different types of associated Legendre functions
%that can be used in the series.
%
%With the notation P^m_n(x) for an associated Legendre function of degree n
%and order m evaluated at the point x, one generally means 
%P^m_n(x)=(-1)^m*(1-x^2)^(m/2)*D_m{P_n(x)}
%where P_n(x) is a Legendre polynomial of degree n evaluated at x and
%D_m{P_n(x)} is the mth derivative of the polynomial. These associated
%Legendre functions are the same as those used in the International
%Geomagnetic Reference Field. On the other hand, with the notation
%P_{nm}(x), one generally means (-1)^mP^m_n(x). This notation is also
%called an associated Legendre function. Matlab's built-in function
%legendre(n,x) returns all P^m_n(x) for m=0 to m=n. On the other hand, the
%notation \bar{P}_{nm}(x), refers to a fully normalized associated Legendre
%function. Multiple definitions of what is normalized exist. In this
%function, the normalization is such that
%\bar{P}_{nm}(x)=P_{nm}(x)*sqrt(k*(2*n+1)*factorial(n-m)/factorial(n+m))
%where k=1 if m=0 and k=2 otherwise. This differs from the normalization
%constant that Matlab uses in the normalized version of its built-in
%associated Legendre function. It is the same as the normalization 
%constant the NGA uses in its coefficients for the EGM96 and the
%EGM2008 gravitation models. When considering Schmidt semi-normalized
%associated Legednre functions, the relation is 
%S^m_n(x)=P_{nm}(x)*sqrt(2*factorial(n-m)/factorial(n+m))*(-1)^m
%Fully normalized associated Legendre functions are most commonly used with
%geomagnetic models. Schmidt semi-normalized associated Legendre functions
%are most commonly used with geomagnetic models.
%
%In electromagnetics, one often uses complex spherical harmonics. The
%complex spherical harmonic model that is usually used in electromagnetics is 
%V=(c/r)*sum_{n=0}^M\sum_{m=-n}^n(a/r)^n*A(n,m)*\bar{P}_{nm}(sin(theta))*exp(1j*m*lambda)
%Initially, it just looks like one has combined the S and C values into one
%complex coefficient A and then multiplied the sin term by 1j. However, if
%that were all that had been done, one wouldn't be able to match many
%models, because if A(n,m)=aR+1j*aI, expanding the term
%A(n,m)*exp(1j*m*lambda), one gets 
%aR*cos(m*lambda)-aI*sin(m*lambda)+1j*(aI*cos(m*lambda)+aR*sin(m*lambda))
%Where one can see that it is not possible to independently set the real
%and the imaginary terms. However, note that the formulation of the inner
%sum used for complex values starts at -n and not at 0. Expand the
%coefficient A(n,-m)=aR1+1j*aI1. Again, expanding the term,
%A(n,-m)*exp(-1j*m*lambda), one gets
%aR1*cos(m*lambda)+aI1*sin(m*lambda)+1j*(aI1*cos(m*lambda)-aR1*sin(m*lambda))
%The values in A(n,-m) and A(n,m) can be combined to get
%(aR+aR1)*cos(m*lambda)+(aI1-aI)*sin(m*lambda)+1j*((aI+aI1)*cos(m*lambda)+(aR-aR1)*sin(m*lambda))
%Consequently, evaluating a complex spherical harmonic sum is equivalent to
%evaluating two real spherical harmonic sums, where
%C1(n,m)=aR+aR1, S1(n,m)=aI1-aI and C2(m,n)=aI+aI1, S2(m,n)=aR-aR1
%The first sum corresponds to the real part of the solution and the second
%sum corresponds to the imaginary part of the solution. This function can
%handle complex spherical harmonics directly, more efficiently than calling
%the function twice. If given complex coefficients for sums starting from
%-m, then the complex coefficients to give this function are C=C1+1j*C2 and
%S=S1*1j*S2.
%
%When spherical harmonic potentials and gradients are desired to a high
%degree and order, one can easily run into the limitations of double
%precision floating point numbers. For that reason, spherical harmonic
%series are often expressed in terms of fully normalized associated
%Legendre functions, rather than the standard Legendre functions that are
%more commonly used in textbooks. However, when the degree and order
%becomes very large, one can no longer directly compute even the fully
%normalized associated Legendre functions at all points, due to precision
%limitations. This means that one can not simply use the built-in
%legendre function in Matlab to obtain the polynomials. Consequently, when
%fully normalized coefficients are passed, the modified forward row (MFR)
%algorithm of [1] is used in this function to avoid precision problems.
%Note that the notation of the normalized Legendre functions in the above
%paper is slightly different than that used here. First, the paper uses a
%colatitude (Let's call it \tilde{\theta}), not a latitude theta
%(elevation) in its spherical coordinate system. Thus, the sine term for
%the fully normalized associated Legendre function in the
%sum becomes a cosine term. That is, \bar{P}_{nm}(cos(\tilde{\theta})).
%However, the paper abbreviates the notation and does not explicitly write
%the cosine operator Thus, one might be led to erroneously believe that the
%paper is evaluating \bar{P}_{nm}(\tilde{\theta}).
%
%The Legendre algorithm of [2], with the modifications to user Horner's
%method as in [1], cannot be used for computing the acceleration due to
%gravity near the poles, because of singularities in the derivatives of
%spherical coordinates. Thus, the method of Pines, as expressed using fully
%-normalized Helmholtz polynomials as described in [2] is the default used
%instead when the latitude is within 2 degrees of the pole and the user
%wants gradV and/or HessianV.
%
%This validity of this function can be tested on the surface of the
%reference ellipsoid using the coefficients for the ellipsoid. That is, at
%any point given in ellipsoidal coordinates with zero height, the
%gravity potential should be a constant and equal that given by the
%ellipParam2Potential function. For example,
%
% M=50;
% lambda=-20*pi/180;
% theta=80*pi/180;%Ellipsoidal latitude
% [C,S]=ellipsGravCoeffs(M);
% %The zero altitude places the point on the reference ellipse.
% pointEllips=[theta;lambda;0];
% point=ellips2Sphere(pointEllips);
% CartPoint=ellips2Cart(pointEllips);
% a=Constants.WGS84SemiMajorAxis;
% c=Constants.WGS84GMWithAtmosphere;
% omega=Constants.WGS84EarthRotationRate;
% [V,gradV]=spherHarmonicEval(C,S,point,a,c);
%
%The value 
% U=V+0.5*omega^2*(CartPoint(1)^2+CartPoint(2)^2);
%should be a constant for all values of lambda and theta above (withing
%precision limits for a 20th-degree series), as they are all points on the
%surface of the reference ellipsoid.
%Also, the value
% g=gradV+omega^2*[CartPoint(1);CartPoint(2);0];
%is the acceleration due to gravity including the centrifugal acceleration
%due to the rotation of the Earth. The values U and g should (within
%precision bounds for a 20th-degree series) be the same as the values of U
%and g returned by
% [U,g]=ellipsParam2Grav(CartPoint);
%
%More information on spherical harmonic representations of the Earth's
%gravitational field is given in [5].
%
%It is recommended that the C++ version of this function be compiled and
%used in place of this as it is thousands of times faster when C and S are
%very large.
%
%EXAMPLE 1:
%Here we validate the analytic gradient and Hessian values here by
%comparing them to numerically computed values when using a magnetic model.
% [C,S,a,c]=getWMMCoeffs();
% xCart=ellips2Cart([0.1;-0.1;0]);
% xCart=xCart/norm(xCart);
% systemType=2;
% xSpher=Cart2Sphere(xCart,systemType);
% 
% [val,gradV,HessianV]=spherHarmonicEval(C,S,xSpher,a,c,systemType);
% epsVal=1e-8;
% xCart1=xCart+[epsVal;0;0];
% xSpher1=Cart2Sphere(xCart1,systemType);
% [valdX,gradVdX]=spherHarmonicEval(C,S,xSpher1,a,c,systemType);
% xCart1=xCart+[0;epsVal;0];
% xSpher1=Cart2Sphere(xCart1,systemType);
% [valdY,gradVdY]=spherHarmonicEval(C,S,xSpher1,a,c,systemType);
% xCart1=xCart+[0;0;epsVal];
% xSpher1=Cart2Sphere(xCart1,systemType);
% [valdZ,gradVdZ]=spherHarmonicEval(C,S,xSpher1,a,c,systemType);
% gradNumDiff=zeros(3,1);
% gradNumDiff(1)=(valdX-val)/epsVal;
% gradNumDiff(2)=(valdY-val)/epsVal;
% gradNumDiff(3)=(valdZ-val)/epsVal;
% HessianNumDiff=zeros(3,3);
% HessianNumDiff(1,1)=(gradVdX(1)-gradV(1))/epsVal;
% HessianNumDiff(1,2)=(gradVdX(2)-gradV(2))/epsVal;
% HessianNumDiff(2,1)=HessianNumDiff(1,2);
% HessianNumDiff(1,3)=(gradVdX(3)-gradV(3))/epsVal;
% HessianNumDiff(3,1)=HessianNumDiff(1,3);
% HessianNumDiff(2,2)=(gradVdY(2)-gradV(2))/epsVal;
% HessianNumDiff(2,3)=(gradVdY(3)-gradV(3))/epsVal;
% HessianNumDiff(3,2)=HessianNumDiff(2,3);
% HessianNumDiff(3,3)=(gradVdZ(3)-gradV(3))/epsVal;
% max(abs((gradNumDiff(:)-gradV(:))./gradV(:)))
% max(abs((HessianNumDiff(:)-HessianV(:))./HessianV(:)))
%One will see agreement in the numerically computed gradient and Hessian
%better than 1e-6.
%
%EXAMPLE 2:
%This is the same as example 1, except we are using spherical derivatives.
% [C,S,a,c]=getWMMCoeffs();
% xCart=ellips2Cart([0.1;-0.1;0]);
% xCart=xCart/norm(xCart);
% systemType=0;
% spherDerivs=true;
% xSpher=Cart2Sphere(xCart,systemType);
% 
% [val,gradV,HessianV]=spherHarmonicEval(C,S,xSpher,a,c,systemType,spherDerivs);
% epsVal=1e-8;
% xSpher1=xSpher+[epsVal;0;0];
% [valdX,gradVdX]=spherHarmonicEval(C,S,xSpher1,a,c,systemType,spherDerivs);
% xSpher1=xSpher+[0;epsVal;0];
% [valdY,gradVdY]=spherHarmonicEval(C,S,xSpher1,a,c,systemType,spherDerivs);
% xSpher1=xSpher+[0;0;epsVal];
% [valdZ,gradVdZ]=spherHarmonicEval(C,S,xSpher1,a,c,systemType,spherDerivs);
% gradNumDiff=zeros(3,1);
% gradNumDiff(1)=(valdX-val)/epsVal;
% gradNumDiff(2)=(valdY-val)/epsVal;
% gradNumDiff(3)=(valdZ-val)/epsVal;
% HessianNumDiff=zeros(3,3);
% HessianNumDiff(1,1)=(gradVdX(1)-gradV(1))/epsVal;
% HessianNumDiff(1,2)=(gradVdX(2)-gradV(2))/epsVal;
% HessianNumDiff(2,1)=HessianNumDiff(1,2);
% HessianNumDiff(1,3)=(gradVdX(3)-gradV(3))/epsVal;
% HessianNumDiff(3,1)=HessianNumDiff(1,3);
% HessianNumDiff(2,2)=(gradVdY(2)-gradV(2))/epsVal;
% HessianNumDiff(2,3)=(gradVdY(3)-gradV(3))/epsVal;
% HessianNumDiff(3,2)=HessianNumDiff(2,3);
% HessianNumDiff(3,3)=(gradVdZ(3)-gradV(3))/epsVal;
% max(abs((gradNumDiff(:)-gradV(:))./gradV(:)))
% max(abs((HessianNumDiff(:)-HessianV(:))./HessianV(:)))
%Again, with spherical coordinates, one will see that the numerical and
%analytic valuesa gree to better than 1e-6.
%
%REFERENCES:
%[1] S. A. Holmes and W. E. Featherstone, "A unified approach to the
%    Clenshaw summation and the recursive computation of very high degree
%    and order normalised associated Legendre functions," Journal of
%    Geodesy, vol. 76, no. 5, pp. 279-299, May 2002.
%[2] E. Fantino and S. Casotto, "Methods of harmonic synthesis for global
%    geopotential models and their first-, second- and third-order
%    gradients," Journal of Geodesy, vol. 83, no. 7, pp. 595-619, Jul.
%    2009.
%[3] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%[4] E. Fantino and S. Casotto, "Erratum to: Methods of harmonic synthesis
%    for global geopotential models and their first-, second- and third-
%    order gradients," Journal of Geodesy, vol. 89, no. 11, pp. 1159-1161,
%    Nov. 2015.
%[5] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%[6] S. A. Holmes and W. E. Featherstone, "Short note: Extending simplified
%    high-degree synthesis methods to second latitude derivatives of
%    geopotential," Journal of Geodesy, vol. 76, no. 8, pp. 447-450, Nov.
%    2002.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(algorithm))
    algorithm=0; 
end

if(nargin<8||isempty(scalFactor))
    scalFactor=10^(-280);
end

if(nargin<7||isempty(spherDerivs))
    spherDerivs=false; 
end

if(nargin<6||isempty(systemType))
    systemType=0; 
end

if(nargin<5||isempty(c))
    if(size(point,1)==2)
        c=1;
    else
        c=Constants.EGM2008GM;
    end
end

if(nargin<4||isempty(a))
    if(size(point,1)==2)
        a=1;
    else
        a=Constants.EGM2008SemiMajorAxis;
    end
end

if(systemType~=0&&systemType~=2)
    error('An unsupported systemType was specified.')
end

M=(1/2)*(sqrt(1+8*length(C))-1)-1;

if(M<3)
    error('The coefficients must be provided to at least degree 3. To use a lower degree, one can insert zero coefficients.');
end

numPoints=size(point,2);
%If we are evaluating terrain heights.
switch(size(point,1))
    case 2
        point=[ones(1,numPoints);point(1,:);point(2,:)];
    case 3
    otherwise
        error('Invalid point length');
end

%Using a CountingClusterSet simplifies the indexation of the coefficients.
C=CountingClusterSet(C);
S=CountingClusterSet(S);

%Preallocate space used by the modified forward row algorithm when
%evaluating over multiple values with the same range and latitude but
%different longitudes.
if(algorithm==0||algorithm==1)
    A=zeros(M+1,1);
    B=zeros(M+1,1);
    if(nargout>1)
        Ar=zeros(M+1,1);
        Br=zeros(M+1,1);
        ATheta=zeros(M+1,1);
        BTheta=zeros(M+1,1);

        if(nargout>2)
            AThetaTheta=zeros(M+1,1);
            BThetaTheta=zeros(M+1,1);

            Arr=zeros(M+1,1);
            Brr=zeros(M+1,1);

            AThetar=zeros(M+1,1);
            BThetar=zeros(M+1,1);
        end
    end
end

V=zeros(numPoints,1);
gradV=zeros(3,numPoints);
HessianV=zeros(3,3,numPoints);

rPrev=Inf;
thetaPrev=Inf;
for curPoint=1:numPoints
    pointCur=point(:,curPoint);
    
    if(systemType==2)
        pointCur(3)=pi/2-pointCur(3); 
    end
    
    r=pointCur(1);
    lambda=pointCur(2);
    thetaCur=pointCur(3);
    
    rChanged=r~=rPrev;
    thetaChanged=thetaCur~=thetaPrev;
    rPrev=r;
    thetaPrev=thetaCur;

    if(rChanged)
        crScal=(c/r)/scalFactor;
        
        %This stores all of the powers of a/r needed for the sum,
        %regardless of which algorithm is used.
        nCoeff=zeros(M+1,1);
        nCoeff(1)=1;
        for n=1:M
            nCoeff(n+1)=nCoeff(n)*(a/r);
        end
    end
    
    switch(algorithm)
        case 0
        %At latitudes that are not near the poles, the Legendre method is
        %used. It cannot be used for the gradient or Hessian near the
        %poles, because of the singularity of the spherical coordinate
        %system.
            useLegendre=abs(thetaCur)<88*(pi/180)||nargout<2;
        case 1
            useLegendre=true;
        case 2
            useLegendre=false;
        otherwise
            error('Unknown Algorithm option specified.');
    end

    if(useLegendre)
        %Compute the sine and cosine terms.
        [SinVec,CosVec]=calcSinCosTerms(lambda,M);
        %The formulae for spherical harmonic synthesis with Legendre's
        %method uses clolatitude, pi/2-elevation
        theta=pi/2-thetaCur;
        u=sin(theta);
        if(thetaChanged)
            if(nargout==2)
                [PBarUVals,dPBarUValsdTheta]=NALegendreCosRat(theta,M,scalFactor);
            elseif(nargout>2)
                [PBarUVals,dPBarUValsdTheta,d2PBarUValsdTheta2]=NALegendreCosRat(theta,M,scalFactor);
            else
                PBarUVals=NALegendreCosRat(theta,M,scalFactor);
            end
        end

        %Evaluate Equation 7 from the Holmes and Featherstone paper.
        if(rChanged||thetaChanged)
            A(:)=0;
            B(:)=0;
            for m=0:M
                for n=m:M
                    CScal=nCoeff(n+1)*C(n+1,m+1);
                    SScal=nCoeff(n+1)*S(n+1,m+1);

                    %From Table 4
                    A(m+1)=A(m+1)+CScal*PBarUVals(n+1,m+1);
                    B(m+1)=B(m+1)+SScal*PBarUVals(n+1,m+1);
                end
            end
            
            %If additional terms should be computed so a gradient can
            %be determined.
            if(nargout>1)
                Ar(:)=0;
                Br(:)=0;
                ATheta(:)=0;
                BTheta(:)=0;
                
                for m=0:M
                    for n=m:M
                        CScal=nCoeff(n+1)*C(n+1,m+1);
                        SScal=nCoeff(n+1)*S(n+1,m+1);

                        %From Table 4
                        Ar(m+1)=Ar(m+1)+(n+1)*CScal*PBarUVals(n+1,m+1);
                        Br(m+1)=Br(m+1)+(n+1)*SScal*PBarUVals(n+1,m+1);

                        %From Table 4
                        ATheta(m+1)=ATheta(m+1)+CScal*dPBarUValsdTheta(n+1,m+1);
                        BTheta(m+1)=BTheta(m+1)+SScal*dPBarUValsdTheta(n+1,m+1);
                    end
                end
                
                %If additional terms should be computed so a Hessian can
                %be determined.
                if(nargout>2)
                    Arr(:)=0;
                    Brr(:)=0;
                    AThetar(:)=0;
                    BThetar(:)=0;
                    AThetaTheta(:)=0;
                    BThetaTheta(:)=0;
                    
                    for m=0:M
                        for n=m:M
                            CScal=nCoeff(n+1)*C(n+1,m+1);
                            SScal=nCoeff(n+1)*S(n+1,m+1);

                            %From Table 5, with the correction from the
                            %erratum.
                            Arr(m+1)=Arr(m+1)+(n+1)*(n+2)*CScal*PBarUVals(n+1,m+1);
                            Brr(m+1)=Brr(m+1)+(n+1)*(n+2)*SScal*PBarUVals(n+1,m+1);

                            %From Table 5
                            AThetar(m+1)=AThetar(m+1)+(n+1)*CScal*dPBarUValsdTheta(n+1,m+1);
                            BThetar(m+1)=BThetar(m+1)+(n+1)*SScal*dPBarUValsdTheta(n+1,m+1);

                            %From Table 5
                            AThetaTheta(m+1)=AThetaTheta(m+1)+CScal*d2PBarUValsdTheta2(n+1,m+1);
                            BThetaTheta(m+1)=BThetaTheta(m+1)+SScal*d2PBarUValsdTheta2(n+1,m+1);
                        end
                    end
                end
            end
        end

        %Use Horner's method to compute V all values.
        V(curPoint)=0;
        for m=M:-1:0
            V(curPoint)=V(curPoint)*u+(A(m+1)*CosVec(m+1)+B(m+1)*SinVec(m+1));
        end

        V(curPoint)=crScal*V(curPoint);

        if(nargout>1)
            %Use Horner's method to compute all dV values.
            dVdr=0;
            dVdLambda=0;
            dVdTheta=0;

            %The following first-order derivative formulae are from
            %Table 1 (expressed using Horner's method).
            for m=M:-1:0
                dVdr=dVdr*u-(Ar(m+1)*CosVec(m+1)+Br(m+1)*SinVec(m+1));
                dVdLambda=dVdLambda*u-m*(A(m+1)*SinVec(m+1)-B(m+1)*CosVec(m+1));
                dVdTheta=dVdTheta*u+(ATheta(m+1)*CosVec(m+1)+BTheta(m+1)*SinVec(m+1));
            end

            dVdr=crScal*dVdr/r;
            dVdLambda=crScal*dVdLambda;
            %The minus sign adjusts for the coordinate system change.
            dVdTheta=-crScal*dVdTheta;
            
            if(spherDerivs)
                %The sign change on the theta terms deals with the
                %different input coordinate system used.
                gradV(:,curPoint)=[dVdr;dVdLambda;dVdTheta];
            else%Convert the derivatives to Cartesian coordinates.
                J=calcSpherConvJacob(pointCur);
                gradV(:,curPoint)=J'*[dVdr;dVdLambda;dVdTheta];
            end
        end

        if(nargout>2)
            %Use Horner's method to compute d2V all values.
            d2VdLambdadLambda=0;
            d2VdLambdadTheta=0;
            d2VdrdLambda=0;
            d2VdThetadTheta=0;
            d2VdrdTheta=0;
            d2Vdrdr=0;

            %The following second-order derivative formulae are from
            %Table 2 (expressed using Horner's method).
            for m=M:-1:0
                d2VdLambdadLambda=d2VdLambdadLambda*u-m^2*(A(m+1)*CosVec(m+1)+B(m+1)*SinVec(m+1));
                d2VdLambdadTheta=d2VdLambdadTheta*u+m*(ATheta(m+1)*SinVec(m+1)-BTheta(m+1)*CosVec(m+1));
                d2VdrdLambda=d2VdrdLambda*u+m*(Ar(m+1)*SinVec(m+1)-Br(m+1)*CosVec(m+1));
                d2VdThetadTheta=d2VdThetadTheta*u+(AThetaTheta(m+1)*CosVec(m+1)+BThetaTheta(m+1)*SinVec(m+1));
                d2VdrdTheta=d2VdrdTheta*u-(AThetar(m+1)*CosVec(m+1)+BThetar(m+1)*SinVec(m+1));
                d2Vdrdr=d2Vdrdr*u+(Arr(m+1)*CosVec(m+1)+Brr(m+1)*SinVec(m+1));
            end
            
            %The minus signs in the following equations adjust for the
            %spherical coordinate system difference.
            d2Vdrdr=crScal*d2Vdrdr/r^2;
            d2VdLambdadLambda=crScal*d2VdLambdadLambda;
            d2VdThetadTheta=crScal*d2VdThetadTheta;
            d2VdrdLambda=crScal*d2VdrdLambda/r;
            d2VdrdTheta=-crScal*d2VdrdTheta/r;
            d2VdLambdadTheta=crScal*d2VdLambdadTheta;

            if(spherDerivs)
                HessianV(1,1,curPoint)=d2Vdrdr;
                HessianV(2,2,curPoint)=d2VdLambdadLambda;
                HessianV(3,3,curPoint)=d2VdThetadTheta;
                HessianV(1,2,curPoint)=d2VdrdLambda;
                HessianV(2,1,curPoint)=HessianV(1,2,curPoint);
                HessianV(1,3,curPoint)=d2VdrdTheta;
                HessianV(3,1,curPoint)=HessianV(1,3,curPoint);
                HessianV(2,3,curPoint)=d2VdLambdadTheta;
                HessianV(3,2,curPoint)=HessianV(2,3,curPoint);
            else%Convert the Hessian to Cartesian coordinates.
                drdx=J(1,1);
                drdy=J(1,2);
                drdz=J(1,3);
                dLambdadx=J(2,1);
                dLambdady=J(2,2);
                dLambdadz=J(2,3);
                dPhidx=J(3,1);
                dPhidy=J(3,2);
                dPhidz=J(3,3);

                H=calcSpherConvHessian(pointCur);

                drdxdx=H(1,1,1);
                drdydy=H(2,2,1);
                drdzdz=H(3,3,1);
                drdxdy=H(1,2,1);
                drdxdz=H(1,3,1);
                drdydz=H(2,3,1);

                dLambdadxdx=H(1,1,2);
                dLambdadydy=H(2,2,2);
                dLambdadzdz=H(3,3,2);
                dLambdadxdy=H(1,2,2);
                dLambdadxdz=H(1,3,2);
                dLambdadydz=H(2,3,2);

                dPhidxdx=H(1,1,3);
                dPhidydy=H(2,2,3);
                dPhidzdz=H(3,3,3);
                dPhidxdy=H(1,2,3);
                dPhidxdz=H(1,3,3);
                dPhidydz=H(2,3,3);

                HessianV(1,1,curPoint)=d2VdLambdadLambda*dLambdadx^2+2*d2VdLambdadTheta*dLambdadx*dPhidx+d2VdThetadTheta*dPhidx^2+2*dPhidx*drdx*d2VdrdTheta+2*dLambdadx*drdx*d2VdrdLambda+dVdLambda*dLambdadxdx+dVdTheta*dPhidxdx+dVdr*drdxdx+drdx^2*d2Vdrdr;
                HessianV(2,2,curPoint)=d2VdThetadTheta*dPhidy^2+2*dLambdady*dPhidy*d2VdLambdadTheta+dVdLambda*dLambdadydy+dVdTheta*dPhidydy+dLambdady^2*d2VdLambdadLambda+drdydy*dVdr+2*dPhidy*drdy*d2VdrdTheta+2*dLambdady*drdy*d2VdrdLambda+drdy^2*d2Vdrdr;
                HessianV(3,3,curPoint)=dVdTheta*dPhidzdz+dPhidz^2*d2VdThetadTheta+dLambdadzdz*dVdLambda+2*dLambdadz*dPhidz*d2VdLambdadTheta+(dLambdadz)^2*d2VdLambdadLambda+drdzdz*dVdr+2*dPhidz*drdz*d2VdrdTheta+2*dLambdadz*drdz*d2VdrdLambda+drdz^2*d2Vdrdr;
                HessianV(1,2,curPoint)=dPhidy*d2VdLambdadTheta*dLambdadx+dLambdady*d2VdLambdadLambda*dLambdadx+d2VdThetadTheta*dPhidy*dPhidx+dLambdady*d2VdLambdadTheta*dPhidx+drdy*dPhidx*d2VdrdTheta+dPhidy*drdx*d2VdrdTheta+dVdLambda*dLambdadxdy+dVdTheta*dPhidxdy+dVdr*drdxdy+drdy*dLambdadx*d2VdrdLambda+dLambdady*drdx*d2VdrdLambda+drdy*drdx*d2Vdrdr;
                HessianV(2,1,curPoint)=HessianV(1,2,curPoint);
                HessianV(1,3,curPoint)=dPhidz*d2VdLambdadTheta*dLambdadx+dLambdadz*d2VdLambdadLambda*dLambdadx+dPhidz*d2VdThetadTheta*dPhidx+dLambdadz*d2VdLambdadTheta*dPhidx+dVdLambda*dLambdadxdz+dVdTheta*dPhidxdz+dVdr*drdxdz+drdz*dPhidx*d2VdrdTheta+dPhidz*drdx*d2VdrdTheta+drdz*dLambdadx*d2VdrdLambda+dLambdadz*drdx*d2VdrdLambda+drdz*drdx*d2Vdrdr;
                HessianV(3,1,curPoint)=HessianV(1,3,curPoint);
                HessianV(2,3,curPoint)=dPhidz*d2VdThetadTheta*dPhidy+dVdLambda*dLambdadydz+dVdTheta*dPhidydz+dPhidz*dLambdady*d2VdLambdadTheta+dLambdadz*dPhidy*d2VdLambdadTheta+dLambdadz*dLambdady*d2VdLambdadLambda+drdydz*dVdr+drdz*dPhidy*d2VdrdTheta+dPhidz*drdy*d2VdrdTheta+drdz*dLambdady*d2VdrdLambda+dLambdadz*drdy*d2VdrdLambda+drdz*drdy*d2Vdrdr;
                HessianV(3,2,curPoint)=HessianV(2,3,curPoint);
            end
        end
    else
        %At latitudes that are near the poles, the non-singular algorithm
        %of Pines using the fully normalized Helmholtz equations from
        %Fantino and Casotto is used. The algorithm has been slightly
        %modified so that the c/r term is out front and the fully
        %normalized Helmholtz polynomials can be scaled. Also, lumped
        %coefficients are not used. The Pines algorithm is generally slower
        %than the algorithm of Holmes and Featherstone and it suffers a
        %loss of precision near the equator. Thus, the Pines algorithm is
        %only used near the poles where the other algorithm has issues with
        %a singularity.
        
        CartPoint=spher2Cart(pointCur);

        x=CartPoint(1);
        y=CartPoint(2);
        z=CartPoint(3);

        %Get the direction cosines used by Pines' algorithm.
        s=x/r;
        t=y/r;
        u=z/r;

        %Compute the fully normalized Helmholtz polynomials.
        if(thetaChanged)
            if(nargout>1)
                [HBar,dHBardu,d2HBardu2]=normHelmholtz(u,M,scalFactor);
            else
                HBar=normHelmholtz(u,M,scalFactor);
            end
        end

        %Recursively compute the rm and im terms for the sums.
        rm=zeros(M+1,1);
        im=zeros(M+1,1);
        rm(0+1)=1;
        im(0+1)=0;
        for m=1:M
            %These are equation 49 in the Fantino and Casotto paper.
            rm(m+1)=s*rm(m-1+1)-t*im(m-1+1);
            im(m+1)=s*im(m-1+1)+t*rm(m-1+1);
        end

        %Perform the sum for the potential from Equation 44 in the Fantino
        %and Casotto paper.
        V(curPoint)=0;
        for n=0:M
            innerTerm=0;
            for m=0:n
                innerTerm=innerTerm+(C(n+1,m+1)*rm(m+1)+S(n+1,m+1)*im(m+1))*HBar(n+1,m+1);
            end
            V(curPoint)=V(curPoint)+nCoeff(n+1)*innerTerm;
        end

        V(curPoint)=crScal*V(curPoint);

        %If only the gradient is desired as the next output
        if(nargout==2)
            %If the gradient and not the Hessian is desired.
            a1=0;
            a2=0;
            a3=0;
            a4=0;
 
            for m=0:M
                A1=0;
                A2=0;
                A3=0;

                B1=0;
                B2=0;
                B3=0;

                %Compute the lumped coefficients for Pine's method from
                %Table 13 for the current m.
                for n=m:M
                    HVal=HBar(n+1,m+1);
                    dHVal=dHBardu(n+1,m+1);
                    
                    %The expressions for Lmn, is from Table 14
                    Lmn=(n+m+1)*HVal+u*dHVal;
                    
                    rhoC=nCoeff(n+1)*C(n+1,m+1);
                    rhoS=nCoeff(n+1)*S(n+1,m+1);

                    A1=A1+rhoC*HVal;
                    A2=A2+rhoC*dHVal;
                    A3=A3+rhoC*Lmn;

                    B1=B1+rhoS*HVal;
                    B2=B2+rhoS*dHVal;
                    B3=B3+rhoS*Lmn;
                end
                if(m>=1)
                    a1=a1+m*(A1*rm(m-1+1)+B1*im(m-1+1));
                    a2=a2+m*(B1*rm(m-1+1)-A1*im(m-1+1));
                end
                a3=a3+(A2*rm(m+1)+B2*im(m+1));
                a4=a4-(A3*rm(m+1)+B3*im(m+1));
            end
            
            a1=a1/r;
            a2=a2/r;
            a3=a3/r;
            a4=a4/r;
            
            dVdx=crScal*(a1+s*a4);
            dVdy=crScal*(a2+t*a4);
            dVdz=crScal*(a3+u*a4);
 
            if(spherDerivs)
                %Convert the derivatives to spherical coordinates.
                J=calcSpherInvJacob(pointCur)';

                gradV(:,curPoint)=J*[dVdx;dVdy;dVdz];
            else%If a gradient in Cartesian coordinates is desired.
                gradV(1,curPoint)=dVdx;
                gradV(2,curPoint)=dVdy;
                gradV(3,curPoint)=dVdz;
            end
        else
            %If the gradient and the Hessian are desired.
            a1=0;
            a2=0;
            a3=0;
            a4=0;
            a11=0;
            a12=0;
            a13=0;
            a14=0;
            a23=0;
            a24=0;
            a33=0;
            a34=0;
            a44=0;

            for m=0:M
                A1=0;
                A2=0;
                A3=0;
                A4=0;
                A5=0;
                A6=0;

                B1=0;
                B2=0;
                B3=0;
                B4=0;
                B5=0;
                B6=0;

                %Compute the lumped coefficients for Pine's method from
                %Table 13 for the current m.
                for n=m:M
                    HVal=HBar(n+1,m+1);
                    dHVal=dHBardu(n+1,m+1);
                    d2HVal=d2HBardu2(n+1,m+1);
                    
                    %The expressions for Lmn, dLmn, and Omn are from
                    %Table 14
                    Lmn=(n+m+1)*HVal+u*dHVal;
                    dLmn=(n+m+2)*dHVal+u*d2HVal;
                    Omn=(n+m+1)*(n+m+2)*HVal+2*u*(n+m+2)*dHVal+u^2*d2HVal;
                    
                    rhoC=nCoeff(n+1)*C(n+1,m+1);
                    rhoS=nCoeff(n+1)*S(n+1,m+1);

                    A1=A1+rhoC*HVal;
                    A2=A2+rhoC*dHVal;
                    A3=A3+rhoC*Lmn;
                    A4=A4+rhoC*d2HVal;
                    A5=A5+rhoC*dLmn;
                    A6=A6+rhoC*Omn;

                    B1=B1+rhoS*HVal;
                    B2=B2+rhoS*dHVal;
                    B3=B3+rhoS*Lmn;
                    B4=B4+rhoS*d2HVal;
                    B5=B5+rhoS*dLmn;
                    B6=B6+rhoS*Omn;
                end
                if(m>=1)
                    a1=a1+m*(A1*rm(m-1+1)+B1*im(m-1+1));
                    a2=a2+m*(B1*rm(m-1+1)-A1*im(m-1+1));
                end
                a3=a3+(A2*rm(m+1)+B2*im(m+1));
                a4=a4-(A3*rm(m+1)+B3*im(m+1));

                if(m>=2)
                    a11=a11+m*(m-1)*(A1*rm(m-2+1)+B1*im(m-2+1));
                    a12=a12+m*(m-1)*(B1*rm(m-2+1)-A1*im(m-2+1));
                end
                if(m>=1)
                    a13=a13+m*(A2*rm(m-1+1)+B2*im(m-1+1));
                    a14=a14-m*(A3*rm(m-1+1)+B3*im(m-1+1));
                    a23=a23+m*(B2*rm(m-1+1)-A2*im(m-1+1));
                    a24=a24-m*(B3*rm(m-1+1)-A3*im(m-1+1));
                end
                a33=a33+(A4*rm(m+1)+B4*im(m+1));
                a34=a34-(A5*rm(m+1)+B5*im(m+1));
                a44=a44+(A6*rm(m+1)+B6*im(m+1));
            end
            
            a1=a1/r;
            a2=a2/r;
            a3=a3/r;
            a4=a4/r;
            a11=a11/r^2;
            a12=a12/r^2;
            a13=a13/r^2;
            a14=a14/r^2;
            a23=a23/r^2;
            a24=a24/r^2;
            a33=a33/r^2;
            a34=a34/r^2;
            a44=a44/r^2;
            a22=-a11;
            
            crScal=(c/r)/scalFactor;
            
            dVdx=crScal*(a1+s*a4);
            dVdy=crScal*(a2+t*a4);
            dVdz=crScal*(a3+u*a4);
 
            if(spherDerivs)
                %Convert the derivatives to spherical coordinates.
                J=calcSpherInvJacob(pointCur)';

                gradV(:,curPoint)=J*[dVdx;dVdy;dVdz];
            else%If a gradient in Cartesian coordinates is desired.
                gradV(1,curPoint)=dVdx;
                gradV(2,curPoint)=dVdy;
                gradV(3,curPoint)=dVdz;
            end
            
            d2Vdxdx=crScal*(a11+2*s*a14+a4/r+s^2*a44-s^2*a4/r);
            d2Vdydy=crScal*(a22+2*t*a24+a4/r+t^2*a44-t^2*a4/r);
            d2Vdzdz=crScal*(a33+2*u*a34+a4/r+u^2*a44-u^2*a4/r);
            
            d2Vdxdy=crScal*(a12+s*t*a44+s*a24+t*a14-s*t*a4/r);
            d2Vdxdz=crScal*(a13+s*u*a44+s*a34+u*a14-s*u*a4/r);
            d2Vdydz=crScal*(a23+t*u*a44+t*a34+u*a24-t*u*a4/r);
            
            if(spherDerivs)
                dxdr=J(1,1);
                dxdAz=J(2,1);
                dxdEl=J(3,1);
                
                dydr=J(1,2);
                dydAz=J(2,2);
                dydEl=J(3,2);
                
                dzdr=J(1,3);
                dzdAz=J(2,3);
                dzdEl=J(3,3);

                H=calcSpherInvHessian(pointCur);
                
                d2xdrdr=H(1,1,1);
                d2xdAzdAz=H(2,2,1);
                d2xdEldEl=H(3,3,1);
                d2xdrdAz=H(1,2,1);
                d2xdrdEl=H(1,3,1);
                d2xdAzdEl=H(2,3,1);
                
                d2ydrdr=H(1,1,2);
                d2ydAzdAz=H(2,2,2);
                d2ydEldEl=H(3,3,2);
                d2ydrdAz=H(1,2,2);
                d2ydrdEl=H(1,3,2);
                d2ydAzdEl=H(2,3,2);
                
                d2zdrdr=H(1,1,3);
                d2zdAzdAz=H(2,2,3);
                d2zdEldEl=H(3,3,3);
                d2zdrdAz=H(1,2,3);
                d2zdrdEl=H(1,3,3);
                d2zdAzdEl=H(2,3,3);
                
                %d2Vdrdr
                HessianV(1,1,curPoint)=d2Vdydy*dydr^2+2*d2Vdydz*dydr*dzdr+d2Vdzdz*dzdr^2+2*dxdr*dzdr*d2Vdxdz+2*dxdr*dydr*d2Vdxdy+dxdr^2*d2Vdxdx+dVdx*d2xdrdr+dVdy*d2ydrdr+dVdz*d2zdrdr;
                %d2VdAzdAz
                HessianV(2,2,curPoint)=d2Vdzdz*dzdAz^2+2*dydAz*dzdAz*d2Vdydz+dydAz^2*d2Vdydy+dVdy*d2ydAzdAz+dVdz*d2zdAzdAz+d2xdAzdAz*dVdx+2*dxdAz*dzdAz*d2Vdxdz+2*dxdAz*dydAz*d2Vdxdy+dxdAz^2*d2Vdxdx;
                %d2VdEldEl
                HessianV(3,3,curPoint)=dzdEl^2*d2Vdzdz+dVdz*d2zdEldEl+d2ydEldEl*dVdy+2*dydEl*dzdEl*d2Vdydz+dydEl^2*d2Vdydy+d2xdEldEl*dVdx+2*dxdEl*dzdEl*d2Vdxdz+2*dxdEl*dydEl*d2Vdxdy+dxdEl^2*d2Vdxdx;
                %d2VdrdAz
                HessianV(1,2,curPoint)=dzdAz*d2Vdydz*dydr+dydAz*d2Vdydy*dydr+d2Vdzdz*dzdAz*dzdr+dydAz*d2Vdydz*dzdr+dzdAz*dxdr*d2Vdxdz+dxdAz*dzdr*d2Vdxdz+dydAz*dxdr*d2Vdxdy+dxdAz*dydr*d2Vdxdy+dVdx*d2xdrdAz+dVdy*d2ydrdAz+dVdz*d2zdrdAz+dxdAz*dxdr*d2Vdxdx;
                HessianV(2,1,curPoint)=HessianV(1,2,curPoint);
                %d2VdrdEl
                HessianV(1,3,curPoint)=dzdEl*d2Vdydz*dydr+dydEl*d2Vdydy*dydr+dzdEl*d2Vdzdz*dzdr+dydEl*d2Vdydz*dzdr+dzdEl*dxdr*d2Vdxdz+dxdEl*dzdr*d2Vdxdz+dVdx*d2xdrdEl+dVdy*d2ydrdEl+dVdz*d2zdrdEl+dydEl*dxdr*d2Vdxdy+dxdEl*dydr*d2Vdxdy+dxdEl*dxdr*d2Vdxdx;
                HessianV(3,1,curPoint)=HessianV(1,3,curPoint);
                %d2VdAzdEl
                HessianV(2,3,curPoint)=dzdEl*d2Vdzdz*dzdAz+dzdEl*dydAz*d2Vdydz+dydEl*dzdAz*d2Vdydz+dVdy*d2ydAzdEl+dVdz*d2zdAzdEl+dydEl*dydAz*d2Vdydy+d2xdAzdEl*dVdx+dzdEl*dxdAz*d2Vdxdz+dxdEl*dzdAz*d2Vdxdz+dydEl*dxdAz*d2Vdxdy+dxdEl*dydAz*d2Vdxdy+dxdEl*dxdAz*d2Vdxdx;
                HessianV(3,2,curPoint)=HessianV(2,3,curPoint);
            else
                HessianV(1,1,curPoint)=d2Vdxdx;
                HessianV(2,2,curPoint)=d2Vdydy;
                HessianV(3,3,curPoint)=d2Vdzdz;

                HessianV(1,2,curPoint)=d2Vdxdy;
                HessianV(2,1,curPoint)=HessianV(1,2,curPoint);
                HessianV(1,3,curPoint)=d2Vdxdz;
                HessianV(3,1,curPoint)=HessianV(1,3,curPoint);
                HessianV(2,3,curPoint)=d2Vdydz;
                HessianV(3,2,curPoint)=HessianV(2,3,curPoint);
            end
        end
    end
end

if(spherDerivs&&systemType==2)
    %Flip signs of the elevation terms reflecting the
    %difference in definition between systems 0 and 2.
    gradV(3,:)=-gradV(3,:);
    HessianV(1,3,:)=-HessianV(1,3,:);
    HessianV(3,1,:)=-HessianV(3,1,:);
    HessianV(2,3,:)=-HessianV(2,3,:);
    HessianV(3,2,:)=-HessianV(3,2,:);
end
end

function [SinVec,CosVec]=calcSinCosTerms(lambda,M)
    %Compute sin(m*lambda) and cos(m*lambda) for m=0 to m=M.
    SinVec=zeros(M+1,1);
    CosVec=zeros(M+1,1);
    %Explicitly set the first two terms.
    SinVec(0+1)=0;
    CosVec(0+1)=1;
    SinVec(1+1)=sin(lambda);
    CosVec(1+1)=cos(lambda);
    %Use a double angle identity to get the second order term.
    SinVec(2+1)=2*SinVec(1+1)*CosVec(1+1);
    CosVec(2+1)=1-2*SinVec(1+1)^2;
    %Use a two-part recursion for the rest of the terms.
    for m=3:M
        SinVec(m+1)=2*CosVec(1+1)*SinVec(m-1+1)-SinVec(m-2+1);
        CosVec(m+1)=2*CosVec(1+1)*CosVec(m-1+1)-CosVec(m-2+1);
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
