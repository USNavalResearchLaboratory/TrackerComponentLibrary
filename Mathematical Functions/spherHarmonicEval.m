function [V,gradV]=spherHarmonicEval(C,S,point,a,c,fullyNormalized,scalFactor)
%%SPHERHARMONICEVAL  Evaluate a potential (e.g. gravitational or magnetic)
%                    and/ or the gradient of a  potential when the
%                    potential is expressed in terms of spherical harmonic
%                    coefficients at a particular point given in spherical
%                    Earth-Centered Earth-Fixed coordinates. Alternatively,
%                    this function can be used to evaluate the type of
%                    spherical harmonic series used to express terrain
%                    heights when given just spherical azimuth and
%                    elevation (but no range).                   
%
%INPUTS:    C   A ClusterSet class holding the coefficient terms that are
%               multiplied by cosines in the harmonic expansion. C(n+1,m+1)
%               is the coefficient of degree n and order m. When a maximum
%               degree of M is used, all C must have values for all n
%               from 0 to M and for all m from 0 to n for each n. If
%               coefficients are not present for certain degrees, then
%               insert a 0. It is assumed that M>=3. For magnetic field
%               models, the units of C should generally be Tesla. For
%               gravitational models, C is generally unitless with the
%               units being determined entirely by the parameters a and c.
%           S   A ClusterSet class holding the coefficient terms that are
%               multiplied by sines in the harmonic expansion. The
%               requirements on S are the same as those on C.
%       point   The 3XN set of N points at which the potential and/or
%               gradient should be evaluated given in SPHERICAL, ECEF 
%               coordinates consisting of [r;azimuth;elevation]; When
%               evaluating points on a grid, the algorithm will be fastest
%               if the points are provided presorted by range and then by
%               azimuth. This reduces the amount of recomputation of
%               certain values. Alternatively, if C and S are for
%               evaluating terrain heights, then points are 2XN having the
%               format [azimuth;elevation] and it is best if the points are
%               sorted by azimuth.
%           a   The numerator in the (a/r)^n term in the spherical harmonic
%               sum. Normally, this is some type of a reference radius. For
%               example, when using most gravitational models, a is the
%               semi-major axis of the reference ellipsoid. If this
%               parameter is omitted, it is assumed that one is using the
%               spherical harmonics with something like the National
%               Geospatial Intelligence Agency's (NGA's) EGM96 or EGM2008
%               models, in which case a=Constants.EGM2008SemiMajorAxis is
%               used. Alternatively, when point is 2D for a terrain height,
%               this parameter is ignored (not necessary) and an empty
%               matrix can just be passed.
%           c   The constant value by which the spherical harmonic series
%               is multiplied. For example, for gravitational potentials,
%               the value is usually GM where G is the universal
%               gravitational constant and  M is the mass of the Earth.
%               When using the International Geomagnetic Reference Field,
%               the constant is a^2, where a is the same as the numerator
%               in the a/r term. If this parameter is omitted, it is
%               assumed that one is using the spherical harmonics with
%               something like the NGA's EGM96 or EGM2008 models, in which
%               case c=Constants.EGM2008GM is used. Alternatively, when
%               point is 2D for a terrain height, this parameter is ignored
%               (not necessary) and an empty matrix can just be passed.
%fullyNormalized A boolean variable indicating whether the coefficients are
%               fully normalized. If false, then it is assumed that the
%               coefficients are Schmidt semi-normalized. If this parameter
%               is omitted,  the default is true. The text below describes
%               the normalization. It is suggested that the coefficients be
%               normalized prior to passing them to this function as the
%               renormalization can be slow and can consume a lot of
%               memory.
%    scalFactor An optional scale factor used in computing the normalized
%               associated Legendre polynomials if fullyNormalized=true.
%               Generally, the default value (is if the scalFactor
%               parameter is omitted) of 10^(-280) that is suggested
%               in the Holmes and Featherstone paper (cited below) is
%               sufficient. When very high-order models are used, this
%               scale factor prevents overflows.
%
%OUTPUTS:  V    The potential as obtained from the spherical harmonic
%               series. When dealing with gravitational models, the SI
%               units are m^2/s^2, when dealing with geomagnetic models,
%               the SI units are T*m.
%      gradV    The  gradient of the potential in Cartesian coordinates as
%               obtained using the spherical harmonic series. If V is the
%               gravitational potential of the Earth, then gradV is the
%               acceleration due to the Earth's gravity. If the point in
%               question is on the rotating Earth, then the total
%               acceleration experienced by the (rotating) point due to
%               gravity and centrifugal force is gradV+omega^2*[x;y;0],
%               where x and y are the x and y Cartesian coordinates of the
%               point, and omega is the rotation rate in radians per second
%               (The rotation is assumed to be about the z-axis). If V is a
%               magnetic potential, then gradV is the negative of the
%               magnetic flux density vector B.
%
%When non-normalized coefficients are used, the spherical harmonic series
%for the potential is assumed to be of the form
%V=(c/r)*sum_{n=0}^M\sum_{m=0}^n(a/r)^n*(C(n+1,m+1)*cos(m*lambda)+S(n+1,m+1)*sin(m*lambda))*P^m_n(sin(theta))
%where r, lambda, and theta are the range, azimuth and elevation in a
%spherical coordinate system and P^m_n(x) is the notation for an associated
%Legendre function of degree n and order m evaluated at the point x. Note
%that theta is essentially a latitude, not a co-latitude, as is commonly
%used in an alternate definition of spherical coordinates. When computing
%terrain heights, which are specified only by spherical latitude and
%longitude, the (c/r) term as well as the (a/r) term disappear. When fully
%normalized coefficients are used (fullyNormalized=true), then P^m_n(x) is
%replaced with fully normalized associated Legendre functions
%\bar{P}_{nm}(x) and when Schmidt semi-normalized coefficients are used
%(fullyNormalized=false), then P^m_n(x) is replaced with S^m_n(x), the
%Schmidt semi-normalized legendre functions. This function does not
%explicitely handle coefficients that have no normalization, as they are
%not commonly used in practice.
%
%Note that models will often separate out the n=0 term, since
%sin(0*lambda)=0, cos(0*lambda)=1 and P^0_0(x)=1 and \bar{P}_{00}(x)=1 for
%any x. In such an instance, a C(0+1,0+1) might need to be explicitely
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
%constant the the NGA uses in its coefficients for the EGM96 and the
%EGM2008 gravitation models. When considering Schmidt semi-normalized
%associated Legednre functions, the relation is 
%S^m_n(x)=P_{nm}(x)*sqrt(2*factorial(n-m)/factorial(n+m))*(-1)^m
%Fully normalized associated Legendre functions are most commonly used with
%geomagnetic models. Schmidt semi-normalized associated Legendre functions
%are most commonly used with geomagnetic models.
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
%However, the paper abbreviates the notation and does not explicitely write
%the cosine operator Thus, one might be led to erroneously believe that the
%paper is evaluating \bar{P}_{nm}(\tilde{\theta}).
%
%The Holmes and Featherstone algorithm cannot be used for computing the
%acceleration due to gravity near the poles, because of singularities in
%the derivatives of spherical coordinates. Thus, the method of Pines, as
%expressed using fully-normalized Helmholtz polynomials as described in [2]
%is used instead when the latitude is within 2 degrees of the pole and the
%user wants gradV.
%
%On the other hand, when Schmidt-normalized coefficients are used, the
%coefficients are first converted to coefficients for fully-normalized
%associated Legendre functions. Note that the passed C and S values are
%duplicated and not modified during the conversion.
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
%gravitational field is given in 
%B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%SpringerWienNewYork, 2006.
%
%If the helper function spherHarmonicEvalCPPInt has been compiled, then
%after parsing the input, the algorithm will be run via that file (a C++
%implementation). This function is thousands of times faster when
%spherHarmonicEvalCPPInt is compiled, especially for high degrees and
%orders.
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
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(scalFactor))
    scalFactor=10^(-280);
end

if(nargin<6||isempty(fullyNormalized))
    fullyNormalized=true;
end

if(nargin<5||isempty(c))
    c=Constants.EGM2008GM;
end

if(nargin<4||isempty(a))
    a=Constants.EGM2008SemiMajorAxis;
end

M=C.numClusters()-1;

if(M<3)
    error('The coefficients must be provided to at least degree 3. To use a lower degree, one can insert zero coefficients.');
end

numPoints=size(point,2);
%If we are evaluating terrain heights.
switch(size(point,1))
    case 2
        a=1;
        c=1;
        point=[ones(1,numPoints);point(1,:);point(2,:)];
    case 3
    otherwise
        error('Invalid point length');
end

%If the coefficients are Schmidt-quasi-normalized, then convert them to
%fully normalized coefficients.
if(fullyNormalized==false)
    %Duplicate the input coefficients so that when they are modified, the
    %orignal values are not changed.
    C=C.duplicate();
    S=S.duplicate();

    for n=0:M
        %Deal with all of the other m values.
        k=1/sqrt(1+2*n);
        for m=0:n
            C(n+1,m+1)=k*C(n+1,m+1);
            S(n+1,m+1)=k*S(n+1,m+1);
        end
    end
end

%%If a compiled C++ implementation exists, then just use that.
if(exist('spherHarmonicEvalCPPInt','file'))
    if(nargout==2)
        [V,gradV]=spherHarmonicEvalCPPInt(C.clusterEls,S.clusterEls,point,a,c,scalFactor);
    else
        V=spherHarmonicEvalCPPInt(C.clusterEls,S.clusterEls,point,a,c,scalFactor);
    end
    return
end

%Preallocate space used by the modified forward row algorithm when
%evaluating over multiple values with the same range and latitude but
%different longitudes.
XC=zeros(M+1,1);
XS=zeros(M+1,1);
if(nargout>1)
    XCdr=zeros(M+1,1);
    XSdr=zeros(M+1,1);
    XCdTheta=zeros(M+1,1);
    XSdTheta=zeros(M+1,1);
end

V=zeros(numPoints,1);
gradV=zeros(3,numPoints);

rPrev=Inf;
thetaPrev=Inf;
for curPoint=1:numPoints
    r=point(1,curPoint);
    lambda=point(2,curPoint);
    thetaCur=point(3,curPoint);
    
    rChanged=r~=rPrev;
    thetaChanged=thetaCur~=thetaPrev;
    rPrev=r;
    thetaPrev=thetaCur;

    if(rChanged)
        %This stores all of the powers of a/r needed for the sum,
        %regardless of which algorithm is used.
        nCoeff=zeros(M+1,1);
        nCoeff(1)=1;
        for n=1:M
            nCoeff(n+1)=nCoeff(n)*(a/r);
        end
    end

    if(abs(thetaCur)<88*pi/180||nargout<2)
        %At latitudes that are not near the poles, the algorithm of Holmes and
        %Featherstone is used. It can not be used for the gradient near the
        %poles, because of the singularity of the spherical coordinate system.
        
        %Compute the sine and cosine terms.
        [SinVec,CosVec]=calcSinCosTerms(lambda,M);
        
        theta=pi/2-thetaCur;
        u=sin(theta);
        if(thetaChanged)
            %The spherical coordinate system used in ellips2Sphere uses azimuth and
            %elevation (latitude). However the formulae for spherical harmonic
            %synthesis in the Holmes and Featherstone paper use, pi/2-elevation
            %(colatitude). Thus, the point must be transformed.
            
            [PBarUVals,dPBarUValsdTheta]=NALegendreCosRat(theta,M,scalFactor);
        end
        
        %Evaluate Equation 7 from the Holmes and Featherstone paper.
        if(rChanged||thetaChanged)
            XC=zeros(M+1,1);
            XS=zeros(M+1,1);
            
            for m=0:M
                for n=m:M
                    XC(m+1)=XC(m+1)+nCoeff(n+1)*C(n+1,m+1)*PBarUVals(n+1,m+1);
                    XS(m+1)=XS(m+1)+nCoeff(n+1)*S(n+1,m+1)*PBarUVals(n+1,m+1);
                end
            end
        end
        
        %Use Horner's method to compute V.
        V(curPoint)=0;
        for m=M:-1:0
            OmegaRat=XC(m+1)*CosVec(m+1)+XS(m+1)*SinVec(m+1);
            V(curPoint)=V(curPoint)*u+OmegaRat;
        end
        
        %Multiply by the constant in front of the sum and get rid of the scale
        %factor.
        V(curPoint)=(c/r)*V(curPoint)/scalFactor;
        
        %If the gradient is desired.
        if(nargout>1)
            dVdr=0;
            dVdLambda=0;
            dVdTheta=0;
            
            if(rChanged||thetaChanged)
                XCdr=zeros(M+1,1);
                XSdr=zeros(M+1,1);
                XCdTheta=zeros(M+1,1);
                XSdTheta=zeros(M+1,1);
                
                %Evaluate Equation 7 from the Holmes and Featherstone paper.
                for m=0:M
                    for n=m:M
                        CScal=nCoeff(n+1)*C(n+1,m+1);
                        SScal=nCoeff(n+1)*S(n+1,m+1);

                        XCdr(m+1)=XCdr(m+1)+(n+1)*CScal*PBarUVals(n+1,m+1);
                        XSdr(m+1)=XSdr(m+1)+(n+1)*SScal*PBarUVals(n+1,m+1);

                        XCdTheta(m+1)=XCdTheta(m+1)+CScal*dPBarUValsdTheta(n+1,m+1);
                        XSdTheta(m+1)=XSdTheta(m+1)+SScal*dPBarUValsdTheta(n+1,m+1);
                    end
                end
            end
            
            for m=M:-1:0
                OmegaRat=XCdr(m+1)*CosVec(m+1)+XSdr(m+1)*SinVec(m+1);
                dVdr=dVdr*u+OmegaRat;

                OmegaRat=m*(-XC(m+1)*SinVec(m+1)+XS(m+1)*CosVec(m+1));
                dVdLambda=dVdLambda*u+OmegaRat;

                OmegaRat=XCdTheta(m+1)*CosVec(m+1)+XSdTheta(m+1)*SinVec(m+1);
                dVdTheta=dVdTheta*u+OmegaRat;
            end
            
            dVdr=-(c/r^2)*dVdr/scalFactor;
            dVdLambda=(c/r)*dVdLambda/scalFactor;
            %The minus sign is because the input coordinate was with respect to
            %latitude, not the co-latitude that the NALegendreCosRat function uses.
            dVdTheta=-(c/r)*dVdTheta/scalFactor;

            gradV(:,curPoint) = calcSpherJacob(point(:,curPoint))'*[dVdr;dVdLambda;dVdTheta];
        end
    else
        %At latitudes that are near the poles, the non-singular algorithm of
        %Pines using the fully normalized Helmholtz equations from Fantino and
        %Casotto is used. The algorithm has been slightly modified so that the
        %c/r term is out front and the fully normalized Helmholtz polynomials
        %can be scaled. Also, lumped coefficients are not used. The Pines
        %algorithm is generally slower than the algorithm of Holmes and
        %Featherstone and it is suffers a loss of precision near the equator.
        %Thus, the Pines algorithm is only used near the poles where the other
        %algorithm has issues with a singularity.
        
        CartPoint=spher2Cart(point(:,curPoint));

        x=CartPoint(1);
        y=CartPoint(2);
        z=CartPoint(3);

        %Get the direction cosines used by Pines' algorithm.
        s=x/r;
        t=y/r;
        u=z/r;

        %Compute the fully normalized Helmholtz polynomials.
        if(thetaChanged)
            [HBar,dHBardu]=normHelmholtz(u,M,scalFactor);
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

        %Perform the sum for the potential from Equation 44 in the Fantino and
        %Casotto paper.
        V(curPoint)=0;
        for n=0:M
            innerTerm=0;
            for m=0:n
                innerTerm=innerTerm+(C(n+1,m+1)*rm(m+1)+S(n+1,m+1)*im(m+1))*HBar(n+1,m+1);
            end
            V(curPoint)=V(curPoint)+nCoeff(n+1)*innerTerm;
        end

        V(curPoint)=(c/r)*V(curPoint)/scalFactor;

        if(nargout>1)
            %Now, compute the gradient.
            a1=0;
            a2=0;
            a3=0;
            a4=0;
            %The equations in these loops are from Table 10.
            for n=0:M
                a1Loop=0;
                a2Loop=0;

                %The m=0 case only applies to a3 and a4.
                m=0;
                HVal=HBar(n+1,m+1);
                dHVal=dHBardu(n+1,m+1);
                CProdMN=C(n+1,m+1)*rm(m+1)+S(n+1,m+1)*im(m+1);

                a3Loop=CProdMN*dHVal;
                Lmn=(n+m+1)*HVal+u*dHVal;
                a4Loop=-CProdMN*Lmn;

                for m=1:n
                    HVal=HBar(n+1,m+1);
                    dHVal=dHBardu(n+1,m+1);

                    a1Loop=a1Loop+m*(C(n+1,m+1)*rm(m-1+1)+S(n+1,m+1)*im(m-1+1))*HVal;
                    a2Loop=a2Loop+m*(S(n+1,m+1)*rm(m-1+1)-C(n+1,m+1)*im(m-1+1))*HVal;

                    CProdMN=C(n+1,m+1)*rm(m+1)+S(n+1,m+1)*im(m+1);

                    a3Loop=a3Loop+CProdMN*dHVal;
                    Lmn=(n+m+1)*HVal+u*dHVal;
                    a4Loop=a4Loop-CProdMN*Lmn;
                end

                a1=a1+nCoeff(n+1)*a1Loop;
                a2=a2+nCoeff(n+1)*a2Loop;
                a3=a3+nCoeff(n+1)*a3Loop;
                a4=a4+nCoeff(n+1)*a4Loop;
            end

            %These are equation 70. However, an additional 1/r term has been added,
            %which the original paper omitted when going from Equation 68 to 70.
            gradV(1,curPoint)=(c/r^2)*(a1+s*a4)/scalFactor;
            gradV(2,curPoint)=(c/r^2)*(a2+t*a4)/scalFactor;
            gradV(3,curPoint)=(c/r^2)*(a3+u*a4)/scalFactor;
        end
    end

end
end

function [SinVec,CosVec]=calcSinCosTerms(lambda,M)
    %Compute sin(m*lambda) and cos(m*lambda) for m=0 to m=M.
    SinVec=zeros(M+1,1);
    CosVec=zeros(M+1,1);
    %Explicitely set the first two terms.
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
