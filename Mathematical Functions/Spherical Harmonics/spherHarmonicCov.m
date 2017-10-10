function [sigma2,Sigma,didOverflow]=spherHarmonicCov(CStdDev,SStdDev,point,a,c,scalFactor)
%%SPHERHARMONICCOV Evaluate the variance of a potential or the covariance
%                  matrix of a gradient that one might compute using a
%                  spherical harmonic coefficient model with
%                  spherHarmonicEval. Here, one provides the standard
%                  deviations associated with the coefficients in the
%                  model, such as those obtained using the function
%                  getEGMGravCoeffs for the EGM2008 gravity model.
%
%INPUTS: CStdDev A length (M+2)*(M+1)/2 array holding the standard
%               deviations of the coefficient terms that are multiplied by
%               cosines in the harmonic expansion. The coefficients must be
%               fully normalized using the type of full normalization that
%               is used in the EGM2008 model. Their normalization type can
%               be changed using the changeSpherHarmonicNorm function. If
%               given to a CountingClusterSet class, C(n+1,m+1) is the
%               coefficient of degree n and order m. When a maximum degree
%               of M is used, all C must have values for all n from 0 to M
%               and for all m from 0 to n for each n. If coefficients are
%               not present for certain degrees, then insert a 0. It is
%               assumed that M>=3. For magnetic field models, the units of
%               C should generally be Tesla. For gravitational models, C is
%               generally unitless with the units being determined entirely
%               by the parameters a and c.
%       SStdDev A length (M+2)*(M+1)/2 array holding the standard
%               deviations of the coefficient terms that are multiplied by
%               sines in the harmonic expansion. The requirements on
%               SStdDev are the same as those on CStdDev.
%         point The 3XN set of N points at which the potential and/or
%               gradient should be evaluated given in SPHERICAL, ECEF 
%               coordinates consisting of [r;azimuth;elevation]; When
%               evaluating points on a grid, the algorithm will be fastest
%               if the points are provided presorted by range and then by
%               azimuth. This reduces the amount of recomputation of
%               certain values. Alternatively, if C and S are for
%               evaluating terrain heights, then points are 2XN having the
%               format [azimuth;elevation] and it is best if the points are
%               sorted by azimuth.
%             a The numerator in the (a/r)^n term in the spherical harmonic
%               sum. Normally, this is some type of a reference radius. For
%               example, when using most gravitational models, a is the
%               semi-major axis of the reference ellipsoid. If this
%               parameter is omitted, it is assumed that one is using the
%               spherical harmonics with something like the National
%               Geospatial Intelligence Agency's (NGA's) EGM96 or EGM2008
%               models, in which case a=Constants.EGM2008SemiMajorAxis is
%               used unless point is 2D, in which case c=1 is used.
%             c The constant value by which the spherical harmonic series
%               is multiplied. For example, for gravitational potentials,
%               the value is usually GM where G is the universal
%               gravitational constant and M is the mass of the Earth. When
%               using the International Geomagnetic Reference Field (IGRF),
%               the constant is a^2, where a is the same as the numerator
%               in the a/r term. If this parameter is omitted, it is
%               assumed that one is using the spherical harmonics with
%               something like the NGA's EGM96 or EGM2008 models, in which
%               case c=Constants.EGM2008GM is used unless point is 2D, in
%               which case c=1 is used.
%    scalFactor An optional scale factor used in computing the normalized
%               associated Legendre polynomials if fullyNormalized=true.
%               Generally, the default value (is if the scalFactor
%               parameter is omitted) of 2^(-470) is sufficient. When very
%               high-order models are used, this scale factor prevents
%               overflows. However overflows (and a loss of precision) are
%               unavoidable when using the full EGM2008 model. These
%               effects are worse near the poles.
%
%OUTPUTS: sigma2 The NX1 vector of variances (squared standard deviations)
%                of the potential estimate at the given points.
%          Sigma The covariance matrix of the gradient of the potential at
%                the given points.
%    didOverFlow This indicates whether during the computation of sigma2 or
%                Sigma there were any overflow errors leading to term
%                being dropped. Note that this does not indicate potential
%                losses of precision due to underflow errors making terms
%                zero, which becomes more common the smaller scalFactor is.
%                Also, if scaling is extremely bad, it is possible for
%                NaNs or Inf terms to still be returned in sigma2 and
%                Sigma.
%
%The algorithm used here is described in [1].
%
%Since Matlab uses double precision arithmetic, when using high degree and
%order models, such as the full 2190 degree EGM2008 model, the precision of
%Sigma will be reduced as higher-order terms can experience overflow
%problems and are thus discarded to avoid NaNs from occurring. Lower-order
%models will not suffer from the same problem.
%
%EXAMPLE:
%Here, we evaluate the function at a point on the reference ellipsoid and
%then at a point at the pole using the full degree 2190 EGM2008 model. Note
%that the Matlab implementation will bee too slow for this high a degree
%and the C++ implementation should be compiled.
% latLongAlt=[-30*(pi/180);45*(pi/180);0];
% spherLoc=Cart2Sphere(ellips2Cart(latLongAlt));
% 
% [C,S,a,c,CStdDev,SStdDev]=getEGMGravCoeffs(2190,false);
% [sigma21,Sigma1,didOverflow1]=spherHarmonicCov(CStdDev,SStdDev,spherLoc,a,c)
% spherLoc(end)=pi/2;%90 degree elevation --the pole.
% [sigma22,Sigma2,didOverflow2]=spherHarmonicCov(CStdDev,SStdDev,spherLoc,a,c)
%One will see that in both instances, finite values are returned, but for
%the point at the pole, overflows are flagged indicating the loss of
%prevision due to some terms being dropped due to overflows. A model with
%fewer spherical harmonic coefficients would not be as susceptible to such
%overflow problems.
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(scalFactor))
    scalFactor=2^(-470);
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

M=(1/2)*(sqrt(1+8*length(CStdDev))-1)-1;

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

CStdDev=CountingClusterSet(CStdDev);
SStdDev=CountingClusterSet(SStdDev);

sigma2=zeros(numPoints,1);
Sigma=zeros(3,3,numPoints);

didOverflow=false;
rPrev=Inf;
thetaPrev=Inf;
for curPoint=1:numPoints
    r=point(1,curPoint);
    thetaCur=point(3,curPoint);
    
    rChanged=r~=rPrev;
    rPrev=r;

    if(rChanged)
        nCoeff=zeros(M+1,1);
        nCoeff(1)=1;
        for n=1:M
            nCoeff(n+1)=nCoeff(n)*(a/r);
        end
    end

    %The non-singular algorithm of Pines using the fully normalized
    %Helmholtz equations from Fantino and Casotto is used. The algorithm
    %has been slightly modified so that the c/r term is out front and the
    %fully normalized Helmholtz polynomials can be scaled. Also, lumped
    %coefficients are not used. The Pines algorithm can suffer a loss of
    %precision near the equator. However, it is simple to just square the
    %terms in the sum.
    CartPoint=spher2Cart(point(:,curPoint));

    x=CartPoint(1);
    y=CartPoint(2);
    z=CartPoint(3);

    %Get the direction cosines used by Pines' algorithm.
    s=x/r;
    t=y/r;
    u=z/r;

    %Compute the fully normalized Helmholtz polynomials.
    if(thetaPrev~=thetaCur)
        [HBar,dHBardu]=normHelmholtz(u,M,scalFactor);
        thetaPrev=thetaCur;
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
    %Casotto paper, but square all of the terms to represent a sigma.
    %All cross terms are expected to be zero, so the sum is very similar to
    %the sum for the potential in spherHarmonicEval.
    sigma2(curPoint)=0;
    for n=0:M
        innerTerm=0;
        for m=0:n
            innerTerm=innerTerm+(CStdDev(n+1,m+1)*rm(m+1)*HBar(n+1,m+1))^2+(SStdDev(n+1,m+1)*im(m+1)*HBar(n+1,m+1))^2;
        end
        if(isfinite(innerTerm))
            sigma2(curPoint)=sigma2(curPoint)+nCoeff(n+1)^2*innerTerm;
        else
            didOverflow=true;
        end
    end

    sigma2(curPoint)=(c/r)^2*sigma2(curPoint)/scalFactor^2;

    %Now, compute the cosigma matrix of the gradient, if requested.
    if(nargout>1)
        a11=0;
        a22=0;
        a33=0;
        a44=0;
        a12=0;
        a13=0;
        a14=0;
        a23=0;
        a24=0;
        a34=0;
    
        %The equations in these loops are from Table 10.
        for n=0:M
            a11Loop=0;
            a22Loop=0;
            a12Loop=0;
            a13Loop=0;
            a14Loop=0;
            a23Loop=0;
            a24Loop=0;

            %The m=0 case only applies to a3 and a4, so that means only to
            %a33, a34, and a44.
            m=0;
            HVal=HBar(n+1,m+1);
            dHVal=dHBardu(n+1,m+1);
            CProdMN=CStdDev(n+1,m+1)^2*rm(m+1)^2+SStdDev(n+1,m+1)^2*im(m+1)^2;
            Lmn=(n+m+1)*HVal+u*dHVal;%Defined in Table 14.
            
            a33Loop=CProdMN*dHVal^2;  
            a44Loop=CProdMN*Lmn^2;
            a34Loop=-CProdMN*Lmn*dHVal;

            for m=1:n
                HVal=HBar(n+1,m+1);
                dHVal=dHBardu(n+1,m+1);
                
                CProdMN=CStdDev(n+1,m+1)^2*rm(m+1)^2+SStdDev(n+1,m+1)^2*im(m+1)^2;
                Lmn=(n+m+1)*HVal+u*dHVal;
                
                %These if-statements are to deal with numerical precision
                %problems near the poles. We want to avoid 0*Inf terms due
                %to limitations in the valid range of double precision
                %numbers. Of course, the loss of the terms where overflow
                %occurs means that the covariance matrix will be
                %underestimated.
                if(isfinite(HVal))
                    a11Loop=a11Loop+m^2*(CStdDev(n+1,m+1)^2*rm(m-1+1)^2+SStdDev(n+1,m+1)^2*im(m-1+1)^2)*HVal^2;
                    a12Loop=a12Loop+m^2*rm(m-1+1)*im(m-1+1)*(SStdDev(n+1,m+1)^2-CStdDev(n+1,m+1)^2)*HVal^2;
                    a22Loop=a22Loop+m^2*(SStdDev(n+1,m+1)^2*rm(m-1+1)^2+CStdDev(n+1,m+1)^2*im(m-1+1)^2)*HVal^2;
                else
                    didOverflow=true;
                end
                
                if(isfinite(Lmn))
                    a44Loop=a44Loop+CProdMN*Lmn^2;
                    if(isfinite(HVal))
                        a14Loop=a14Loop-m*(CStdDev(n+1,m+1)^2*rm(m-1+1)*rm(m+1)+SStdDev(n+1,m+1)^2*im(m-1+1)*im(m+1))*HVal*Lmn;
                        a24Loop=a24Loop-m*(-CStdDev(n+1,m+1)^2*im(m-1+1)*rm(m+1)+SStdDev(n+1,m+1)^2*rm(m-1+1)*im(m+1))*HVal*Lmn;
                    end
                    if(isfinite(dHVal))
                        a34Loop=a34Loop-CProdMN*Lmn*dHVal;
                    end
                else
                    didOverflow=true;
                end
                
                if(isfinite(dHVal))
                    a33Loop=a33Loop+CProdMN*dHVal^2;
                    if(isfinite(HVal))
                        a13Loop=a13Loop+m*(CStdDev(n+1,m+1)^2*rm(m-1+1)*rm(m+1)+SStdDev(n+1,m+1)^2*im(m-1+1)*im(m+1))*HVal*dHVal;
                        a23Loop=a23Loop+m*(-CStdDev(n+1,m+1)^2*im(m-1+1)*rm(m+1)+SStdDev(n+1,m+1)^2*rm(m-1+1)*im(m+1))*HVal*dHVal;
                    end
                else
                    didOverflow=true;
                end
            end
            
            a11=a11+nCoeff(n+1)^2*a11Loop;
            a22=a22+nCoeff(n+1)^2*a22Loop;
            a33=a33+nCoeff(n+1)^2*a33Loop;
            a44=a44+nCoeff(n+1)^2*a44Loop;
            a12=a12+nCoeff(n+1)^2*a12Loop;
            a13=a13+nCoeff(n+1)^2*a13Loop;
            a14=a14+nCoeff(n+1)^2*a14Loop;
            a23=a23+nCoeff(n+1)^2*a23Loop;
            a24=a24+nCoeff(n+1)^2*a24Loop;
            a34=a34+nCoeff(n+1)^2*a34Loop;
        end

%These are based on squaring the terms in equation 70, removing cross
%terms. However, an additional 1/r (squared) term has been added,
%which the original paper omitted when going from Equation 68 to 70.
        
        s11=a11+2*s*a14+s^2*a44;
        s12=a12+s*a24+t*a14+s*t*a44;
        s13=a13+s*a34+u*a14+s*u*a44;
        s22=a22+2*t*a24+t^2*a44;
        s23=a23+t*a34+u*a24+t*u*a44;
        s33=a33+2*u*a34+u^2*a44;
        
        temp=c/(r^2*scalFactor);

        Sigma(:,:,curPoint)=temp*(temp*[s11,s12,s13;
                                        s12,s22,s23;
                                        s13,s23,s33]);
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
