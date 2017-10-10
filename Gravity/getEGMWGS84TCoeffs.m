function [C,S,a,c]=getEGMWGS84TCoeffs(M,useNGAApprox,modelType)
%%GETEGMWGS84TCOEFFS Obtain spherical harmonic coefficients for computing
%            the disturbing potential T at a point when using the Earth
%            Gravitational Model 2008 (EGM2008) or Earth Gravitational
%            Model 1996 (EGM96) with the WGS-84 reference ellipsoid. These
%            coefficients can be used, with other information, to determine
%            geoid heights above the WGS-84 ellipsoid.
%
%INPUTS: M The integer maximum order of the spherical harmonic coefficients
%          obtained. This is a value between 2 and 2190 for the EGM2008
%          model and betwen 2 and 360 for the EGM96 model. If this
%          parameter is omitted or an empty matrix is passed, the default
%          value being the total number of coefficients in the model is
%          used.
% useNGAApprox If true, use the same approximation as the National
%          Geospatial Intelligence Agency (NGA) in their public code for
%          geoid undulations, which omits certain terms. If this parameter
%          is omitted, then the default value of true is used.
% modelType An optional parameter specifying coefficient model to use.
%          Possible values are
%          0 (The default if omitted) Use the EGM2008 model.
%          1 Use the EGM96 model.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the harmonic expansion. This can be given to a
%           CountingClusterSet class such that C(n+1,m+1) is the
%           coefficient of degree n and order m. When a maximum degree of
%           M is used, all C have values for all n from 0 to M and for all
%           m from 0 to n for each n. The coefficients are unitless.
%         S An array holding the coefficient terms that are multiplied by
%           sines in the harmonic expansion. The format of S is the same as
%           that of C.
%         a The numerator in the (a/r)^n term in the spherical harmonic
%           sum, having units of meters.
%         c The constant value by which the spherical harmonic series is
%           multiplied, having units of m^3/s^2.
%
%The coefficients are meant for use in the function
%getEGMGeoidHeight(); However, they can be used in general for
%computing the disturbing potential when using the function
%spherHarmonicEval. Note that corrections for the values of the universal
%gravitational constant times the mass of the Earth and the semi-major axis
%of the Earth are included when useNGAApprox=false. When useNGAApprox=true,
%the inclusion occurs when using the EGM2008 model, but not when using the
%EGM96 model as the reference code 
%
%The disturbing potential is described in Chapter 2.12 of [1].
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        modelType=0;%EGM2008 Model
    end
    
    if(nargin<2)
        useNGAApprox=true;
    end
    
    if(nargin<1||isempty(M))
        if(modelType==0)
            M=2190;%EGM2008 Model
        else
            M=360;%EGM96 Model
        end
    end

    GMEllipse=Constants.WGS84GMWithAtmosphere;
    aEllipse=Constants.WGS84SemiMajorAxis;
    omegaEllipse=Constants.WGS84EarthRotationRate;
    fEllipse=Constants.WGS84Flattening;
    GMPotential=Constants.EGM2008GM;
    aPotential=Constants.EGM2008SemiMajorAxis;

    %Load the coefficients for the EGM2008 or EGM96 geopotential model.
    [C,S,a,c]=getEGMGravCoeffs(M,true,modelType);
    
%Create the coefficients for the disturbing potential as described in
%"There is no such thing as 'The' EGM96 geoid: Subtle points on the use of
%a global geopotential model" when considering the fact that aPotential
%and GMPotential in the models might be different. However, to get the
%highest agreement with the results of the NGA's EGM96 code, the
%corrections for GM and a are omitted. Also, the NGA only uses the first
%10 terms of ellipsGravCoeffs for the differences its implementation.

    if(useNGAApprox&&modelType==1)
        GMRat=1;
        aRat=1;
        M=10;

        %The S terms in the ellipsoidal model should all be zero.
        CEllips=ellipsGravCoeffs(M,true,omegaEllipse,aEllipse,fEllipse,GMEllipse);
    else
        GMRat=GMEllipse/GMPotential;
        aRat=aEllipse/aPotential;
        %The S terms in the ellipsoidal model should all be zero.
        CEllips=ellipsGravCoeffs(M,true,omegaEllipse,aEllipse,fEllipse,GMEllipse);
    end
    
    C=CountingClusterSet(C);
    CEllips=CountingClusterSet(CEllips);
    aRecur=1;
    for n=0:M
        C(n+1,:)=C(n+1,:)-GMRat*aRecur*CEllips(n+1,:);
        aRecur=aRecur*aRat;
    end
    clear CEllips;
    C=C.clusterEls;
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
