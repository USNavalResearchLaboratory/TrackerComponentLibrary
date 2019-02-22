function [deltaC,deltaS]=getGravCoeffOffset4Drift(TT1,TT2,model)
%%GETGRAVCOEFFOFFSET4DRIFT Get the offset of the coefficients for a
%                   particular spherical harmonic model of the Earth's
%                   gravitation due to drift over time. The drifts over
%                   time occur for many reasons, such as because of the
%                   rebound of the tectonic plates after the glaciers
%                   receded at the end of the last ice age.
%
%INPUTS:  TT1,TT2 Two parts of a Julian date given in terrestrial time
%                 (TT).
%                    The units of the date are days. The full date is
%                   the sum of both terms. The date is broken into two
%                   parts to provide more bits of precision, though the
%                   drift models are so low-firdelity that it does not
%                   matter in this instance.
%             model An optional parameter indicating the gravitational
%                   model for which the drift terms should be computed. The
%                   default if omitted is 0. Possible values are
%                   0) The EGM2008 model is used, so the drift terms as
%                      given by the IERS 2010 conventions are returned.
%                   1) The EGM96 gravitational model is used, so the drift
%                      terms given in the readme filte for that model are
%                      returned.
%
%OUTPUTS: deltaC An array holding the offsets to the fully normalized
%               coefficient terms that are multiplied by cosines in the
%               spherical harmonic expansion of the gravitational
%               potential. If passed to a CountingClusterSet class, then
%               C(n+1,m+1) is the coefficient of degree n and order m. The
%               offsets only go up to the degree and order of the terms in
%               the gravitational model that would have to be modified.
%        deltaS An array holding the coefficient terms that are multiplied
%               by sines in the spherical harmonic expansion. The format of
%               S is the same as that of C.
%
%Both of the drift models are given in terms of a rate per years after an
%epoch date. The temporal coordinate system is not  very well defined, but
%one would assume that it should be in Julian years TT after the epoch
%date. A Julian year has exactly 365.25 days of exactly 86400 seconds.
%
%The drift terms for the low-order coefficients in the EGM2008 model are
%documented in Table 6.2 on page 80 of [1]; the adjustments for polar
%motion are from page 81.
%
%The drift terms for the EGM96 model were included in the readme file for
%the model when it was distributed by NASA. Given newer data it is likely
%that these drift terms are incorrect/ out of date since newer models do
%not put drifts on the S terms.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    model=0;
end

M=4;%Coefficients up to order 4 can be modified.
totalNumCoeffs=(M+1)*(M+2)/2;
deltaC=zeros(totalNumCoeffs,1);
deltaS=zeros(totalNumCoeffs,1);

if(model==0)%If using the EGM2008 model.
    %The epoch Date is J2000.0
    [TTEpoch1,TTEpoch2]=Cal2TT(2000,1,1,0,0,0);
    if(TT1>TT2)
        deltaT=(TT1-TTEpoch1)+(TT2-TTEpoch2);
    else
        deltaT=(TT2-TTEpoch1)+(TT1-TTEpoch2);
    end
    
    %Convert the time difference from Julian days to Julian years.
    deltaT=deltaT/365.25;
    
    %J2 drift rate, adjusts C(2+1,0+1)
    deltaC(4)=Constants.EGM2008C20BarDot*deltaT;
    %J3 drift rate, adjusts C(3+1,0+1)
    deltaC(7)=Constants.EGM2008C30BarDot*deltaT;
    %J4 drift rate, adjusts C(4+1,0+1)
    deltaC(11)=Constants.EGM2008C40BarDot*deltaT;
else
    %The epoch Date is 1986.0
    [TTEpoch1,TTEpoch2]=Cal2TT(1986,1,1,0,0,0);
    if(TT1>TT2)
        deltaT=(TT1-TTEpoch1)+(TT2-TTEpoch2);
    else
        deltaT=(TT2-TTEpoch1)+(TT1-TTEpoch2);
    end
    
    %Convert the time difference from Julian days to Julian years.
    deltaT=deltaT/365.25;
    
    %J2 drift rate, adjusts C(2+1,0+1)
    deltaC(4)=Constants.EGM96C20BarDot*deltaT;
    %C21 drift rate, adjusts C(2+1,1+1)
    deltaC(5)=Constants.EGM96C21BarDot*deltaT;
    %S21 drift rate, adjusts S(2+1,1+1)
    deltaS(5)=Constants.EGM96S21BarDot*deltaT;
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
