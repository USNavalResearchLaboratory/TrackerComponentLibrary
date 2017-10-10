function xpBarypBar=meanRotPoleLoc(TT1,TT2,deltaTTUT1,clockLoc)
%%MEANROTPOLELOC Compute the location of the mean rotation pole according
%                to the model of the IERS 2010 Conventions.
%
%INPUTS: TT1,TT2 Two parts of a Julian date given in terrestrial time (TT). 
%                The units of the date are days. The full date is the sum
%                of both terms. The date is broken into two parts to
%                provide more bits of precision. It does not matter how the
%                date is split.
%     deltaTTUT1 An optional parameter specifying the offset between TT and
%                UT1 in seconds. If this parameter is omitted or an empty
%                matrix is passed, then the value of the function getEOP
%                will be used.
%       clockLoc An optional 3X1 vector specifying the location of the
%                clock in WGS-84 ECEF Cartesian [x;y;z] coordinates with
%                units of meters. Due to relativistic effects, clocks that
%                are synchronized with respect to TT are not synchronized
%                with respect to TDB. If this parameter is omitted, then a
%                clock at the center of the Earth is used.
%
%OUTPUTS: xpBarypBar xpBarypBar=[xpBar;ypBar] are the location of the mean
%                    rotation pole as defined in Section 7.1.4 of the IERS
%                    2010 conventions, given in radians.
%
%The 2010 International Earth Rotation System Service (IERS) conventions
%are [1]. This models the secular change in the polar motion variables. The
%mean rotation pole comes up when dealing with pole tides and when
%computing the potential due to the rotation of the Earth.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    deltaTTUT1=[];
end

if(nargin<4)
    clockLoc=[];
end

%These are the coefficients in Table 7.7 for computing the IERS 2010 mean
%pole model.
CoeffMat=[0, 55.974,    346.346,   23.513,  358.891;
          1, 1.8243,     1.7896,   7.6141,  -0.6287;
          2, 0.18413,   -0.10729,  0.0,      0.0;
          3, 0.007024,  -0.000908, 0.0,      0.0];

T0=2000.0;

%Get the Besselian epoch for the TT Julian date.
T=TT2BesselEpoch(TT1,TT2,deltaTTUT1,clockLoc);
tDiff=T-T0;

numRows=size(CoeffMat,1);
xpBar=0;
ypBar=0;
for curRow=1:numRows
    %The degree
    i=CoeffMat(curRow,1);
    
    if(T<=2010)
        xpBari=CoeffMat(curRow,2);
        ypBari=CoeffMat(curRow,3);
    else
        xpBari=CoeffMat(curRow,4);
        ypBari=CoeffMat(curRow,5);
    end
    
    xpBar=xpBar+tDiff^i*xpBari;
    ypBar=ypBar+tDiff^i*ypBari;
end

%The constant to convert milliarcseconds to radians. The multiplications go
%mas->as->arcminutes->deg->rad
mas2Rad=(1/1000)*(1/60)*(1/60)*(pi/180);

xpBarypBar=[xpBar;ypBar]*mas2Rad;

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
