function vec=DelaunayVar(TDB1,TDB2)
%%DELAUNAYVAR Get the Delaunay variables that are needed for tide models
%             involving Doodson numbers using the equations in the IERS
%             conventions.
%
%INPUTS: T The time as measured in Julian centuries Barycentric dynamical
%          time (TDB).
%
%OUTPUTS: vec A vector of the Delaunay variables in the order with all
%             units in RADIANS.
%             vec(1) Mean anomaly of the moon: l
%             vec(2) Mean anomaly of the sun: l'
%             vec(3) F=L-Omega, The mean longitude of the moon minus the
%                    mean longitude of the ascending node of the moon.
%             vec(4) Mean Elongation of the moon from the sun: D
%             vec(5) Mean Longitude of the Ascending Node of the Moon:
%                    Omega
%
%This function implements Equation 5.43 in Section 5.7.2 of [1]. The units
%in those equations are degrees (and arcseconds). However, in Equation 6.8
%in Section 6.2, it looks like the angles should be in radians, since one
%is using trigonometric functions. Thus, the return value here is put into
%radians to be convenient.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The time as measured in Julian centuries Barycentric dynamical time
    %(TDB). The precision of the model is low enough that it shouldn't
    %matter that the two parts are added together.
    T=(TDB1+TDB2)/36525;

    vec=zeros(5,1);

    %The conversion from arcseconds to degrees
    sec2Deg=1/3600;

    %Mean anomaly of the moon: l
    vec(1)=134.96340251+sec2Deg*(1717915923.2178*T+31.8792*T^2+0.051635*T^3-0.00024470*T^4);
    %Mean anomaly of the sun: l'
    vec(2)=357.52910918+sec2Deg*(129596581.0481*T-0.5532*T^2+0.000136*T^3-0.00001149*T^4);
    %L-Omega, The mean longitude of the moon minus the mean longitude of the
    %ascending node of the moon.
    vec(3)=93.27209062+sec2Deg*(1739527262.8478*T-12.7512*T^2-0.001037*T^3+0.00000417*T^4);
    %Mean Elongation of the moon from the sun: D
    vec(4)=297.85019547+sec2Deg*(1602961601.2090*T-6.3706*T^2+0.006593*T^3-0.00003169*T^4);
    %Mean Longitude of the Ascending Node of the Moon: Omega
    vec(5)=125.04455501+sec2Deg*(-6962890.5431*T+7.4722*T^2+0.007702*T^3-0.00005939*T^4);

    %Convert from degrees to radians.
    vec=mod(vec*(pi/180),2*pi);
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
