function pointsHarmon=Cart2EllipsHarmon(cartPoints,E)
%%CART2ELLIPSHARMON Convert from Cartesian coordinates to ellipsoidal
%                   harmonic coordinates. Note that ellipsoidal harmonic
%                   coordinates are not the same as ellipsoidal
%                   coordinates.
%
%INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
%                   that are to be transformed into ellipsoidal
%                   coordinates. Each columns of cartPoints is of the
%                   format [x;y;z].
%                 E The linear eccentricity defining the ellipsoidal
%                   harmonic coordinate system. If this parameter is
%                   omitted, then the linear eccentricity of the WGS84
%                   reference ellipsoid is used.
%
%OUTPUTS: pointsHarmon A matrix of the converted points. Each column of the
%                   matrix has the format
%                   [reduced latitude;longitude;semiminor axis], with
%                   reduced latitude and longitude given in radians.
%
%The ellipsoidal harmonic coordinate system is described in Chapter 1.15 of
%[1]. Note that some folks use the complement of the reduced latitude in
%place of the reduced latitude. The complement is pi/2-the reduced
%latitude.
%
%The conversion should work for all points that are not at the origin. For
%points that are particularly close to the origin, (generally, points deep
%within the Earth) the term
%t=x.^2+y.^2+z.^2-E^2;
%will be negative. If t>0, it can be easily verified that
%u=sqrt((t/2)*(1+sqrt(1+4*E^2*z^2/t^2)));
%On the other hand, for t<0, it can be verified (using, for example,
%Mathematica) that the term
%sqrt((abs(t)/2).*(1+sqrt(1+4*E^2*z.^2./t.^2)));
%is equal to E*abs(cos(phi)). Thus, a different conversion is used
%depending on the value of t. The conversion follows from solving the
%equations in the aforementioned book by Hofmann-Wellenhof and Moritz.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(E))
        a=Constants.WGS84SemiMajorAxis;
        f=Constants.WGS84Flattening;
        b=a*(1 - f);
        E=sqrt(a^2-b^2);
    end

    numPoints=size(cartPoints,2);

    x=cartPoints(1,:);
    y=cartPoints(2,:);
    z=cartPoints(3,:);

    t=x.^2+y.^2+z.^2-E^2;
    uPossible=sqrt((abs(t)/2).*(1+sqrt(1+4*E^2*z.^2./t.^2)));
    
    %This deals with the case where t is very close to zero so division by
    %t can be problematic.
    if(~isfinite(uPossible))
        uPossible=0;
    end
    
    lambda=atan2(y,x);

    u=zeros(1,numPoints);
    phi=zeros(1,numPoints);
    for curPoint=1:numPoints
        if(t(curPoint)>0)
            u(curPoint)=uPossible(curPoint);
            rat=z(curPoint)/u(curPoint);
            %Deal with precision issues
            if(abs(rat)>1.0)
                rat=sign(z);
            end
            
            phi(curPoint)=acos(rat);
        else
            phi(curPoint)=acos(sign(z(curPoint))*uPossible(curPoint)/E);
            u(curPoint)=z(curPoint)/cos(phi(curPoint));
        end
    end

    beta=pi/2-phi;
    pointsHarmon=[beta;lambda;u];
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
