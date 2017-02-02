function J=calcPolarJacob(point,systemType)
%%CALCPOLARJACOB  Compute the Jacobian matrix for a point in polar
%                 [range;azimuth] coordinates.  The angle can be measured
%                 either counterclockwise from the x-axis, which is
%                 standard in mathematics, or clockwise from the y axis,
%                 which is more common in navigation.
%
%INPUTS: point   A point in the format [range;azimuth], where the angle is
%                given in radians.
%      systemType   An optional parameter specifying the axes from which
%                   the azimuth angle is measured. Possible vaues are
%                   0 (The default if omitted) The azimuth angle is
%                     counterclockwise from the x axis.
%                   1 The azimuth angle is measured clockwise from the y
%                     axis.
%
%OUTPUTS: J     The 2X2 Jacobian matrix. Each row is a components of range,
%               and azimuth (in that order by row) with derivatives taken
%               with respect to [x,y] by column.
%
%May 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    systemType=0;
end

CartPoint=pol2Cart(point,systemType);
x=CartPoint(1);
y=CartPoint(2);

r=point(1);
J=zeros(2,2);

switch(systemType)
    case 0
        %Derivatives with respect to x.
        J(1,1)=x/r;
        J(2,1)=-y/(x^2+y^2);

        %Derivatives with respect to y.
        J(1,2)=y/r;
        J(2,2)=x/(x^2+y^2);
    case 1
        %Derivatives with respect to x.
        J(1,1)=x/r;
        J(2,1)=y/(x^2+y^2);
        
        %Derivatives with respect to y.
        J(1,2)=y/r;
        J(2,2)=-x/(x^2+y^2);
    otherwise
        error('Invalid system type specified.')
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
