function plhPoints=perspectiveProj2Ellipse(xyhPoints,plhCamera,a,f)
%%PERSPECTIVEPROJ2ELLIPSE Convert points from a perspective projection with
%        respect to a particular camera location to ellipsoidal
%        coordinates. The camera is facing straight down at the Earth.
%
%INPUTS:  xyh The 3XN set of [x;y;h] points of the perspective projection.
%     plhCamera The 3X1 [latitude;longitudeheight] reference point for the
%               camera. h must be finite.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84Flattening is used.
%
%OUTPUTS: plhPoints A 3XN set of [latitude;longitude;height] of converted
%               points to with latitude and longitude given in radians.
%
%This implements the general perspective projection for an ellipsoid that
%is given in Chapter 23 of [1].
%
%EXAMPLE:
%The example on page 323 of Appendix A, (without h0) with a second point
%added to convert to show that the conversions work with more than one
%point. Here, we convert the points and then convert back to show that the
%conversions are consistent with each other. The relative error is on the
%order of finite precision limitations, as one might expect.
% a=6378206.4;
% e2=0.00676866;
% f=1-sqrt(1-e2);
% plhCamera=[deg2rad([39;-77]);500e3+200];
% plhPoint=[[deg2rad([41;-74]);100],[deg2rad([42;-75]);200]];
% xyh=ellips2PerspectiveProj(plhPoint,plhCamera,a,f);
% plhBack=perspectiveProj2Ellipse(xyh,plhCamera,a,f);
% RelErr=max(abs((plhBack(:)-plhPoint(:))./plhPoint(:)))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%It seems to always converge in under 10 iterations.
maxIter=100;

e2=f*(2-f);%Squared eccentricity.

if(~isfinite(plhCamera(3)))
    P=Inf;
else
    lCamera=ellips2Cart(plhCamera,a,f);
    P=norm(lCamera)/a;
end

%Camera coordinates:
phi1=plhCamera(1);
lambda0=plhCamera(2);
H=plhCamera(3);

x=xyhPoints(1,:);
y=xyhPoints(2,:);
h=xyhPoints(3,:);
numPts=size(x,2);

sinPhi1=sin(phi1);
cosPhi1=cos(phi1);

%Equation 8-23
N1=a./sqrt(1-e2.*sinPhi1.^2);
%Equation 23-17
phig=phi1-asin(N1.*e2.*sinPhi1.*cosPhi1./(P.*a));
%Equation 23-22
B=P.*cos(phi1-phig);
%Equation 23-23
D=P.*sin(phi1-phig);
%Equation 23-24
L=1-e2.*cosPhi1.^2;
%Equation 23-25
G=1-e2.*sinPhi1.^2;
J=2.*e2.*sinPhi1.*cosPhi1;
u=-2.*B.*L.*H-2.*D.*G.*y+B.*J.*y+D.*H.*J;
v=L.*H.^2+G.*y.^2-H.*J.*y+(1-e2).*x.^2;

cosPhig=cos(phig);
plhPoints=zeros(3,numPts);
for curPt=1:numPts
    E=1;%Initial estimate.

    EPrev=Inf;
    phiPrev=Inf;
    for curIter=1:maxIter
        t=P^2*(1-e2*cosPhig^2)-E*(1-e2);
        Kp=(-u(curPt)+sqrt(u(curPt)^2-4*t*v(curPt)))/(2*t);
        X=a*((B-H/Kp)*cosPhi1-(y(curPt)/Kp-D)*sinPhi1);
        Y=a*x(curPt)/Kp;
        S=(y(curPt)/Kp-D)*cosPhi1+(B-H/Kp)*sinPhi1;
        lambda=lambda0+atan2(Y,X);
        if(curIter==1)
            phi=asin(S);
            if(h(curPt)==0)
                break;
            end
        end
        sinPhi=sin(phi);
        phi=asin(S/((1-e2)/sqrt(1-e2*sinPhi.^2)+h(curPt)/a));
        E=(1/sqrt(1-e2*sinPhi^2)+h(curPt)/a)^2-e2*sinPhi^2*(1/(1-e2*sinPhi^2)-h(curPt)^2/(a^2-a^2*e2));

        if(abs(phi-phiPrev)==0&&abs(E-EPrev)==0)
            %If it converged.
            break;
        else
            EPrev=E;
            phiPrev=phi;
        end
    end
    plhPoints(:,curPt)=[phi;lambda;h(curPt)];
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
