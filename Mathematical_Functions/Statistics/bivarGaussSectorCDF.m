function M=bivarGaussSectorCDF(theta1,theta2,R)
%%BIVARGAUSSECTORCDF Find the value of the cumulative distirbution function
%               (CDF) for an angular region spanning from theta1 to theta2
%               over the zero-mean normal distribution with covariance
%               matrix R. The sector considered always starts from the
%               origin. The angular region spanned goes starts at theta 1
%               and goes in the direction of increasing angle to theta2.
%               Numerically, if crossing a -pi/pi boundary, it is okay if
%               theta2 is less than theta1. It is possible to span a region
%               more than pi radians in width.
%
%INPUTS: theta1, theta2 The starting and ending angles of the region over
%               which the CDF value is desired. The values must be between
%               in the range of -pi to pi radians.
%             R The 2X2 matrix of the distribution.
%
%OUTPUTS: M The probability within the specified region.
%
%This implements the solution from Theorem 1 of [1]. What appears to be a
%mistake in the sign of a term in Equation 11 has been corrected.
%
%EXAMPLE:
%Here, we compute the probability using this function as well as the
%probability using a Monte Carlo method. The probabilities will tend to be
%close.
% theta1=-pi/2;
% theta2=pi/3;
% rho=-0.6;
% R=[2,rho;
%    rho,6];
% S=chol(R,'lower');
% prob=bivarGaussSectorCDF(theta1,theta2,R)
% numPts=1e7;
% z=S*randn(2,numPts);
% theta=atan2(z(2,:),z(1,:));
% MCProb=sum(pointIsInPolAngSpan(theta,theta1,theta2))/numPts
%
%REFERENCES:
%[1] V. Savaux and L. Le Magoarou, "On the computation of integrals of
%    bivariate Gaussian distribution," in Proceedings of the IEEE Symposium
%    on Computers and Communications, Rennes, France, 7-10 Jul. 2020.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(theta2<theta1)
    thetaDiff=(theta2+2*pi)-theta1;
else
    thetaDiff=theta2-theta1;
end

isLongWayAround=thetaDiff>pi;

R12=R(1,2);
sigmaX2=R(1,1);
sigmaY2=R(2,2);

%Equation 8
if(sigmaX2>sigmaY2)
    phi=(1/2)*atan(2*R12/(sigmaX2-sigmaY2));
elseif(sigmaX2<sigmaY2)
    phi=pi/2+(1/2)*atan(2*R12/(sigmaX2-sigmaY2));
else
    phi=pi/4;
end

sigmaXi=sqrt(sigmaX2*cos(phi)^2+sigmaY2*sin(phi)^2+R12*sin(2*phi));
sigmaNu=sqrt(sigmaX2*sin(phi)^2+sigmaY2*cos(phi)^2-R12*sin(2*phi));

p=sigmaNu/sigmaXi;

%From the definition in Fig. 2.
theta1=wrapRange(theta1-phi,-pi,pi);
theta2=wrapRange(theta2-phi,-pi,pi);

thetam=min(abs(theta1),abs(theta2));
thetaM=max(abs(theta1),abs(theta2));

theta1InRange1=(theta1>=-pi/2&&theta1<=pi/2);
theta2InRange1=(theta2>=-pi/2&&theta2<=pi/2);

if(theta1InRange1&&theta2InRange1)
    if(sign(theta1)==sign(theta2))
        Q=uFunc(thetaM,p)-uFunc(thetam,p);
    else
        Q=uFunc(thetaM,p)+uFunc(thetam,p);
    end
else
    theta1InRange2=(theta1>=-pi&&theta1<=-pi/2)||(theta1>=pi/2&&theta1<=pi);
    theta2InRange2=(theta2>=-pi&&theta2<=-pi/2)||(theta2>=pi/2&&theta2<=pi);
    if(theta1InRange2&&theta2InRange1||theta2InRange2&&theta1InRange1)
        if(sign(theta1)==sign(theta2))
            Q=pi/2-uFunc(pi-thetaM,p)-uFunc(thetam,p);
        else
            %The paper seemed to have a typo on the sign of the first u
            %function for this one.
            Q=pi/2+uFunc(pi-thetaM,p)-uFunc(thetam,p);
        end
    else
        if(sign(theta1)==sign(theta2))
            Q=uFunc(pi-thetam,p)-uFunc(pi-thetaM,p);
        else
            Q=pi-uFunc(pi-thetaM,p)-uFunc(pi-thetam,p);
        end
    end
end

M=Q/pi;

%The paper does not consider the case of a sector going more than pi
%around. This covers it.
if(isLongWayAround)
    M=1-M;
end
end

function val=uFunc(theta,p)
    val=(1/2)*atan(1/p*tan(theta));
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
