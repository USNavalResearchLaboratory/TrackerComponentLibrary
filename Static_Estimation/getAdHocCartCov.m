function V = getAdHocCartCov(bandwidth,beamwidth,snr,x,dim)
%%GETADHOCCARTCOV Computes a covariance matrix in 2D or 3D using radar
%                 sensor parameters and optional target location.
%
%INPUT:
% bandwidth: A scalar bandwidth for the radar sensor.
% beamwidth: A scalar beamwidth for the radar sensor (both azimuth and
%            elevation beamwidth are equal) or a length 2 vector with the
%            beamwidth ordering [azimuth, elevation].
% snr: A scalar signal-to-noise ratio expressed in dB.
% x: A 2-by-1 or 3-by-1 vector representing the estimated target location
%    in Cartesian coordinates. Defaults to [1;0;0];
% dim: An integer, either 2 or 3, which indicates whether the radar
%      measures in polar (range and azimuth) or spherical (range, azimuth,
%      and elevation) coordinates. Note that x is expected to be in
%      Cartesian coordinates. Defaults to the dimensionality of x.
%
%OUTPUT:
% V: A 2-by-2 or 3-by-3 covariance matrix.
%
%EXAMPLE: Constructs an ad-hoc measurement ellipse for a estimated target
%         location given bandwidth, beamwidth, and SNR.
% bandwidth = 5e6;
% beamwidth = [deg2rad(2),deg2rad(10)];
% snr = 10;
% x = [1e3;1e3;1e3];
% V = getAdHocCartCov(bandwidth,beamwidth,snr,x);
% drawEllipse(x,V\eye(3))
% hold on
% plot3([0,x(1)],[0,x(2)],[0,x(3)],'r','LineWidth',2)
% hold off
% axis('equal')
%
%December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
if nargin<4
    x = [1;0;0];
end
if nargin<5
    dim = size(x,1);
end
if length(beamwidth)==2
    azBeamwidth = beamwidth(1);
    elBeamwidth = beamwidth(2);
else
    azBeamwidth = beamwidth(1);
    elBeamwidth = azBeamwidth;
end
r = norm(x(1:dim));
rangeRes = Constants.speedOfLight/(2*bandwidth*sqrt(2*snr));
azAngleRes = 2*r*sin(azBeamwidth/2);
elAngleRes = 2*r*sin(elBeamwidth/2);
R2 = rotAxis2Vec(x(1:dim),'x');
if dim==3
    V = [(rangeRes/2)^2,0,0;
        0,(azAngleRes/2)^2,0;
        0,0,(elAngleRes/2)^2];
    R1 = rotMat2D(pi-atan2(R2(3,2),R2(3,3)));
    V(2:3,2:3) = R1*V(2:3,2:3)*R1';
elseif dim==2
    V = [(rangeRes/2)^2,0;
        0,(azAngleRes/2)^2];
end
V = R2*V*R2';
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