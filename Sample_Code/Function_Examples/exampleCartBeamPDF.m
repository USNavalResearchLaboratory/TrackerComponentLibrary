function pxy = exampleCartBeamPDF(rmin,rmax,bmean,bstd,x,y)
%%EXAMPLECARTBEAMPDF An example probability density which roughly mimics
%                    a range/bearing sensor measurement using a smoothed
%                    uniform density on a range interval and Gaussian
%                    bearing density which are assumed uncorrelated.
%
%INPUTS:
% rmin: The range interval's minimum value. This must be greater than or
%       equal to zero, finite, and less than rmax.
% rmax: The range interval's maximum value. This must be greater than rmin
%       and finite.
% bmean: The mean bearing for the beam. This must be finite.
% bstd: The standard deviation for the beam. This must be positive and
%       finite.
% x: An array of x-coordinates. The size must agree with y for
%    element-wise operations.
% y: An array of y-coordinates. The size must agree with x for
%    element-wise operations.
%
%OUTPUTS:
% pxy: An array of probability density values at the points (x,y). The
%      size of the array is size(x).
%
%EXAMPLE: Plots the PDF for a sensor beam.
% sensor = [-75;-50];
% sensorbeam = @(x,y) exampleCartBeamPDF(50,150,deg2rad(45),deg2rad(5),x-sensor(1),y-sensor(2));
% x = linspace(-100,100,1e3);
% [X,Y] = meshgrid(x);
% P = sensorbeam(X,Y);
% contourf(X,Y,P)
% hold on
% plot(sensor(1),sensor(2),'ro',"MarkerFaceColor",'r')
% hold off
%
%REFERENCES:
%[1] W. J. Farrell and C. Ganesh, "Generalized Chernoff fusion
%    approximation for practical distributed data fusion," in 2009 12th
%    International Conference on Information Fusion, 2009, pp. 555-562.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if rmin<0 || ~isfinite(rmin)
    error("rmin must be a finite positive number.")
end
if rmin>rmax || ~isfinite(rmax)
    error("rmax must be a finite number greater than rmin.")
end
if ~isfinite(bmean)
    error("bmean must be finite.")
end
if bstd<=0 || ~isfinite(bstd)
    error("bstd must be a finite positive number.")
end
if ~all(size(x)==size(y))
    error("The shape of x must be the same as the shape of y.")
end

r = sqrt(x.^2+y.^2);
b = mod(atan2(y,x),2*pi);

pxy = beamRangePDF(rmin,rmax,r).*beamBearingPDF(bmean,bstd,b);
end

function pr = beamRangePDF(rmin,rmax,r)
%%BEAMRANGEPDF Computes the probability density for the range component
%              assuming a smoothed uniform distribution as defined in [1].
%
%INPUTS:
% rmin: The range interval's minimum value. This must be greater than or
%       equal to zero, finite, and less than rmax.
% rmax: The range interval's maximum value. This must be greater than rmin
%       and finite.
% r: The array of range values at which the density should be evaluated.
%
%OUTPUTS:
% pr: The density values corresponding to the given values in r.
%
%REFERENCES:
%[1] W. J. Farrell and C. Ganesh, "Generalized Chernoff fusion
%    approximation for practical distributed data fusion," in 2009 12th
%    International Conference on Information Fusion, 2009, pp. 555-562.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    pr = (atan(rmax-r)+atan(r-rmin))/((rmin+rmax)*pi);
end

function pb = beamBearingPDF(bmean,bstd,b)
%%BEAMBEARINGPDF Computes the probability density for the bearing component
%                assuming a Gaussian distribution as defined in [1].
%
%INPUTS:
% bmean: The mean bearing for the beam. This must be finite.
% bstd: The standard deviation for the beam. This must be positive and
%       finite.
% b: The array of bearing values at which the density should be evaluated.
%
%REFERENCES:
%[1] W. J. Farrell and C. Ganesh, "Generalized Chernoff fusion
%    approximation for practical distributed data fusion," in 2009 12th
%    International Conference on Information Fusion, 2009, pp. 555-562.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
pb = WrappedNormalD.PDF(b(:)',bmean,bstd.^2);
pb = reshape(pb,size(b));
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
