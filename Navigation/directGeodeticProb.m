function [latLonEnd,azEnd]=directGeodeticProb(latLonStart,azStart,distVal,a,f)
%%DIRECTGEODETICPROB Solve the direct geodetic problem. That is, given an
%                    initial point and an initial bearing on an ellipsoidal
%                    Earth, find the end point and final bearing if one
%                    were to travel one or more given distances along a
%                    geodesic curve (the shortest curve between two points
%                    on a curved surface).
%
%INPUTS: latLonStart The 2X1 initial point given in geodetic latitude and
%                   longitude in radians of the format
%                   [latitude;longitude]. The latitude must be between
%                   -pi/2 and pi/2 radians and the longitude between -pi
%                   and pi radians.
%           azStart The forward azimuth (initial heading) at the starting
%                   point in radians East of true North on the reference
%                   ellipsoid.
%              dist An NX1 or 1XN vector of the distances in meters 
%                   that will be traveled on the geodesic curves starting
%                   at latLonStart with initial heading azStart where
%                   solutions are desired.
%                 a The semi-major axis of the reference ellipsoid (in
%                   meters). If this argument is omitted or an empty matrix
%                   is passed, the value in Constants.WGS84SemiMajorAxis is
%                   used.
%                 f The flattening factor of the reference ellipsoid. If
%                   this argument is omitted or an empty matrix is passed,
%                   the value in Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonEnd A 2XN matrix of geodetic latitude and longitudes of the
%                   final points of the geodesic trajectory given in
%                   radians as [latitude;longitude].
%             azEnd An NX1 vector of the forward azimuth (bearing) at the
%                   ending points in radians East of true North on the
%                   reference ellipsoid.
%
%The algorithm of [1] is used. In [1], series terms up to order 6 are given
%for a number of variables. In this implementation, values to the tenth
%order are given, though with double precision arithmetic it doesn't make a
%difference. Note that values to the 30th order have been posted online at
%https://geographiclib.sourceforge.io/html/geodseries30.html
%However, there is generally no point in using higher orders unless one
%uses higher than double-precision floating point arithmetic. The series
%are only valid for small values of f, which is appropriate when dealing
%with planets such as the Earth that are approximated with ellipsoids
%having a low eccentricity.
%
%EXAMPLE:
%We show that indirectGeodeticProb and directGeodeticProb are consistent
%with each toher.
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% [azStart,dist]=indirectGeodeticProb(latLonStart,latLonEnd);
% latLonPoint=directGeodeticProb(latLonStart,azStart,dist);
% max(abs(latLonPoint-latLonEnd))
%One will see that the two values are the same within finite-precision
%limits, as expected.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<4||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

[latLonEnd,azEnd]=directGeodeticKarney(latLonStart,azStart,distVal,a,f);

end

function [latLonEnd,azEnd]=directGeodeticKarney(latLonStart,azStart,distVal,a,f)
%%DIRECTGEODETICKARNEY This function implements the algorithm to solve the
%            direct geodetic problem as given in [1] with the series
%            expansions extended to the tenth order. This method should be
%            numerically stable everywhere.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 1
b=a*(1-f);
%Equation 2
n=f/(2-f);
%From Equations 24 and 25.
[A3Coeffs,C3Coeffs]=computeA3C3Coeffs(n);
%Equation 3
e2=f*(2-f);
%Equation 4
ePrime2=e2/(1-e2);

N=length(distVal);
latLonEnd=zeros(2,N);
azEnd=zeros(N,1);

%Change to the notation in [1].
phi1=latLonStart(1);
alpha1=azStart;
sinAlpha1=sin(alpha1);
cosAlpha1=cos(alpha1);

%Equation 6.
%We could just do
%beta1=atan((1-f)*tan(phi1))
%and 
%sinBeta1=sin(beta1);
%cosBeta1=cos(beta1);
%but this is numerically stabler around +/- pi/2.
%sinBeta1 and cosBeta1 are here a non-normalized vector for the sine and
%cosine.
sinBeta1=(1-f)*sin(phi1);
cosBeta1=max(sqrt(realmin),cos(phi1));
normVal=hypot(sinBeta1,cosBeta1);
sinBeta1=sinBeta1/normVal;
cosBeta1=cosBeta1/normVal;

%Solve triangle NEA
%Equation 10
cosAlpha0=abs(cosAlpha1+1j*sinAlpha1*sinBeta1);
sinAlpha0=sinAlpha1*cosBeta1;

%Equation 11
cosSigma1=cosAlpha1*cosBeta1;
sinSigma1=sinBeta1;
normVal=norm([cosSigma1;sinSigma1]);
cosSigma1=cosSigma1/normVal;
sinSigma1=sinSigma1/normVal;

sigma1=atan2(sinSigma1,cosSigma1);
%Use a double angle identity
cos2Sigma1=1-2*sinSigma1^2;
sin2Sigma1=2*sinSigma1*cosSigma1;

%Equation 12
omega1=angle(cosSigma1+1j*sinAlpha0*sinSigma1);

%Determine sigma2
%Equation 9
k2=ePrime2*cosAlpha0^2;
%Equation 16
temp=sqrt(1+k2);
epsilon=(temp-1)/(temp+1);

%Equations 17 and 18
[A1,C1]=computeA1C1(epsilon);
%Equation 15
I1Sigma1=A1*(sigma1+evalSinCosSeries(C1,cos2Sigma1,sin2Sigma1));
%Equation 21
C1p=computeC1p(epsilon);
for curPoint=1:N
    s12=distVal(curPoint);

    %Equation 7
    s1=I1Sigma1*b;
    s2=s1+s12;

    tau2=s2/(b*A1);

    %Equation 20
    sigma2=tau2+evalSinCosSeries(C1p,cos(2*tau2),sin(2*tau2));
   
    %Note that the above solution for sigma2 could also have been solved
    %using an elliptic integral of the second kind and an iterative
    %solution as follows:
%     %Directly evaluate the elliptic integral in equation 7.
%     I1Sigma1=ellipIntInc2Kind(sigma1,-k2);
% 
%     %Equation 7
%     s1=I1Sigma1*b;
%     s2=s1+s12;
% 
%     %Use Newton's method to invert Equation 7 for a value sigma2 that
%     %results in the solution s2/b. The cost function that we want to
%     %zero is
%     %(s2-b*I1(sigma2))^2
%     %We just iterate until convergence.
% 
%     sigma2=sigma1;
%     curIter=0;
%     while(1)
%         I1Sigma2=ellipIntInc2Kind(sigma2,-k2);
%         fVal=(s2-I1Sigma2*b)^2;
%         derivVal=-2*b*sqrt(1+k2*sin(sigma2)^2)*(s2-b*I1Sigma2);
%         deltaSigma2=-fVal/derivVal;
%         sigma2=sigma2+deltaSigma2;
% 
%         curIter=curIter+1;
%         if(abs(deltaSigma2)<=eps(sigma2))
%             break;
%         end
%     end

    sinSigma2=sin(sigma2);
    cosSigma2=cos(sigma2);
    %Use a double angle identity
    cos2Sigma2=1-2*sinSigma2^2;
    sin2Sigma2=2*sinSigma2*cosSigma2;

    %Solve triangle NEB
    %Equation 14
    %alpha2=angle(cosAlpha0*cosSigma2+1j*sinAlpha0);
    alpha2=atan2(sinAlpha0,cosAlpha0*cosSigma2);
    %Equation 13
    %beta2=angle(abs(cosAlpha0*cosSigma2+1j*sinAlpha0)+1j*cosAlpha0*sinSigma2);
    beta2=atan2(cosAlpha0*sinSigma2,sqrt((cosAlpha0*cosSigma2)^2+sinAlpha0^2));
    %Equation 12
    %omega2=angle(cosSigma2+1j*sinAlpha0*sinSigma2);
    omega2=atan2(sinAlpha0*sinSigma2,cosSigma2);
    
    %Determine lambda12
    [A3,C3]=computeA3C3(A3Coeffs,C3Coeffs,epsilon);
    %Equation 23
    I3Sigma1=A3*(sigma1+evalSinCosSeries(C3,cos2Sigma1,sin2Sigma1));
    I3Sigma2=A3*(sigma2+evalSinCosSeries(C3,cos2Sigma2,sin2Sigma2));
    %Equation 8
    lambda1=omega1-f*sinAlpha0*I3Sigma1;
    lambda2=omega2-f*sinAlpha0*I3Sigma2;
    lambda12=lambda2-lambda1;

    %Solution
    %Equation 6
    phi2=atan(tan(beta2)/(1-f));

    latLonEnd(:,curPoint)=[phi2;wrapRange(latLonStart(2)+lambda12,-pi,pi)];
    
    azEnd(curPoint)=alpha2;
end
end

function [A1,C1]=computeA1C1(e)
%COMPUTEA1C1 This function computes A1 and C1 as in Equation 17 of [1], but
%            the series have been extended to the tenth order.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Coefficients for A1 as in Equation 17, but extended to the tenth order.
A1Coeffs=[1;
          1/4;
          1/64;
          1/256;
          25/16384;
          49/65536];
e2=e*e;
A1=polyValGen(A1Coeffs,e2,1)/(1-e);

%Coefficients for the C1 terms as in Equation 17, but extended to the tenth
%order.
C1Coeffs=[-1/2,    -1/16,   -1/48,     -5/512,    -7/1280, -7/2048,     -33/14336, -429/262144, -715/589824, -2431/2621440;
           3/16,    1/32,    3/256,     3/512,     7/2048,  9/4096,      99/65536,  143/131072,  0,           0;
          -1/32,   -9/2048, -3/2048,   -11/16384, -3/8192, -117/524288,  0,         0,           0,           0;
           19/2048, 7/4096,  17/24576,  3/8192,    0,       0,           0,         0,           0,           0;
          -3/4096,  1/6553,  0,         0,         0,       0,           0,         0,           0,           0];
      
C1=zeros(10,1);
ePow=e;
C1(1)=ePow*polyValGen(C1Coeffs(:,1),e2,1);
ePow=e2;
C1(2)=ePow*polyValGen(C1Coeffs(:,2),e2,1);
ePow=ePow*e;
C1(3)=ePow*polyValGen(C1Coeffs(1:4,3),e2,1);
ePow=ePow*e;
C1(4)=ePow*polyValGen(C1Coeffs(1:4,4),e2,1);
ePow=ePow*e;
C1(5)=ePow*polyValGen(C1Coeffs(1:3,5),e2,1);
ePow=ePow*e;
C1(6)=ePow*polyValGen(C1Coeffs(1:3,6),e2,1);
ePow=ePow*e;
C1(7)=ePow*polyValGen(C1Coeffs(1:2,7),e2,1);
ePow=ePow*e;
C1(8)=ePow*polyValGen(C1Coeffs(1:2,8),e2,1);
ePow=ePow*e;
C1(9)=ePow*polyValGen(C1Coeffs(1,9),e2,1);
ePow=ePow*e;
C1(10)=ePow*polyValGen(C1Coeffs(1,10),e2,1);
end

function C1p=computeC1p(e)
%COMPUTEC1p This function computes C1' as in Equation 21 of [1], but the
%           series have been extended to the tenth order.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

C1pCoeffs=[ 1/2,          5/16,             29/96,          539/1536,        3467/7680,       38081/61440,       459485/516096,  109167851/82575360, 83141299/41287680, 9303339907/2972712960;
           -9/32,        -37/96,           -75/128,        -2391/2560,      -28223/18432,    -733437/286720,    -709743/163840, -550835669/74317824, 0,                 0;
            205/1536,     1335/4096,        2901/4096,      1082857/737280,  1361343/458752,  10820079/1835008,  0,              0,                  0,                 0;
           -4879/73728,  -86171/368640,    -443327/655360, -2722891/1548288, 0,               0,                 0,              0,                  0,                 0;
            9039/327680,  4119073/28311552, 0,              0,               0,               0,                 0,              0,                  0,                 0];

e2=e*e;

C1p=zeros(10,1);
ePow=e;
C1p(1)=ePow*polyValGen(C1pCoeffs(:,1),e2,1);
ePow=ePow*e;
C1p(2)=ePow*polyValGen(C1pCoeffs(:,2),e2,1);
ePow=ePow*e;
C1p(3)=ePow*polyValGen(C1pCoeffs(1:4,3),e2,1);
ePow=ePow*e;
C1p(4)=ePow*polyValGen(C1pCoeffs(1:4,4),e2,1);
ePow=ePow*e;
C1p(5)=ePow*polyValGen(C1pCoeffs(1:3,5),e2,1);
ePow=ePow*e;
C1p(6)=ePow*polyValGen(C1pCoeffs(1:3,6),e2,1);
ePow=ePow*e;
C1p(7)=ePow*polyValGen(C1pCoeffs(1:2,7),e2,1);
ePow=ePow*e;
C1p(8)=ePow*polyValGen(C1pCoeffs(1:2,8),e2,1);
ePow=ePow*e;
C1p(9)=ePow*polyValGen(C1pCoeffs(1,9),e2,1);
ePow=ePow*e;
C1p(10)=ePow*polyValGen(C1pCoeffs(1,10),e2,1);
end

function [A3Coeffs,C3Coeffs]=computeA3C3Coeffs(n)
%COMPUTEA3C3CCOEFFS This function computes the coefficients of the
%            polynomials in terms of epsilon for A3 and C3 as in Equations
%            24 and 25 of [1]. The series have been extended to the tenth
%            order.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A3nCoeffs=-[-1,         0,           0,            0,          0,          0,            0,            0,           0,             0,            0;
            1/2,       -1/2,         0,            0,          0,          0,            0,            0,           0,             0,            0;
            1/4,        1/8,        -3/8,          0,          0,          0,            0,            0,           0,             0,            0;
            1/16,       3/16,        1/16,        -5/16,       0,          0,            0,            0,           0,             0,            0;
            3/64,       1/32,        5/32,         5/128,     -35/128,     0,            0,            0,           0,             0,            0;
            3/128,      5/128,       5/256,        35/256,     7/256,     -63/256,       0,            0,           0,             0,            0;
            5/256,      15/1024,     35/1024,      7/512,      63/512,     21/1024,     -231/1024,     0,           0,             0,            0;
            25/2048,    35/2048,     21/2048,      63/2048,    21/2048,    231/2048,     33/2048,     -429/2048,    0,             0,            0;
            175/16384,  35/4096,     63/4096,      63/8192,    231/8192,   33/4096,      429/4096,     429/32768,  -6435/32768,    0,            0;
            245/32768,  315/32768,   105/16384,    231/16384,  99/16384,   429/16384,    429/65536,    6435/65536,  715/65536,    -12155/65536,  0;
            441/65536,  735/131072,  1155/131072,  165/32768,  429/32768,  1287/262144,  6435/262144,  715/131072,  12155/131072,  2431/262144, -46189/262144]';

%Equation 24
A3Coeffs=zeros(11,1);
for k=1:11
    A3Coeffs(k)=polyValGen(A3nCoeffs(1:k,k),n,1);
end

C3Coeffs=zeros(10,10);
C31nCoeffs=[1/4,         -1/4,         0,            0,          0,           0,          0,           0,           0,            0,            0;
            1/8,          0,          -1/8,          0,          0,           0,          0,           0,           0,            0,            0;
            3/64,         3/64,       -1/64,        -5/64,       0,           0,          0,           0,           0,            0,            0;
            5/128,        1/64,        1/64,        -1/64,      -7/128,       0,          0,           0,           0,            0,            0;
            3/128,        11/512,      3/512,        1/256,     -7/512,      -21/512,     0,           0,           0,            0,            0;
            21/1024,      5/512,       13/1024,      1/512,     -1/1024,     -3/256,     -33/1024,     0,           0,            0,            0;
            243/16384,    189/16384,   83/16384,     127/16384,  3/16384,    -51/16384,  -165/16384,  -429/16384,   0,            0,            0;
            435/32768,    109/16384,   1/128,        45/16384,   39/8192,    -11/16384,  -33/8192,    -143/16384,  -715/32768,    0,            0;
            345/32768,    953/131072,  259/65536,    365/65536,  95/65536,    47/16384,  -143/131072, -143/32768,  -1001/131072, -2431/131072,  0;
            2511/262144,  317/65536,   1355/262144,  165/65536,  531/131072,  89/131072,  107/65536,  -169/131072, -1157/262144, -221/32768,   -4199/262144]';
for k=1:10
    C3Coeffs(k,1)=polyValGen(C31nCoeffs(1:(k+1),k),n,1);
end

C32nCoeffs=[1/16,        -3/32,         1/32,        0,           0,             0,             0,            0,            0,           0,            0;
            3/64,        -1/32,        -3/64,        1/32,        0,             0,             0,            0,            0,           0,            0;
            3/128,        1/128,       -9/256,      -3/128,       7/256,         0,             0,            0,            0,           0,            0;
            5/256,        1/256,       -1/128,      -7/256,      -3/256,         3/128,         0,            0,            0,           0,            0;
            27/2048,      69/8192,     -39/8192,    -47/4096,    -41/2048,      -45/8192,       165/8192,     0,            0,           0,            0;
            187/16384,    39/8192,      31/16384,   -63/8192,    -185/16384,    -119/8192,     -33/16384,     143/8192,     0,           0,            0;
            287/32768,    47/8192,      31/65536,   -3/2048,     -537/65536,    -41/4096,      -693/65536,    0,            1001/65536,  0,            0;
            255/32768,    249/65536,    43/16384,   -119/65536,  -25/8192,      -507/65536,    -35/4096,     -507/65536,    39/32768,    221/16384,    0;
            1675/262144,  2127/524288,  753/524288,  357/524288, -3109/1048576, -3873/1048576, -1821/262144, -1885/262144, -2977/524288, 1989/1048576, 12597/1048576]';
for k=1:9
    C3Coeffs(k,2)=polyValGen(C32nCoeffs(1:(k+2),k),n,1);
end

C33nCoeffs=[5/192,       -3/64,        5/192,    -1/192,        0,            0,            0,            0,            0,            0,            0;
            3/128,       -5/192,      -1/64,      5/192,       -1/128,        0,            0,            0,            0,            0,            0;
            7/512,       -1/384,      -77/3072,   5/3072,       65/3072,     -9/1024,       0,            0,            0,            0,            0;
            3/256,       -1/1024,     -71/6144,  -47/3072,      9/1024,       25/1536,     -55/6144,      0,            0,            0,            0;
            139/16384,    143/49152,  -383/49152, -179/16384,  -121/16384,    547/49152,    605/49152,   -143/16384,    0,            0,            0;
            243/32768,    95/49152,   -41/16384,  -147/16384,  -389/49152,   -109/49152,    557/49152,    455/49152,   -273/32768,    0,            0;
            581/98304,    377/131072, -33/16384,  -907/196608, -515/65536,   -1937/393216,  89/98304,     2093/196608,  455/65536,   -1547/196608,  0;
            1383/262144,  103/49152,  -17/262144, -127/32768,  -3853/786432, -25/4096,     -2011/786432,  265/98304,    1895/196608,  85/16384,    -969/131072]';
for k=1:8
    C3Coeffs(k,3)=polyValGen(C33nCoeffs(1:(k+3),k),n,1);
end

C34nCoeffs=[7/512,       -7/256,     5/256,       -7/1024,       1/1024,      0,            0,            0,            0,            0,           0;
            7/512,       -5/256,    -7/2048,       9/512,       -21/2048,     1/512,        0,            0,            0,            0,           0;
            9/1024,      -43/8192,  -129/8192,     39/4096,      91/8192,    -91/8192,      11/4096,      0,            0,            0,           0;
            127/16384,   -23/8192,  -165/16384,   -47/8192,      213/16384,   11/2048,     -175/16384,    13/4096,      0,            0,           0;
            193/32768,    3/8192,   -505/65536,   -227/32768,    75/65536,    801/65536,    165/131072,  -637/65536,    455/131072,   0,           0;
            171/32768,    25/65536, -259/65536,   -471/65536,   -351/131072,  605/131072,   41/4096,     -189/131072,  -1127/131072,  119/32768,   0;
            1121/262144, 339/262144, -801/262144, -2525/524288, -2519/524288, 73/131072,    1539/262144,  1989/262144, -1633/524288, -3927/524288, 969/262144]';
for k=1:7
    C3Coeffs(k,4)=polyValGen(C34nCoeffs(1:(k+4),k),n,1);
end

C35nCoeffs=[21/2560,     -9/512,       15/1024,      -7/1024,    9/5120,       -1/5120,       0,             0,            0,           0,          0;
            9/1024,      -15/1024,     3/2048,        57/5120,  -5/512,         9/2560,      -1/2048,        0,            0,           0,          0;
            99/16384,    -91/16384,   -781/81920,     883/81920, 319/81920,    -783/81920,    387/81920,    -13/16384,     0,           0,          0;
            179/32768,   -55/16384,   -79/10240,     -27/81920,  461/40960,    -139/81920,   -65/8192,       441/81920,   -35/32768,    0,          0;
            141/32768,   -109/131072, -217/32768,    -219/65536, 1559/327680,   5431/655360, -203/40960,    -1943/327680,  369/65536,  -85/65536,   0;
            1013/262144, -15/32768,   -5399/1310720, -199/40960, 1267/1310720,  1007/163840,  6277/1310720, -527/81920,   -659/163840,  459/81920, -969/655360]';
for k=1:6
    C3Coeffs(k,5)=polyValGen(C35nCoeffs(1:(k+5),k),n,1);
end

C36nCoeffs=[11/2048,     -99/8192,       275/24576,     -77/12288,     9/4096,        -11/24576,      1/24576,        0,             0,             0,            0;
            99/16384,    -275/24576,     55/16384,       167/24576,   -407/49152,      35/8192,      -55/49152,       1/8192,        0,             0,            0;
            143/32768,   -253/49152,    -1105/196608,    481/49152,   -73/196608,     -169/24576,     1067/196608,   -11/6144,       15/65536,      0,            0;
            33/8192,     -221/65536,    -23/4096,        457/196608,   267/32768,     -329/65536,    -69/16384,       375/65536,    -77/32768,      17/49152,     0;
            1711/524288, -4333/3145728, -16885/3145728, -1343/1572864, 17381/3145728,  8519/2097152, -42985/6291456, -4885/3145728,  8549/1572864, -5797/2097152, 969/2097152]';
for k=1:5
    C3Coeffs(k,6)=polyValGen(C36nCoeffs(1:(k+6),k),n,1);
end

C37nCoeffs=[429/114688, -143/16384,  143/16384, -91/16384,   39/16384,    -11/16384,       13/114688,    -1/114688,       0,            0,           0;
            143/32768,  -143/16384,  65/16384,   65/16384,  -109/16384,    507/114688,    -27/16384,      39/114688,     -1/32768,      0,           0;
            429/131072, -299/65536, -13/4096,    269/32768, -601/229376,  -989/229376,     9475/1835008, -4667/1835008,  1157/1835008, -17/262144,   0;
            403/131072, -13/4096,   -521/131072, 393/114688, 1209/229376, -11001/1835008, -3979/3670016,  8821/1835008, -833/262144,    429/458752, -57/524288]';
for k=1:4
    C3Coeffs(k,7)=polyValGen(C37nCoeffs(1:(k+7),k),n,1);
end  

C38nCoeffs=[715/262144, -429/65536,     455/65536,   -637/131072,   315/131072,   -55/65536,      13/65536,      -15/524288,    1/524288,      0,           0;
            429/131072, -455/65536,     1053/262144,  35/16384,    -1361/262144,   69/16384,     -2095/1048576,   77/131072,   -105/1048576,   1/131072,    0;
            663/262144, -4173/1048576, -1717/1048576, 3485/524288, -3825/1048576, -9469/4194304,  18469/4194304, -6137/2097152, 4455/4194304, -885/4194304, 19/1048576]';
for k=1:3
    C3Coeffs(k,8)=polyValGen(C38nCoeffs(1:(k+8),k),n,1);
end  

C39nCoeffs=[2431/1179648, -663/131072,  1105/196608, -833/196608,  153/65536, -187/196608,   221/786432,   -15/262144,     17/2359296, -1/2359296,  0;
            663/262144,   -1105/196608, 1003/262144,  187/196608, -391/98304,  1003/262144, -3425/1572864,  1921/2359296, -13/65536,    17/589824, -1/524288]';
for k=1:2
    C3Coeffs(k,9)=polyValGen(C39nCoeffs(1:(k+9),k),n,1);
end

C310nCoeffs=[4199/2621440, -4199/1048576, 4845/1048576, -969/262144, 2907/1310720, -10659/10485760, 741/2097152, -95/1048576, 17/1048576, -19/10485760, 1/10485760]';
for k=1:1
    C3Coeffs(k,10)=polyValGen(C310nCoeffs(1:(k+10),k),n,1);
end

end

function [A3,C3]=computeA3C3(A3Coeffs,C3Coeffs,e)
%COMPUTEA3C3 This function computes A3 and C3 as in Equations 24 and 25 of
%            [1], but the series have been extended to the tenth order. The
%            coefficients for A3 and C3 in terms of n must have been
%            already evaluated.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A3=polyValGen(A3Coeffs,e,1);

C3=zeros(10,1);
ePow=e;
C3(1)=ePow*polyValGen(C3Coeffs(1:10,1),e,1);
ePow=ePow*e;
C3(2)=ePow*polyValGen(C3Coeffs(1:9,2),e,1);
ePow=ePow*e;
C3(3)=ePow*polyValGen(C3Coeffs(1:8,3),e,1);
ePow=ePow*e;
C3(4)=ePow*polyValGen(C3Coeffs(1:7,4),e,1);
ePow=ePow*e;
C3(5)=ePow*polyValGen(C3Coeffs(1:6,5),e,1);
ePow=ePow*e;
C3(6)=ePow*polyValGen(C3Coeffs(1:5,6),e,1);
ePow=ePow*e;
C3(7)=ePow*polyValGen(C3Coeffs(1:4,7),e,1);
ePow=ePow*e;
C3(8)=ePow*polyValGen(C3Coeffs(1:3,8),e,1);
ePow=ePow*e;
C3(9)=ePow*polyValGen(C3Coeffs(1:2,9),e,1);
ePow=ePow*e;
C3(10)=ePow*polyValGen(C3Coeffs(1:1,10),e,1);
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
