function demoMeasConvFunctions
%%DEMOMEASCONVFUNCTIONS Demonstrate functions that take Gaussian-corrupted
%           measurements (mean and covariance matrix) and obtain the mean
%           and covariance matrix of the measurements converted into
%           Cartesian coordinates. Two types of conversion approaches
%           dominate the literature: Those using Taylor series
%           approximations and those using cubature integration. The
%           cubature integration methods are theoretically more accurate.
%           However, at low noise levels, the two approaches are
%           essentially equivalent.
%
%The problem of obtaining the mean and covariance matrix of a target
%measurement in Cartesian coordinates can be easily solved using cubature
%integration, as dicussed in [1]. However, more traditional approaches to
%the problem utilize Taylor expansions. Two examples of using Taylor series
%expansions for measurement conversion are in Chapter 10.4.3 of [2].
%
%This function compares various Taylor series expansions (the traditional
%approach) to cubature integration utilizing fifth-order cubature points.
%Three scenarios are considered: polar measurements, monostatic spherical
%measurements, and monostatic r-u-v measurements. This function essentially
%compares the performance of the functions pol2CartTaylor,
%pol2CartCubature, spher2CartTaylor, spher2CartCubature,
%monostatRuv2CartTaylor, and ruv2CartCubature.
%
%As the root mean squared error (RMSE) for the cubature method tends to be
%lowest, plots of RMSE versus the RMSE of the cubature method are generated
%along with plots of the  average normalized estimation error squared
%(NEES) are shown. The NEES is a method of quanifying how well the
%covariance matrices agree with the observed accuracy. Bars for the 95%
%confidence region are drawn.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('Example 1: Polar Measurement Conversion')

%We consider measurements at 45 degrees azimuth at ranges from 1 to 70km
%from the sensor.
numPoints=50;
numRuns=200;

%50 meter range standard deviation
sigmaR=50;
%10 degree angular standard deviation.
sigmaTheta=10*(pi/180);
SR=diag([sigmaR,sigmaTheta]);
R=SR*SR';%Measurement variance

az=pi/4;
dirTarget=[cos(az);sin(az)];
tarDist=linspace(1e3,70e3,numPoints);
locList=[dirTarget(1)*tarDist;dirTarget(2)*tarDist];

numMethods=5;
RMSE=zeros(numPoints,numMethods);
NEES=zeros(numPoints,numMethods);

[xi,w]=fifthOrderCubPoints(2);

for curPoint=1:numPoints
    zCartTrue=locList(:,curPoint);
    z=Cart2Pol(zCartTrue);
    
    MSECur=zeros(1,numMethods);
    NEESCur=zeros(1,numMethods);
    for curRun=1:numRuns 
        zMeas=z+SR*randn(2,1);

        [zCart,RCart]=pol2CartTaylor(zMeas,R,[],[],0);
        diff=zCartTrue-zCart;
        MSECur(1)=MSECur(1)+diff'*diff;
        NEESCur(1)=NEESCur(1)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=pol2CartTaylor(zMeas,R,[],[],1);
        diff=zCartTrue-zCart;
        MSECur(2)=MSECur(2)+diff'*diff;
        NEESCur(2)=NEESCur(2)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=pol2CartTaylor(zMeas,R,[],[],2);
        diff=zCartTrue-zCart;
        MSECur(3)=MSECur(3)+diff'*diff;
        NEESCur(3)=NEESCur(3)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=pol2CartTaylor(zMeas,R,[],[],3);
        diff=zCartTrue-zCart;
        MSECur(4)=MSECur(4)+diff'*diff;
        NEESCur(4)=NEESCur(4)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=pol2CartCubature(zMeas,SR,[],[],[],[],[],xi,w);
        diff=zCartTrue-zCart;
        MSECur(5)=MSECur(5)+diff'*diff;
        NEESCur(5)=NEESCur(5)+diff'*inv(RCart)*diff;
    end
    RMSE(curPoint,:)=sqrt(MSECur/numRuns);
    NEES(curPoint,:)=NEESCur/(numRuns*2);
end

figure(1)
clf
hold on
plot(tarDist/1e3,RMSE(:,1)./RMSE(:,5),'-r','linewidth',2)%Standard conversion
plot(tarDist/1e3,RMSE(:,2)./RMSE(:,5),'--b','linewidth',2)%Additive
plot(tarDist/1e3,RMSE(:,3)./RMSE(:,5),':g','linewidth',2)%Multiplicative
plot(tarDist/1e3,RMSE(:,4)./RMSE(:,5),'-c','linewidth',4)%Modified Multiplicative
plot(tarDist/1e3,RMSE(:,5)./RMSE(:,5),'--m','linewidth',2)%Cubature

h1=xlabel('Distance (km)');
h2=ylabel('RMSE Ratio');
title('RMSE Ratio of Conversion to Cubature Conversion (Polar Meas.)')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('Standard Conversion','Additive','Multiplicative', 'Modified Multiplicative','Cubature','Location','NorthEast') 

figure(2)
clf
hold on
plot(tarDist/1e3,NEES(:,1),'-r','linewidth',2)%Standard conversion
plot(tarDist/1e3,NEES(:,2),'--b','linewidth',2)%Additive
plot(tarDist/1e3,NEES(:,3),':g','linewidth',2)%Multiplicative
plot(tarDist/1e3,NEES(:,4),'-c','linewidth',4)%Modified Multiplicative
plot(tarDist/1e3,NEES(:,5),'--m','linewidth',2)%Cubature

%The upper and lower bounds marking a 95% confidence region.
UB=ChiSquareD.invCDF(1-0.05/2,numRuns*2)/(numRuns*2);
LB=ChiSquareD.invCDF(0.05/2,numRuns*2)/(numRuns*2);

%Draw the bounds
plot([tarDist(1)/1e3,tarDist(end)/1e3],[UB,UB],'-k')
plot([tarDist(1)/1e3,tarDist(end)/1e3],[LB,LB],'-k')

axis([tarDist(1)/1e3, tarDist(end)/1e3, 0 3])

h1=xlabel('Distance (km)');
h2=ylabel('NEES');
title('NEES Comparison of Polar Conversions')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('Standard Conversion','Additive','Multiplicative', 'Modified Multiplicative','Cubature','Location','NorthEast') 

disp('Example 2: Spherical Measurement Conversion')

%We consider measurements at 45 degrees azimuth, 15 degrees elevation at
%ranges from 1 to 70km from the sensor.
numPoints=50;
numRuns=200;

%50 meter range standard deviation
sigmaR=50;
%10 degree angular standard deviation.
sigmaTheta=10*(pi/180);
sigmaEta=10*(pi/180);
SR=diag([sigmaR,sigmaTheta,sigmaEta]);
R=SR*SR';%Measurement variance

az=60*(pi/180);
el=30*(pi/180);

dirTarget=[cos(az)*cos(el);sin(az)*cos(el);sin(el)];
tarDist=linspace(1e3,70e3,numPoints);
locList=bsxfun(@times,dirTarget,tarDist);

numMethods=3;
RMSE=zeros(numPoints,numMethods);
NEES=zeros(numPoints,numMethods);

[xi,w]=fifthOrderCubPoints(3);

for curPoint=1:numPoints
    zCartTrue=locList(:,curPoint);
    z=Cart2Sphere(zCartTrue,0,true);
    
    MSECur=zeros(1,numMethods);
    NEESCur=zeros(1,numMethods);
    for curRun=1:numRuns 
        zMeas=z+SR*randn(3,1);

        [zCart,RCart]=spher2CartTaylor(zMeas,R,[],[],0);
        diff=zCartTrue-zCart;
        MSECur(1)=MSECur(1)+diff'*diff;
        NEESCur(1)=NEESCur(1)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=spher2CartTaylor(zMeas,R,[],[],1);
        diff=zCartTrue-zCart;
        MSECur(2)=MSECur(2)+diff'*diff;
        NEESCur(2)=NEESCur(2)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=spher2CartCubature(zMeas,SR,0,true,[],[],[],xi,w);
        diff=zCartTrue-zCart;
        MSECur(3)=MSECur(3)+diff'*diff;
        NEESCur(3)=NEESCur(3)+diff'*inv(RCart)*diff;
    end
    RMSE(curPoint,:)=sqrt(MSECur/numRuns);
    NEES(curPoint,:)=NEESCur/(numRuns*3);
end

figure(3)
clf
hold on
plot(tarDist/1e3,RMSE(:,1)./RMSE(:,3),'-r','linewidth',2)%Multiplicative
plot(tarDist/1e3,RMSE(:,2)./RMSE(:,3),'--b','linewidth',2)%Modified Multiplicative
plot(tarDist/1e3,RMSE(:,3)./RMSE(:,3),':g','linewidth',2)%Cubature

h1=xlabel('Distance (km)');
h2=ylabel('RMSE Ratio');
title('RMSE Ratio of Conversion to Cubature Conversion (Spherical Meas.)')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('Multiplicative', 'Modified Multiplicative','Cubature','Location','NorthEast') 

figure(4)
clf
hold on
plot(tarDist/1e3,NEES(:,1),'-r','linewidth',2)%Multiplicative
plot(tarDist/1e3,NEES(:,2),'--b','linewidth',2)%Modified Multiplicative
plot(tarDist/1e3,NEES(:,3),':g','linewidth',2)%Cubature

%The upper and lower bounds marking a 95% confidence region.
UB=ChiSquareD.invCDF(1-0.05/2,numRuns*3)/(numRuns*3);
LB=ChiSquareD.invCDF(0.05/2,numRuns*3)/(numRuns*3);

%Draw the bounds
plot([tarDist(1)/1e3,tarDist(end)/1e3],[UB,UB],'-k')
plot([tarDist(1)/1e3,tarDist(end)/1e3],[LB,LB],'-k')

axis([tarDist(1)/1e3, tarDist(end)/1e3, 0 3])

h1=xlabel('Distance (km)');
h2=ylabel('NEES');
title('NEES Comparison of Spherical Conversions')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('Multiplicative', 'Modified Multiplicative','Cubature','Location','NorthEast') 

disp('Example 3: Monostatic r-u-v Measurement Conversion')

%We consider measurements at 45 degrees azimuth, 15 degrees elevation at
%ranges from 1 to 70km from the sensor.
numPoints=50;
numRuns=200;

%10 meter range standard deviation (one-way)
sigmaR=10;
useHalfRange=true;
%Direction cosine standard deviations
sigmaU=5e-3;
sigmaV=5e-3;
SR=diag([sigmaR,sigmaU,sigmaV]);
R=SR*SR';%Measurement variance

%We offset the target in two angles from the look direction, which is the
%z-axis. This is not the same spherical coordinate system as used in the
%previous example.
az=10*(pi/180);
el=10*(pi/180);
dirTarget=zeros(3,1);
dirTarget(1)=sin(el)*sin(az);
dirTarget(2)=sin(el)*cos(az);
dirTarget(3)=sqrt(1-sum(dirTarget(1:2).^2));

tarDist=linspace(100e3,400e3,numPoints);
locList=bsxfun(@times,dirTarget,tarDist);

numMethods=3;
RMSE=zeros(numPoints,numMethods);
NEES=zeros(numPoints,numMethods);

[xi,w]=fifthOrderCubPoints(3);

for curPoint=1:numPoints
    zCartTrue=locList(:,curPoint);
    z=Cart2Ruv(zCartTrue,true);
    
    MSECur=zeros(1,numMethods);
    NEESCur=zeros(1,numMethods);
    for curRun=1:numRuns 
        zMeas=z+SR*randn(3,1);

        [zCart,RCart]=monostatRuv2CartTaylor(zMeas,R,useHalfRange,[],[],0);
        diff=zCartTrue-zCart;
        MSECur(1)=MSECur(1)+diff'*diff;
        NEESCur(1)=NEESCur(1)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=monostatRuv2CartTaylor(zMeas,R,useHalfRange,[],[],1);
        diff=zCartTrue-zCart;
        MSECur(2)=MSECur(2)+diff'*diff;
        NEESCur(2)=NEESCur(2)+diff'*inv(RCart)*diff;
        
        [zCart,RCart]=ruv2CartCubature(zMeas,SR,useHalfRange,[],[],[],xi,w);
        diff=zCartTrue-zCart;
        MSECur(3)=MSECur(3)+diff'*diff;
        NEESCur(3)=NEESCur(3)+diff'*inv(RCart)*diff;
    end
    RMSE(curPoint,:)=sqrt(MSECur/numRuns);
    NEES(curPoint,:)=NEESCur/(numRuns*3);
end

figure(5)
clf
hold on
plot(tarDist/1e3,RMSE(:,1)./RMSE(:,3),'-r','linewidth',2)%CM1 
plot(tarDist/1e3,RMSE(:,2)./RMSE(:,3),'--b','linewidth',2)%CM2
plot(tarDist/1e3,RMSE(:,3)./RMSE(:,3),':g','linewidth',2)%Cubature

h1=xlabel('Distance (km)');
h2=ylabel('RMSE Ratio');
title('RMSE Ratio of Conversion to Cubature Conversion (r-u-v Meas.)')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('CM1', 'CM2','Cubature','Location','NorthEast') 

figure(6)
clf
hold on
plot(tarDist/1e3,NEES(:,1),'-r','linewidth',2)%CM1
plot(tarDist/1e3,NEES(:,2),'--b','linewidth',2)%CM2
plot(tarDist/1e3,NEES(:,3),':g','linewidth',2)%Cubature

%The upper and lower bounds marking a 95% confidence region.
UB=ChiSquareD.invCDF(1-0.05/2,numRuns*3)/(numRuns*3);
LB=ChiSquareD.invCDF(0.05/2,numRuns*3)/(numRuns*3);

%Draw the bounds
plot([tarDist(1)/1e3,tarDist(end)/1e3],[UB,UB],'-k')
plot([tarDist(1)/1e3,tarDist(end)/1e3],[LB,LB],'-k')

axis([tarDist(1)/1e3, tarDist(end)/1e3, 0 2])

h1=xlabel('Distance (km)');
h2=ylabel('NEES');
title('NEES Comparison of r-u-v Conversions')
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
legend('CM1', 'CM2','Cubature','Location','NorthEast') 

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
