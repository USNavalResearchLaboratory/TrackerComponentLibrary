function [g,coeffs,mu]=TaylorTapering(nBar,sidelobedB,xyPoints,a)
%%TAYLORTAPERING The Taylor tapering is a set of amplitude weights for a
%          continuous circular (narrowband) aperture antenna that will
%          reduce the height of nBar of the close-in sidelobes at a cost of
%          widening the beam. Such a tapering can be discretized and
%          applied to the elements in a circular phased array (An array of
%          antenna elements can be viewed as a discrete approximation to a
%          continuous aperture). This function will provide the tapering
%          values at a set of discrete points given by xyPoints (the origin
%          is taken to be the center of the aperture). The radius of the
%          aperture can either be provided or is taken as value of the
%          farthest point provided. This function also returns coefficients
%          for one to efficiently compute the tapering at points on their
%          own. If xyPoints is empty, only the coefficients are returned.
%
%INPUTS: nBar The number of terms to use in the expansion. High values can
%             produce undesirable results.
%  sidelobedB The number of decibels of the ratio of the close-in sidelobe
%             voltages to the main lobe voltage. This must be a negative
%             number. A typical value is -30.
%    xyPoints A 2XN set of N points in the aperture plane at which the
%             tapering values should be evaluated. If this parameter is
%             omitted or an empty matrix is passed, then an empty matrix is
%             returned for the output g. The center of the aperture is
%             taken to be the origin. 
%           a The radius of the aperture. Tapering weights for points in
%             xyPoints outside of the aperture are taken to be zero. If
%             this parameter is omitted or an empty matrix is passed, then
%             the radius is taken to be the distance of the farthest point
%             from the origin in xyPoints.
%
%OUTPUTS: g The NX1 set of discretized Taylor tapering values evaluated at
%           the points given in xyPoints. If xyPoints is omitted, then this
%           is an empty matrix. All Taylor tapering values are positive and
%           real. The coefficients are not normalized.
% coeffs, mu These two outputs can be used to evaluate the Taylor tapering
%           values at arbitrary points. coeffs is an (nBar-1)X1 vector and
%           mu is an nBarX1 vector. Given a 2X1 point xy and a radius of
%           the aperture, set the normalized radius to p=pi*norm(xy)/a and
%           the Taylor tapering weight g at the point is
%           g=2/pi^2+sum(coeffs.*besselj(0,mu(1:(end-1))*p));
%        
%This implements the algorithm of [1]. Running
% xyPoints=[(0:20)/20;
%           zeros(1,21)];
% g=TaylorTapering(3,-25,xyPoints,1)
%provides the points that are in the first column of Table III in [2].
%
%Now, one must consider how large an aperture (in wavelengths to use to be
%able to get a desired beamwidth and widelobe level. In the following, the
%beamwidth is given in terms of "standard beamwidths". One standard
%beamwidth is lambda/(2*a) where lambda is the beamwidth and a is the
%radius of the aperture. As noted in [1], the beamwidth in radius between
%-3dB points for a uniformly illuminated circular aperture (no tapering) is
%1.029 standard beamwidths.
%
%The expression for the aperture size is given in terms of a value A which
%is given in terms of the sidelobe voltage ratio eta where
% eta=10^(sldelobedB/(-20));
%The value A is defined to be
% A=acosh(eta)/pi;
%The ideal space factor in radians (standard beamwidths) of a one-
%wavelength source is 
% beta0=(2/pi)*sqrt(acosh(eta)^2-acosh(eta/sqrt(2))^2);
%We also need the beam broadening factor sigma. This is
% sigma=mu(end)/sqrt(A^2+(nBar-1/2)^2);
%where mu is one of the outputs of the TaylorTapering function. mu can also
%be found as
% mu=BesselJZeros(1,nBar)/pi;
%Thus, the ratio of the radius a of the antenna to the wavelength lambda of
%the narrowband signal to be sent/received necessary to obtain a desired
%standard beamwidth beta in radians is
% a/lambda=beta0/(2*beta)*sigma
%Using the example from [2], to get a three degree beamwidth with -30db
%sidelobes with nBar=4, one has
% sidelobedB=-30;
% eta=10^(-sidelobedB/20)=31.6228
% A=acosh(eta)/pi=1.3200
% beta0=(2/pi)*sqrt(acosh(eta)^2-acosh(eta/sqrt(2))^2)=1.0569
% nBar=4;
% mu=BesselJZeros(1,nBar)/pi;
% sigma=mu(end)/sqrt(A^2+(nBar-1/2)^2)=1.1338
% betaVal=3*(pi/180);%3 degrees in radians
% a/lambda=beta0/(2*betaVal)*sigma=11.4427
%Thus an aperture radius of at least 11.4427 wavelengths is needed to
%achieve the desired performance. In practice, however, these numbers are
%only guidelines and one can often go lower than the designed level for a
%particular aperture.
%
%EXAMPLE 1:
%Here, we evaluate the tapering values for 30dB down on a fine grid of
%points to plot what the tapering gradient looks like.
% numPoints=300;
% points1D=linspace(-1,1,numPoints);
% [X,Y]=meshgrid(points1D,points1D);
% %The Taylor tapering weights, evaluated across the aperture. Points
% %outside the aperture are assigned a weight of 0. All Taylor weights are
% %real.
% xyPoints=[X(:)';Y(:)'];
% a=1;%Aperture radius=1.
% gTaylor=TaylorTapering(4,-30,xyPoints,a);
% gTaylor=reshape(gTaylor,numPoints,numPoints);
% 
% figure(1)
% clf
% surface(X,Y,gTaylor,'EdgeColor','None')
% colormap(jet(256))
% colorbar()
% view(45,45)
% light()
% axis square
% h1=xlabel('x');
% h2=ylabel('Taper Value');
% title('Taylor Tapering Weight')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here, we consider the array response when using tapering values for 30dB
%sidelobes on a circular array with lambda/2 spacing between elements.
% %First, we create a circular array. The element locations are given in
% %terms of the wavelength lambda, so lambda will not appear in the
% %equations for the sum beam.
% xyVals=getShaped2DLattice([25;25],'circular');
% 
% %Get the tapering.
% nBar=4;
% sldelobedB=-30;
% g=TaylorTapering(nBar,sldelobedB,xyVals);
% 
% T=diag(g);
% [Rsp,U,V]=standardUVBeamPattern(T,xyVals,'NormPowGain');
% 
% %Display the magnitude of the array response in decibels.
% figure(1)
% clf
% surface(U,V,10*log10(Rsp),'EdgeColor','None')
% colormap(jet(256));
% caxis([-70,0])
% colorbar()
% view(45,45)
% light()
% axis square
% h1=xlabel('u');
% h2=ylabel('v');
% title('Taylor Weighted Array Response (Decibels)')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] T. T. Taylor, "Design of circular apertures for narrow beamwidth and
%    low sidelobes," IRE Transactions on Antennas and Propagation, vol. 8,
%    no. 1, pp. 17-22, Jan. 1960.
%[2] R. C. Hansen, "Tables of Taylor distributions for circular aperture
%    antennas," IRE Transactions on Antennas and Propagation, vol. 8, no.
%    1, pp. 23-26, Jan. 1960.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Compute the first nBar-1 roots u of besselj(1,pi*u)
%These are defined by Equation 13 in [1]. 
mu=BesselJZeros(1,nBar)/pi;

%Equation 18 of [1] relates A to the sidelobe ratio. Here, it is just
%modified so that the side lobe level can be given in decibels.
eta=10^(sidelobedB/(-20));
A=acosh(eta)/pi;

%The square of the beam broadening factor in Equation 15 of [1].
sigma2=mu(nBar)^2/(A^2+(nBar-(1/2))^2);

b0Terms=besselj(0,pi*mu(1:(nBar-1)));

%Compute the space factor for m>0 at all of the mu points using Equation
%29 in [1].
F=zeros(nBar-1,1);
for m=1:(nBar-1)
    n=(1:(nBar-1));
    num=prod(1-mu(m)^2./(sigma2*(A^2+(n-(1/2)).^2)));
    denom=prod(1-mu(m)^2./mu([1:m-1,(m+1):(nBar-1)]).^2);
    
    F(m)=-b0Terms(m)*num/denom;
end

%We now build the coefficients for m=1:(nBar-1) in the sum in Equation 30
%of [1].
coeffs=zeros(nBar-1,1);
for m=1:(nBar-1)
    coeffs(m)=(2/pi^2)*F(m)/(b0Terms(m)^2);
end

%If discretized tapering values are desired.
if(nargin>2&&~isempty(xyPoints))
    numPoints=size(xyPoints,2);
    
    if(nargin<4||isempty(a))
        %The maximum distance from the origin to a point is taken to be the
        %radius of the aperture.
        a2=max(sum(xyPoints.^2,1));
        a=sqrt(a2);
    end
    
    g=zeros(numPoints,1);
    for curPoint=1:numPoints
        rho2=sum(xyPoints(:,curPoint).^2,1);
        
        if(rho2<=a2)
            rho=sqrt(rho2);
            %The normalized radius at this point.
            p=pi*rho/a;

            %Terms for m>=0
            g(curPoint)=sum(coeffs.*besselj(0,mu(1:(end-1))*p));
            %Add the m=0 term.
            g(curPoint)=g(curPoint)+2/pi^2;
        end
    end
else%If no tapering values are requested, return an empty matrix.
    g=[];
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
