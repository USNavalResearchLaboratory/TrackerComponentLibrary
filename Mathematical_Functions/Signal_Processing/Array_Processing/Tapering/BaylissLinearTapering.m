function g=BaylissLinearTapering(sidelobedB,N,xPoints,a)
%%BAYLISSLINEARTAPERING The Bayliss tapering is a set of complex amplitude
%          weights for a continuous linear (narrowband) aperture antenna
%          that will form a difference beam (odd symmetry about an axis)
%          and hold the four closest sidelobes to a desired level. Such a
%          tapering can be discretized and applied to the elements in a
%          circular phased array (An array of antenna elements can be
%          viewed as a discrete approximation to a continuous aperture).
%          This function will provide the tapering values at a set of
%          discrete points given by xPoints (the origin is taken to be the
%          center of the aperture). The radius of the aperture can either
%          be provided or is taken as value of the farthest point provided.
%
%INPUTS: sidelobedB The number of decibels of the ratio of the close-in
%             sidelobe voltages to the main lobe voltage. This must be a 
%             negative number. A typical value is -30.
%           N The Bayliss tapering is computed using a certain number of
%             terms. Using too many terms can be undesirable as edge
%             illumination increases, as noted in [1]. If this parameter is
%             omitted or an empty matrix is passed, then the default of 17
%             is used. In [1], it is suggested that N be chosen to be
%             <2*a/lambda, where a is the radius of the aperture and
%             lambda the wavelength.
%    xyPoints An NX1 or 1XN set of N points at which the tapering values
%             should be evaluated. The center of the aperture is taken to
%             be the origin. 
%           a The radius of the aperture. Tapering weights for points in
%             xPoints outside of the aperture are taken to be zero. If
%             this parameter is omitted or an empty matrix is passed, then
%             the radius is taken to be the distance of the farthest point
%             from the origin in xPoints.
%
%OUTPUTS: g The NX1 set of discretized Bayliss tapering values evaluated at
%           the points given in xPoints. All Bayliss tapering values are
%           imaginary. The coefficients are not normalized.
%
%This function implements the algorithm given in the appendix of in [1]
%using the polynomial interpolation values in the table below Figure 4.
%This approximation means that low sidelobe patterns (-45 dB and below)
%will not have good fidelity sidelobes.
%
%EXAMPLE 1:
%Here, we evaluate the tapering values for 30dB down on a fine grid of
%points to plot what the amplitude and phase of the tapering weights looks
%like.
% numPoints=300;
% xPoints=linspace(-1,1,numPoints);
% %The Bayliss tapering weights, evaluated across the aperture.
% a=1;%Aperture radius=1.
% %gTest=BaylissLinear();
% gBayliss=BaylissLinearTapering(-30,17,xPoints,a);
% 
% figure(1)
% clf
% plot(xPoints,abs(gBayliss),'-b','linewidth',2)
% h1=xlabel('x');
% h2=ylabel('Imaginary Weight');
% title('Bayliss Tapering Amplitude')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% figure(2)
% clf
% plot(xPoints,angle(gBayliss),'-b','linewidth',2)
% h1=xlabel('x');
% h2=ylabel('Imaginary Weight');
% title('Bayliss Tapering Phase')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here, we consider the array response when using tapering values for 30dB
%sidelobes on a linear array with lambda/2 spacing between elements.
% N=17;
% sidelobedB=-30;
% Nx=41;%There are 2*Nx+1 points total.
% %Generate points symmetric about the origin.
% xPoints=(-(Nx-1)/2):1/2:((Nx-1)/2);
% g=BaylissLinearTapering(sidelobedB,N,xPoints);
% 
% %Now, display the response with the tapering
% T=diag(g);
% [Rsp,U]=standardUVBeamPattern(T,xPoints,'NormRealVal');
% 
% figure(2)
% clf
% plot(U,Rsp,'-b','LineWidth',2);
% h1=xlabel('u');
% h2=ylabel('Array Response');
% title('Bayliss Weighted Array Response')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] E. T. Bayliss, "Design of monopulse antenna difference patterns with
%    low sidelobes," The Bell System Technical Journal, vol. 47, no. 5, pp.
%    623-650, May-Jun. 1968.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(N))
   N=17; 
end

%Defined below Equation 41 in [1].
mu=((0:(N-1))+1/2)';

%This holds the coefficients for the interpolating polynomials given below
%Figure 4 in [1]. The polynomials all take the desired sidelobe level in
%decibels as an input parameter. The first row is for the term A, which
%is a translation of the SNR parameter to a parameter in the paper. The
%next four rows are xi_1 to xi_4, which are the locations of the first four
%zeroes in the modified pattern. The final row is for p_0, which is related
%to the point at which the peak of the asymptotic difference pattern=1.
polyCoeffTable=[0.30387530,-0.05042922,-0.00027989,-0.00000343,-0.00000002;
                0.98583020,-0.03338850, 0.00014064, 0.00000190, 0.00000001;
                2.00337487,-0.01141548, 0.00041590, 0.00000373, 0.00000001;
                3.00636321,-0.00683394, 0.00029281, 0.00000161, 0;
                4.00518423,-0.00501795, 0.00021735, 0.00000088, 0;
                0.47972120,-0.01456692,-0.00018739,-0.00000218,-0.00000001];

A  =polyCoeffTable(1,1)+sidelobedB*(polyCoeffTable(1,2)+sidelobedB*(polyCoeffTable(1,3)+sidelobedB*(polyCoeffTable(1,4)+sidelobedB*polyCoeffTable(1,5))));
xi1=polyCoeffTable(2,1)+sidelobedB*(polyCoeffTable(2,2)+sidelobedB*(polyCoeffTable(2,3)+sidelobedB*(polyCoeffTable(2,4)+sidelobedB*polyCoeffTable(2,5))));
xi2=polyCoeffTable(3,1)+sidelobedB*(polyCoeffTable(3,2)+sidelobedB*(polyCoeffTable(3,3)+sidelobedB*(polyCoeffTable(3,4)+sidelobedB*polyCoeffTable(3,5))));
xi3=polyCoeffTable(4,1)+sidelobedB*(polyCoeffTable(4,2)+sidelobedB*(polyCoeffTable(4,3)+sidelobedB*(polyCoeffTable(4,4)+sidelobedB*polyCoeffTable(4,5))));
xi4=polyCoeffTable(5,1)+sidelobedB*(polyCoeffTable(5,2)+sidelobedB*(polyCoeffTable(5,3)+sidelobedB*(polyCoeffTable(5,4)+sidelobedB*polyCoeffTable(5,5))));
%p0 =polyCoeffTable(6,1)+sldelobedB*(polyCoeffTable(6,2)+sldelobedB*(polyCoeffTable(6,3)+sldelobedB*(polyCoeffTable(6,4)+sldelobedB*polyCoeffTable(6,5))));

Z=zeros(N+1,1);
%Equation 15
Z(1)=0;%The Z(0) term
%Now, the moved zeros
Z(2)=xi1;
Z(3)=xi2;
Z(4)=xi3;
Z(5)=xi4;
for k=5:N
    %The location of the non-moved zeros as given by Equation 13.
    Z(k+1)=sqrt(A^2+k^2);
end

sigma=(N+1/2)/Z(N+1);

B=zeros(N,1);
for m=0:(N-1)
    %This loop implements Equation 47.
    n=1:(N-1);
    num=prod(1-((m+1/2)./(sigma*Z(1+n))).^2);
    l=[0:(m-1),(m+1):(N-1)];
    denom=prod(1-((m+1/2)./(l+1/2)).^2);
    %We just use C=1.
    B(m+1)=1/(2*1j)*(-1)^m*(m-1/2)^2*num/denom;
end

numPoints=length(xPoints);

if(nargin<4||isempty(a))
    %The maximum distance from the origin to a point is taken to be the
    %radius of the aperture.
    a=sqrt(max(sum(xPoints.^2,1)));
end

%The loop below implements Equation 41 in [1].
g=zeros(numPoints,1);
for curPoint=1:numPoints
    if(abs(xPoints(curPoint))<=a)
        %The normalized radius at this point.
        p=pi*xPoints(curPoint)/a;
        g(curPoint)=g(curPoint)+sum(B.*sin(mu*p));
    end
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
