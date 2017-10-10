function g=TaylorLinearTapering(nBar,sidelobedB,xPoints,a)
%%TAYLORLINEARTAPERING The Taylor tapering is a set of amplitude weights
%          for a continuous linear (narrowband) aperture antenna that will
%          reduce the height of nBar of the close-in sidelobes at a cost of
%          widening the beam. Such a tapering can be discretized and
%          applied to the elements in a linear phased array (An array of
%          antenna elements can be viewed as a discrete approximation to a
%          continuous aperture). This function will provide the tapering
%          values at a set of discrete points given by xPoints (the origin
%          is taken to be the center of the aperture). The radius of the
%          aperture can either be provided or is taken as value of the
%          farthest point provided.
%
%INPUTS: nBar The number of terms to use in the expansion. High values can
%             produce undesirable results.
%  sidelobedB The number of decibels of the ratio of the close-in sidelobe
%             voltages to the main lobe voltage. This must be a negative
%             number. A typical value is -30.
%    xyPoints An NX1 or 1XN set of N points at which the tapering values
%             should be evaluated. The center of the aperture is taken to
%             be the origin. 
%           a The radius of the aperture. Tapering weights for points in
%             xPoints outside of the aperture are taken to be zero. If
%             this parameter is omitted or an empty matrix is passed, then
%             the radius is taken to be the distance of the farthest point
%             from the origin in xPoints.
%
%OUTPUTS: g The NX1 set of discretized Taylor tapering values evaluated at
%           the points given in xPoints. All Taylor tapering values are
%           positive and real. The coefficients are not normalized.
%
%This function implements the algorithm in [1].
%
%EXAMPLE 1:
%Here, we evaluate the tapering values for 30dB down on a large array of
%points to see what the tapering looks like.
% nBar=4;
% sidelobedB=-30;
% numPoints=100;
% %Generate points symmetric about the origin. This makes the aperture size
% %1.
% k=0:(numPoints-1);
% xPoints=(k-(1/2)*numPoints+1/2)/numPoints;
% g=TaylorLinearTapering(nBar,sidelobedB,xPoints);
% 
% figure(1)
% clf
% plot(xPoints,g)
% h1=xlabel('x');
% h2=ylabel('y');
% title('Taylor Tapering Weight')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here, we consider the array response when using tapering values for 30dB
%sidelobes on a linear array with lambda/2 spacing between elements.
% nBar=4;
% sidelobedB=-30;
% numPoints=30;
% %Generate points symmetric about the origin.
% k=0:(numPoints-1);
% xPoints=(k-(1/2)*numPoints+1/2)/2;
% g=TaylorLinearTapering(nBar,sidelobedB,xPoints);
% 
% %Tapering matrix
% T=diag(g);
% 
% %Now, display the response with the tapering
% [Rsp,U]=standardUVBeamPattern(T,xPoints,'NormPowGain');
% 
% figure(2)
% clf
% plot(U,10*log10(Rsp),'-b','LineWidth',2);
% axis([-1, 1, -40 1])
% h1=xlabel('u');
% h2=ylabel('Amplitude Response (decibels)');
% title('Taylor Weighted Array Response')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] T. T. Taylor, "Design of line-source antennas for narrow beamwidth and
%    low side lobes," IRE Transactions on Antennas and Propagation, vol. 3,
%    no. 1, pp. 16-28, Jan. 1955.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%After Equation 35 in [1], it is noted that cosh(pi*A) is the side-lobe
%ratio. Thus, we get the following expression for A in terms of the
%sidelobe level in decibels:
eta=10^(sidelobedB/(-20));
A=acosh(eta)/pi;

%The square of Sigma in Equation 41 in [1].
sigma2=nBar^2/(A^2+(nBar-1/2)^2);

%The suggestion for C after Equation 48 in [1].
C=cosh(pi*A);

%Fm holds values of F in Equation 48 for z being an integer>=1.
Fm=zeros(nBar-1,1);
for m=1:(nBar-1)
%Equation 48 contains singularities when z=m is an integer. To evaluate it
%when z is an integer, we first simplify Q in Equation 47 to eliminate the
%singularity.

    %Equation 47 in [1], taking the limit of z to an integer to eliminate
    %the singularity.
    np=[1:(m-1),(m+1):(nBar-1)];
    Qm=(-1)^(m+1)/(2*prod(1-(m./np).^2));

    n=1:(nBar-1);
    Fm(m)=C*Qm*prod(1-m^2./(sigma2*(A^2+(n-1/2).^2)));
end
F0=C;%The case where z=0.

numPoints=length(xPoints);

if(nargin<4||isempty(a))
    %The maximum distance from the origin to a point is taken to be the
    %radius of the aperture.
    a=sqrt(max(sum(xPoints.^2,1)));
end

%The loop below implements Equation 63 in [1].
g=zeros(numPoints,1);
g(:)=F0/(2*pi);
for curPoint=1:numPoints
    rho=norm(xPoints(curPoint));

    if(rho<=a)
        %The normalized radius at this point.
        p=pi*rho/a;

        m=(1:(nBar-1))';
        g(curPoint)=g(curPoint)+(1/pi)*sum(Fm.*cos(m*p));
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
