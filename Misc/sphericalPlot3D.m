function h=sphericalPlot3D(f,azRange,elRange,numPoints,systemType,varargin)
%%SPHERICALPLOT3D Given a function that provides a range as a function of
%                 spherical coordinates, create a three dimensional surface
%                 plot of the function over a particular range.
%
%INPUTS: f A handle to a function that takes points in spherical
%          coordinates. The function is called as f(AzEl), where AzEl(1) is
%          azimuth and AzEl(2) is elevation. The function returns a one-way
%          range. The range should be >=0.
%  azRange, elRange These are 2X1 vectors with bounds in azimuth and
%          elevation. The first entry in each is the lower bound. The
%          defaults for omitted inputs or empty matrices passed are
%          respectively azRange=[-pi;pi]; and elRange=[-pi/2;pi/2];
% numPoints A 2X1 vector containing how many points to use respectively in
%           azimuth and in elevation. The default if omitted or an empty
%           matrix is passed is [50;25];
% systemType An optional parameter specifying the axes from which
%           the angles are measured. Possible values are
%           0 (The default if omitted) Azimuth is measured counterclockwise
%             from the x-axis in the x-y plane. Elevation is measured up
%             from the x-y plane (towards the z-axis). This is consistent
%             with common spherical coordinate systems for specifying
%             longitude (azimuth) and geocentric latitude (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one desires the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%  varargin Sets of values that should be passed to the plot function to
%           format the ellipses or that will be passed to the surf function
%           to format the surface. For example, one often wants to pass
%           'EdgeColor','None' so that edges are not drawn.
%
%OUTPUTS: h The handle to the chart surface graphic that was produced by
%           this function.
%
%EXAMPLE 1:
%Here we draw what might be an antenna response.
% f=@(AzEl)(1+2*cos(2*(pi/2-AzEl(2))));
% figure(1)
% clf
% sphericalPlot3D(f,[],[],[200;100],[],'EdgeColor','None');
% axis([-3 3 -3 3 -3 3])
% light();
% view(45,20)
% h1=xlabel('x');
% h2=ylabel('y');
% h3=zlabel('z');
% title('Plot of 1+2 cos(2 (\pi/2-elevation))')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%Here we draw a sphere.
% f=@(AzEl)(1);
% figure(2)
% clf
% sphericalPlot3D(f);
% axis square
% h1=xlabel('x');
% h2=ylabel('y');
% h3=zlabel('z');
% title('Plot of a Unit Sphere')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Defaults
if(nargin<2||isempty(azRange))
    azRange=[-pi;pi];
end
if(nargin<3||isempty(elRange))
    elRange=[-pi/2;pi/2];
end
if(nargin<4||isempty(numPoints))
    numPoints=[50;25];
end
if(nargin<5||isempty(systemType))
    systemType=0;
end

az=linspace(azRange(1),azRange(2),numPoints(1));
el=linspace(elRange(1),elRange(2),numPoints(2));
[Az,El]=meshgrid(az,el);

R=zeros(numPoints(1),numPoints(2));
for curPoint=1:prod(numPoints)
    AzEl=[Az(curPoint);El(curPoint)];
    R(curPoint)=f(AzEl);
end
points=[R(:)';Az(:)';El(:)'];

useHalfRange=true;
points=spher2Cart(points,systemType,useHalfRange);

X=reshape(points(1,:),numPoints(1),numPoints(2));
Y=reshape(points(2,:),numPoints(1),numPoints(2));
Z=reshape(points(3,:),numPoints(1),numPoints(2));

h=surf(X,Y,Z,varargin{:});

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
