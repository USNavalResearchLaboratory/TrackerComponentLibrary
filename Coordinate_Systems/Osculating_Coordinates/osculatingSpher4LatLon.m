function [r,spherCent]=osculatingSpher4LatLon(latLon,a,f,is2D,useGeoid,geoidOptions)
%%OSCULATINGSPHER4LATLON Given the latitude and longitude of a point on a
% reference ellipsoid, determine the radius and center location of an
% osculating sphere going through that point. If using the WGS-84 reference
% ellipsoid, one can optionally shift the osculating sphere so that the the
% sphere goes through the same point projected onto the geoid (a
% theoretical surface of constant gravitational potential that is used to
% define mean sea level) instead of a point on the ellipsoid. The
% osculating sphere is such that the azimuth equals the longitude and the
% elevation equals the latitude of the reference point on the ellipsoid.
% Additionally, the local tangent planes at the point (on the ellipsoid and
% the sphere) are the same. The radius of curvature of the sphere at that
% point should nominally equal the radius of curvature of the ellipsoid at
% that point. However, since the ellipsoid has two radii of curvature, the
% geometric mean is used (the mean radius of curvature). On the other hand,
% if is2D is true, then just the radius of curvature in the meridian will
% be used, which is what perfectly matches a cut of the ellipsoid through
% the origin and pole (A 2D ellipse in the x-z plane).
%
%INPUTS: latLon The 2XN set of points of the format [latitude;longitude] in
%               radians where the osculating spheres should touch the
%               ellipsoid.
%             a The semi-major axis length of the ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the ellipsoid. If this argument is
%               omitted or an empty matrix is passed, the value in
%               Constants.WGS84Flattening is used.
%          is2D Indicates whether one just cares about a 2D cut of the
%               ellipsoid. If so, then the radius used will be the radius
%               of curvature in the meridian. The default if omitted or an
%               empty matrix is passed is false.
%      useGeoid If this is true, then the spherCent output is offset so
%               that a point at latLon with a height that equals the height
%               of the geoid at that point will have zero height with
%               respect to the osculating sphere. Alternatively, if this is
%               false, then the osculating sphere will touch the reference
%               ellipse at latLon. If this is true, and one doesn't pass
%               the geoid heights in GeoidOptions, then a and f must
%               either be empty matrices or equal to
%               Constants.WGS84SemiMajorAxis and Constants.WGS84Flattening.
%               The default if omitted or an empty matrix is passed is
%               false.
%  geoidOptions If the geoid heights at the points in latLon are known and
%               useGeoid is true, this is a length N vector of those
%               heights. Otherwise, this is an optional structure that is
%               only used if useGeoid is true and that selects what geoid
%               is used. Possible members are tideSys, useNGAApprox,
%               modelType, and coeffData, which correspond to the same-
%               named inputs to getEGMGeoidHeight. The mean-tide geoid is
%               used. The default if omitted or an empty matrix is passed
%               are the defaults of getEGMGeoidHeight except here,
%               tideSys=1, which selects the mean-tide geoid.
%
%OUTPUTS: r The 1XN set of radiii, one for each osculating sphere.
% spherCent The 3XN locations of the centers of the osculating spheres with
%           respect to the center of the reference ellipsoid.
%
%The formulae in Section 3 of [1] are used.
%
%Note that if the osculating sphere shifted to tough the geoid is desired,
%this function will be very slow if one hasn't called CompileCLibraries to
%compile the spherical harmonic synthesis functions.
%
%EXAMPLE 1:
%Here, we consider a 2D slice of an ellipsoid that goes through
%longitude=0, we choose a point. We then find the osculating sphere. The
%ellipse cut of the ellipsoid and the circle cut of the sphere are then
%plotted with the point. We ALSO plot the circle cut obtained with
%is2D=true. Since we are only considering a 2D cut, one will see that that
%circle is actually better aligned with the curvature of the ellipse at
%that point. However, if one were to take a 2D cut in other directions, it
%would appear to be a worse fit.
% a=20;%Semi-major axis.
% f=0.5;%Flattening factor.
% b=a*(1-f);%The semi-minor axis of the ellipsoid.
% 
% figure(1)
% clf
% hold on
% %Draw a 2D cut of the ellipsoid.
% A=inv([a^2,0;
%        0,b^2]);
% drawEllipse([0;0],A,1,'linewidth',4)
% latLon=[50;0]*(pi/180);
% oscPt=ellips2Cart([latLon;0],a,f);
% [r,spherCent]=osculatingSpher4LatLon(latLon,a,f,false);
% [r2D,spherCent2D]=osculatingSpher4LatLon(latLon,a,f,true);%is2D=true
% %Draw the both 2D cuts of the sphere:
% A=inv([r^2,0;
%        0,r^2]);
% drawEllipse([spherCent(1);spherCent(3)],A,1,'-.','linewidth',4)
% A=inv([r2D^2,0;
%        0,r2D^2]);
% drawEllipse([spherCent2D(1);spherCent2D(3)],A,1,'--','linewidth',2)
% 
% %Note that spherCent(2)==0.
% scatter(oscPt(1),oscPt(3),100,'filled')
% legend('Ellipsoid Cut','Osculating Sphere','Osculating Sphere 2D','location','southwest')
% h1=xlabel('x');
% h2=ylabel('y');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%In this example, we demonstrate that a point placed on the geoid will have
%a zero (within finite precision limits) height when placed on the
%osculating sphere given the useGeoid option.
% latLonPt=deg2rad([19.7216;-155.0849]);%Hilo, Hawaii.
% tideSys=1;
% useNGAApprox=false;
% modelType=0;
% geoidHeight=getEGMGeoidHeight(latLonPt,tideSys,useNGAApprox,modelType);
% is2D=false;
% useGeoid=true;
% geoidOptions.tideSys=tideSys;
% geoidOptions.modeTYpe=modelType;
% geoidOptions.useNGAApprox=useNGAApprox;
% [r,spherCent]=osculatingSpher4LatLon(latLonPt,[],[],is2D,useGeoid,geoidOptions);
% 
% llhOnGeoid=[latLonPt;geoidHeight];
% latLonHOsc=ellips2SpecOscCoords(llhOnGeoid,r,spherCent);
% %The height of a point on the geoid when converted to the osculating sphere
% %should be 0.
% heightOnSphere=latLonHOsc(3)
%
%REFERENCES:
%[1] P. Williams and D. Last, "On Loran-C time-difference to co-ordinate
%    converters," in Proceedings of the 32nd Annual Convention & Technical
%    Symposium of the International Loran Association, Boulder, CO, 3-7
%    Nov. 2003.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(useGeoid))
    useGeoid=false;
end

if(nargin<4||isempty(is2D))
    is2D=false;
end

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<2||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%The squared eccentricity of the ellipsoid.
e2=f*(2-f);

lat=latLon(1,:);%Latitude.
lon=latLon(2,:);

Phi=lat;

denomTerm=sqrt(1-e2*sin(Phi).^2);
%Equation 3.17 in [1]. The radius of curvature in the meridian.
M0=a*(1-e2)./denomTerm.^3;

if(is2D)
    r=M0;
else
    %Equation 3.18 in [1]. The radius of curvature in the prime vertical.
    N0=a./denomTerm;

    %Equation 3.16 in [1].
    r=sqrt(M0.*N0);
end

numPts=size(latLon,2);
xCartEllips=ellips2Cart([latLon(1:2,:);zeros(1,numPts)],a,f);
xCartSpher=spher2Cart([r;lon;lat]);

spherCent=-(xCartSpher-xCartEllips);

if(useGeoid)
    if(nargin>5&&ismatrix(geoidOptions))
        %If the geoid heights are explicitly given.
        geoidHeight=geoidOptions(:).';
    else
        if(a~=Constants.WGS84SemiMajorAxis||f~=Constants.WGS84Flattening)
            error('The geoid option is only available when using the WGS84 reference ellipsoid.')
        end
        %Default geoid options.
        useNGAApprox=false;
        modelType=[];%EGM2008 Model
        tideSys=1;%Mean tide.
        coeffData=[];
        if(nargin>5&&~isempty(geoidOptions))
            if(isfield(geoidOptions,'useNGAApprox'))
                useNGAApprox=geoidOptions.useNGAApprox;
            end
            if(isfield(geoidOptions,'modelType'))
                modelType=geoidOptions.modelType;
            end
            if(isfield(geoidOptions,'coeffData'))
                coeffData=geoidOptions.coeffData;
            end
        end
    
        geoidHeight=getEGMGeoidHeight(latLon(1:2,:),tideSys,useNGAApprox,modelType,coeffData).';
    end
    %Get the "Up" direction by which to move the center of the sphere so
    %that it touches the geoid (not the reference ellipsoid).
    uUp=zeros(3,numPts);
    for curPt=1:numPts
        uUp(:,curPt)=getENUAxes(latLon(1:2,curPt),true,a,f);
    end
    spherCent=spherCent+bsxfun(@times,geoidHeight,uUp);
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
