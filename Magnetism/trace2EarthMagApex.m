function [apexPoints,signVals,exitCode]=trace2EarthMagApex(zCart,modSelParam,a,f,RelTol,AbsTol,maxSteps)
%TRACE2EARTHMAGAPEX Given a point in the in the International Terrestial
%           Reference System (ITRS), a type of Earth-Centered Earth-Fixed
%           (ECEF) system trace along the magnetic field line of the Earth
%           until the apex is found. The apex is the highest point of the
%           magnetic field line above the reference ellipsoid. Such tracing
%           is used for determining apex coordinates and quasi-dipole
%           coordinates, both of which play a role in ionospheric analyses.
%
%INPUTS: zCart One or more points given in Cartesian coordinates in the
%              ITRS with units of meters. zCart is a 3XN matrix with each
%              column having the format [x;y;z]. One would generally 
%              assume that points are not underground.
%  modSelParam An optional parameter selecting or providing the magnetic
%              field model to use. If omitted or an empty matrix is
%              passed, the International Geomagnetic Reference Field
%              (IGRF) at the latest epoch of the model is used via the
%              function getIGRFCoeffs. Possible other values are
%              1) year The parameters of the IGRF for the specified epoch 
%                 year are used. The year is in the Gregorian calendar and
%                 is specified as noted in the comments to the
%                 getIGRFCoeffs  function.
%              2) modSelParam is a structure where modSelParam.algorithm
%                 selects the algorithm to use Possible values for
%                 modSelParam.algorithm are 'IGRF', 'WMM' and 'preloaded'.
%                 If modSelParam.algorithm is 'IGRF or 'WMM', then the IGRF
%                 or the World Magnetic Model (WMM) is used at the
%                 fractional year specified by modSelParam.year (using the
%                 function getIGRFCoeffs or getWMMCoeffs). If
%                 modSelParam.algorithm is preloaded, then the parameters
%                 for the model are given by the modSelParam.C,
%                 modSelParam.S, modSelParam.a, and modSelParam.c, where
%                 the parameters have the same format as the fully
%                 normalized outputs of getIGRFCoeffs or getWMMCoeffs. See
%                 comments below for the use of custom coeffients.
%            a The semi-major axis of the reference ellipsoid. The apex is
%              defined in terms of a maximum magnetic field line height
%              above the reference ellipsoid. If this argument is omitted
%              or an empty matrix is passed, the value in
%              Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%       RelTol The maximum relative error tolerance allowed for the
%              Runge-Kutta algorithm to trace the path towards the apex.
%              If omitted or an empty matrix is passed, the default value
%              of 1e-6 is used.
%       AbsTol The absolute error tolerance allowed, a positive scalar.
%              If omitted or an empty matrix is passed, the default value
%              of 1e-9 is used. This is used in the Runge-Kutta algorithm
%              to trace along the path towards the apex, and it is used
%              for the fminbnd function to find the apex location once it
%              has been bounded.
%     maxSteps The maximum allowable number of steps to perform the
%              adaptive Runge-Kutta integration along a magnetic field
%              line to find the apex. If omitted or an empty matrix is
%              passed, the default of 1024 is used.
%
%OUTPUS: apexPoints A 3XN matrix of the apex points corresponding to the 
%               points in zCart. The apex points are obtained by tracing the
%               magnetic field line starting from zCart until reaching the
%               point farthest from the reference ellipsoid. Near the
%               geomagnetic poles, these points will be very far from the
%               Earth. Thus, if the point ends up being 1e37 m or more away
%               from the Earth, then integration is stopped and the
%               components of apexPoints are just set to Inf (rather than
%               returning an error). If  exitCode~=0, for any point, an
%               empty matrix is returned.
%      signVals An NX1 vector of values (+1, or -1) indicating whether the
%               direction to the apex point along the field line was
%               obtained by tracing in the direction of the magnetic field
%               or opposite to it. This can be used to determine which
%               magnetic hemisphere the points in zCart are in, because
%               +1 means that one starts closer to the magnetic North pole
%               (geographic South pole), and -1 means that one starts
%               closer to the magnetic South pole (geographic North pole).
%               If an error occurs, an empty matrix is returned.
%      exitCode A NX1 vector of codes indicating whether an error occurred.
%               If an error occurs, then the algorithm is halted at the
%               point causing the error.
%               0: Integration was successful.
%               1: Unable to get a small enough step size.
%               2: Apex not found within the maximum number of iterations.
%               3: Non-finite number encountered.
%               4: Problems performing a line search to the peak point.
%
%Basically, starting at the given point, one traces along a magnetic field
%line until reaching the highest point above the reference ellipsoid as
%defined in [2]. The tracing is done using an order 5(4) adaptive
%Runge-Kutta method until the direction of the field line is descending, in
%which case the peak is bounded between the last two points on the field
%line. An interpolation routine is then used with fminbnd to find the peak.
%This assumes that there is only one peak, which is safe to assume with the
%IMM and WMM. On the other hand, if using a custom magnetic field model,
%there might be more peaks and this function will not return the correct
%solution.
%
%Note that in [1], the definition of the apex is given in terms of a
%maximum height or a magnetic field line above the geoid, not above the
%reference ellipsoid (which is an approximation of the geoid). The
%definition with respect to the geoid seems to not be used in the
%literature. Moreover, the fast undulations of the geoid would make a
%magnetic field line have more than one peak, requiring that one find every
%peak to determine where the apex is.
%
%REFERENCES:
%[1] T. E. VanZandt, W. L. Clark, and J. M. Warnock, "Magnetic apex
%    coordinates: A magnetic coordinate system for the ionospheric f2
%    layer," Journal of Geophysical Research, vol. 77, no. 13, pp. 2406-
%    2411, 1 May 1972.
%[2] A. D. Richmond, "Ionospheric electrodynamics using magnetic apex
%    coordinates," Journal of Geomagnetism and Geoelectricity, vol. 47,
%    no. 2, pp. 191-212, 1995.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(modSelParam))
    %If no magnetic model is specified, use the IGRF at the reference
    %epoch.
    [C,S,aMagMod,cMagMod]=getIGRFCoeffs([],true);
elseif(isa(modSelParam,'double'))
    %If just a scalar is passed, assume that it is the epoch year for the
    %IGRF.
    [C,S,aMagMod,cMagMod]=getIGRFCoeffs(modSelParam,true);
else
    %Otherwise, a structrue should be passed indicating the model and the
    %year.
    switch(modSelParam.algorithm)
        case 'WMM'%World Magnetic Model
            [C,S,aMagMod,cMagMod]=getWMMCoeffs(modSelParam.year,true);
        case 'IGRF'%Intergnational Geomagnetic Reference Field
            [C,S,aMagMod,cMagMod]=getIGRFCoeffs(modSelParam.year,true);
        case 'preloaded'%Custom coefficients.
            C=modSelParam.C;
            S=modSelParam.S;
            aMagMod=modSelParam.a;
            cMagMod=modSelParam.c;
        otherwise
            error('Unknown model specified')
    end
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(RelTol))
    RelTol=1e-6;
end

if(nargin<6||isempty(AbsTol))
    AbsTol=1e-9;
end

if(nargin<7||isempty(maxSteps))
    maxSteps=1024;
end

%The parameters selecting the Runge-Kutta algorithm to use for the
%numerical integration.
order=5;
solutionChoice=0;

numPoints=size(zCart,2);
apexPoints=zeros(3,numPoints);
signVals=zeros(numPoints,1);
exitCode=zeros(numPoints,1);
for curPoint=1:numPoints
    zITRS=zCart(:,curPoint);

    %The initial stepsize is arbitrarily set to 10km.
    deltaS=10e3;
    %The maximum stepsize is arbitrarily set very high, so that it can be
    %quickly detected when the apex is extremely far away, because the
    %chosen point is very near the magnetic pole.
    deltaSMaxMag=1e36;

    xCur=zITRS;%The initial Cartesian location.
    sCur=0;%Distance traveled from the initial point is zero.

    %First, we have to decide which direction to go. We want to go in the
    %direction of increasing ellipsoidal height.
    dxdsCur=dxds(xCur,1,C,S,aMagMod,cMagMod);%The positive direction.

    uVert=getEllipsVert(zITRS,a,f);

    %If the angle between the direction and the local vertical is >90
    %degrees, then it is a decreasing direction. If it is greater, then it
    %is an increasing direction. This is the same as checking whether the
    %dot product is negative (the dot product is proportional to the cosine
    %of the angle).
    if(dot(uVert,dxdsCur)<0)
        signVal=-1;%Go the other way.
    else
        signVal=1;
    end
    signVals(curPoint)=signVal;

    %The path direction function.
    derivFun=@(x,s)dxds(x,signVal,C,S,aMagMod,cMagMod);

    %Now, we integrate along the path until the next point takes us to a
    %location where the direction along the magnetic field line is no longer
    %going in a direction of increasing ellipsoidal height. That means that we
    %have passed the apex. Thus, the apex height can be interpolated beteen
    %those final two points.
    apexPoint=[];
    for curStep=1:maxSteps
        xPrev=xCur;
        sPrev=sCur;

        %Arbitrary minimum step size.
        deltaSMinMag=2^4*eps(sCur);

        [deltaS,xCur,sCur,k,dxdsCur,exitCode(curPoint)]=performOneAdaptiveRKStep(xCur,sCur,derivFun,deltaS,deltaSMinMag,deltaSMaxMag,dxdsCur,order,solutionChoice,AbsTol,RelTol);

        if(exitCode(curPoint)~=0)
            signVals=[];
            apexPoints=[];
            return;
        end

        %Determine whether the new point is still a point of increasing
        %ellipsoidal height.
        uVert=getEllipsVert(xCur,a,f);

        %If the peak has been reached.
        if(dot(uVert,dxdsCur)<0)
            %Interpolate to find the peak.

            %Get a Hermite interpolating polynomial over the step. Note that the
            %interpolating polynomial takes a parameter from 0->1 indicating the
            %fraction of the distance traveled between sPrev and sCur.
            [interpPolyA,interpPolyC]=RKInterpPolys(xPrev,sPrev,xCur,sCur,derivFun,order,solutionChoice,k);

            %We will now do a simple line search to find where the direction of
            %the magnetic field vector is 90 degrees offset.
            costFun=@(sFracEst)levelCostFun(sFracEst,interpPolyA,interpPolyC,a,f);
            try
                options=optimset('TolX',AbsTol);
                minFrac=fminbnd(costFun,0,1,options);
            catch%If some error in fminbnd occurred.
                apexPoints=[];
                exitCode(curPoint)=4;
                return
            end

            apexPoint=polyValNewton(minFrac,interpPolyA,interpPolyC);

            break;
        elseif(Cart2Ellipse(xCur,[],a,f)>=1e37)
            %If the ellipsoidal height exceeds 1e37m, then just declare it
            %to be at infinity. This helps near the magnetic poles.
            apexPoint=[Inf;Inf;Inf];
        end
    end

    if(isempty(apexPoint))
       %If the maximum number of iterations was achieved without hitting the
       %apex point, then return with an error.
        apexPoints=[];
        signVals=[];
        exitCode(curPoint)=2;
        return;
    end

    %We have the apexPoint. We can now compute the apex coordinates.
    apexPoints(:,curPoint)=apexPoint;
    exitCode(curPoint)=0;
end

end

function uVert=getEllipsVert(xCart,a,f)
%Get a vector pointing in the direction of the ellipsoidal vertical. This
%requires the position to be converted into ellipsoidal coordinates.

    ellipsCur=Cart2Ellipse(xCart,[],a,f);
    justVertical=true;
    uVert=getENUAxes(ellipsCur,justVertical,a,f);
end

function costVal=levelCostFun(sFracEst,interpPolyA,interpPolyC,a,f)
%The cost function to determine when the angle between the local vertical
%and the magnetic field is 90 degrees. We are trying to drive the dot
%product to zero.

    xValCart=polyValNewton(sFracEst,interpPolyA,interpPolyC);
    dxdsFracCart=polyDerValNewton(sFracEst,interpPolyA,interpPolyC);
    
    uVert=getEllipsVert(xValCart,a,f);
    
    costVal=abs(dot(uVert,dxdsFracCart));
end

function derivVal=dxds(xCart,signVal,C,S,aMagMod,cMagMod)
%The function providing the derivative of position with respect to
%arclength for tracing out the field line.

    pointSpher=Cart2Sphere(xCart);
    [~,gradV]=spherHarmonicEval(C,S,pointSpher,aMagMod,cMagMod);
    B=-gradV;%The magnetic flux

    %Integration is performed in space, so the normalized B vector is the
    %position derivative for the differential equation. The problem of running
    %into a point where B=0 should not arise.
    derivVal=signVal*B/norm(B);

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
