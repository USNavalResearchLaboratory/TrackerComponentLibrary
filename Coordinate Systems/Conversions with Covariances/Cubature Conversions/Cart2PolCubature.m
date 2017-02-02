function [zPol,RPol]=Cart2PolCubature(cartPoint,SR,zRx,xi,w)
%%CART2POLCUBATURE Use cubature integration to approximate the moments of a
%          noise-corrupted Cartesian location converted into 2D polar
%          coordinates. A monostatic, one-way range-rate component can also
%          be computed, if desired. The polar angle (azimuth) is measured
%          from the x-axis counterlockwise. The Cartesian points must be in
%          2D dimensions. This function ignores all propagation effects,
%          such as atmospheric refraction. The range rate does not include
%          special relativistic effects.
%
%INPUTS: cartPoints A 2XN set of N Cartesian points of the form [x;y] that
%               are to be converted into polar coordinates. If range rate
%               is desired, then this should be a 4XN vector of the format
%               [x;y;xdot;ydot].
%            SR The lower-triangular square root of the covariance matrix
%               of cartPoints for positions, if 2X2XN in size, or a 4X4XN
%               covariance matrix including the velocity components if
%               range rate is desired. If all of the covariance matrices
%               are the same, then one can just pass a single 2X2 or 4X4
%               matrix.
%           zRx The 2X1 [x;y] location vector of the receiver in
%               Cartesian coordinates. If range-rate is desired, then this
%               should be a vector of the format [x;y;xdot;ydot]. If
%               this parameter is omitted or an empty matrix is passed,
%               then the receiver is made stationary at the origin.
%           xi  A 2 X numCubaturePoints (for position-only) or
%               4 X numCubaturePoints (if range rate is desired) matrix of
%               cubature points for the numeric integration. If this and
%               the final parameter are omitted or empty matrices are
%               passed, then fifthOrderCubPoints is used to generate
%               cubature points.
%           w   A numCubaturePoints X 1 vector of the weights associated
%               with the cubature points.
%
%OUTPUTS: zPol   The approximate mean of the PDF of the polar-converted
%                measurement in [r;theta] polar coordinates or as
%                [r;theta;rangeRate], if a monostatic range rate is
%                requested. The angle theta is wrapped to remain in the
%                range of -pi to pi.
%         RPol   The approximate 2 X 2 covariance matrix of the PDF of
%                the polar-converted measurement, or the 3X3 covariance
%                matrix of the polar converted measurement with range rate,
%                if the monostatic range rate is requested.
%
%Details of the basic numerical integration used in the conversion are
%given in [1]. However, a few changes had to be made to handle the circular
%nature of the measurements. The weighed average of the angular component
%was done using the meanAng function, and the differences in angle for
%computing the covariance matrix were also wrapped to +/-pi to handle the
%circular nature of the angular component.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numDim=size(cartPoint,1);
    numPoints=size(cartPoint,2);
    
    if(numDim~=2&&numDim~=4)
        error('The Cartesian points have the wrong dimensionality')
    end
    
    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numPoints]);
    end

    if(nargin<3||isempty(zRx))
        zRx=zeros(numDim,1);
    end
    
    if(nargin<5||isempty(xi))
        [xi,w]=fifthOrderCubPoints(numDim,1);
    end
    
    numCubPoints=length(w);
    
    %Allocate space for the return variables.
    if(numDim==2)%I there is no range rate.
        zPol=zeros(2,numPoints);
        RPol=zeros(2,2,numPoints);
        
        %The function to transform measurement differences, accounting
        %for the circular nature of the measurements.
        measDiffWrap=@(deltaZ)[deltaZ(1,:);
                       wrapRange(deltaZ(2,:),-pi,pi)];
    else%If there is range rate.
        zPol=zeros(3,numPoints);
        RPol=zeros(3,3,numPoints);
        %The function to transform measurement differences, accounting
        %for the circular nature of the measurements, with range rate.
        measDiffWrap=@(deltaZ)[deltaZ(1,:);
                               wrapRange(deltaZ(2,:),-pi,pi);
                               deltaZ(3,:)];
    end
    
    for curPoint=1:numPoints
        %Deal with the offset from the origin, including motion.
        cartPoint(:,curPoint)=cartPoint(:,curPoint)-zRx;

        %Transform the cubature points to match the given Gaussian. 
        cubPoints=transformCubPoints(xi,cartPoint(:,curPoint),SR(:,:,curPoint));

        %Convert all of the points into polar coordinates.
        if(numDim==2)
            polPoints=Cart2Pol(cubPoints(1:2,:));
        else%If range rate is requested.
            %Get a one-way, monostatic range rate.
            rrPoints=getRangeRate(cubPoints,true,zeros(4,numCubPoints),zeros(4,numCubPoints),2);
            polPoints=[Cart2Pol(cubPoints(1:2,:));rrPoints];
            
            zPol(3,curPoint)=calcMixtureMoments(polPoints(3,:),w);
        end

        %Extract the first two moments of the transformed points, taking
        %into account the circular nature of the angular measurement.
        zPol(1,curPoint)=calcMixtureMoments(polPoints(1,:),w);
        zPol(2,curPoint)=meanAng(polPoints(2,:),w');
        
        zDiffPoints=bsxfun(@times,measDiffWrap(bsxfun(@minus,polPoints,zPol(:,curPoint))),sqrt(w)');
        
        RPol(:,:,curPoint)=zDiffPoints*zDiffPoints';
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
