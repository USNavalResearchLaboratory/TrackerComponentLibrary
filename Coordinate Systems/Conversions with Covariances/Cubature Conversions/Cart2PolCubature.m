function [zPol,RPol]=Cart2PolCubature(x,S,systemType,useHalfRange,zTx,zRx,M,xi,w)
%%CART2POLCUBATURE Use cubature integration to approximate the moments of
%          a noise-corrupted Cartesian location converted into 2D polar
%          coordinates. This function ignores all propagation effects,
%          such as atmospheric refraction.
%
%INPUTS: x A 2XN set of N Cartesian points of the form [x;y] that are to be
%          converted into polar coordinates.
%        S The 2X2XN lower-triangular square root of the covariance
%          matrices of x. If all of the covariance matrices are the same,
%          then one can just pass a single 2X2 matrix.
% systemType An optional parameter specifying the axes from which the
%          angles are measured. Possible values are
%          0 (The default if omitted) The azimuth angle is counterclockwise
%            from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is true.
%      zTx The 2XN [x;y] location vectors of the transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an
%          empty matrix is passed, then the transmitters are assumed to be
%          at the origin. If only a single vector is passed, then the
%          transmitter location is assumed the same for all of the target
%          states being converted.
%      zRx The 2XN [x;y] location vectors of the receivers in Cartesian
%          coordinates. If this parameter is omitted or an empty matrix is
%          passed, then the receivers are assumed to be at the origin. If
%          only a single vector is passed, then the receiver location is
%          assumed the same for all of the target states being converted.
%        M A 2X2 rotation matrix to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used.
%       xi A 2 X numCubaturePoints matrix of cubature points for the
%          numeric integration. If this and the final parameter are omitted
%          or empty matrices are passed, then fifthOrderCubPoints is used
%          to generate cubature points.
%        w A numCubaturePoints X 1 vector of the weights associated with
%          the cubature points.
%
%OUTPUTS: zPol The 2XN approximate means of the PDFs of the polar-converted
%              measurements in [r;theta] polar coordinates. The angle theta'
%              is wrapped to remain in the range of -pi to pi.
%         RPol The approximate 2X2XN set of covariance matrices of the PDF
%              of the polar-converted measurements.
%
%This function just calls the function state2PolRRCubature in such a way
%that range rate is not computed.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(x,2);
    
    if(nargin<8||isempty(xi))
    	[xi,w]=fifthOrderCubPoints(2);
    end

    if(nargin<7||isempty(M))
        M=eye(2);
    end

    if(nargin<6||isempty(zRx))
        zRx=zeros(2,1);
    end

    if(nargin<5||isempty(zTx))
        zTx=zeros(2,1);
    end

    if(nargin<4||isempty(useHalfRange))
        useHalfRange=true;
    end

    if(nargin<3||isempty(systemType))
        systemType=0; 
    end
    
    if(size(S,3)==1)
        S=repmat(S,[1,1,numPoints]);
    end
    
    [zPol,RPol]=state2PolRRCubature(x(1:2,:),S,systemType,useHalfRange,zTx(1:2,:),zRx(1:2,:),M,xi,w);
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
