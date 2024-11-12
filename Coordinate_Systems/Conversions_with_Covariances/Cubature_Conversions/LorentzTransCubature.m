function [x,P]=LorentzTransCubature(tPosvState,SR,vectorType,c,xi,w)
%%LORENTZTRANSCUBATURE Given a state consisting of a time offset and a
%     position offset as well as a velocity in 3D, transform the time and
%     position using the velocity via a Lorentz transform. See the comments
%     to the LorentzTrans function for more information on Lorentz
%     transformations.
%
%INPUTS: tPosvState The 7X1 noisy [time offset;position offset;velocity]
%            values to convert. The [time;position] part is a common event
%            as viewed by the reference inertial coordinate system (time is
%            a scalar offset). The exact format depends on the vectorType
%            input. The 3X1 velocity is the velocity of the origin of the
%            second inertial coordinate system measured with respect to the
%            first inertial coordinate system. Note that norm(vVec)<c,
%            where c is the speed of light in a vacuum.
%         SR The 7X7 lower-triangular square root covariance matrix
%             associated with tPosvState.
% vectorType A string specifying the type of Lorentz transform matrix to
%            obtain. This can be
%            'Real' (The default if this parameter is omitted or an empty
%                   matrix is passed). All of the entries in the Lorentz
%                   transform matrix are real and the interval being
%                   transformed is assumed to have the real form
%                   xTimePos=[c*Delta t;Delta x']'
%            'RealAsymmetric' All of the entries in the Lorentz transform
%                   matrix are real and the interval being transformed is
%                   assumed to have the real form
%                   xTimePos=[Delta t;Delta x']'
%          c The speed of light in a vacuum. If omitted or an empty matrix
%            is passed,t he default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%  algorithm An optional parameter specifying the algorithm to use. This
%           corresponds to the same-named input of debiasedEstimatorAndCov.
%           If omited or an empty matrix is used, the default in the
%           debiasedEstimatorAndCov is used.
%       xi A 7XnumCubaturePoints matrix of cubature points for the numeric
%          integration. If this and the final parameter are omitted or
%          empty matrices are passed, then fifthOrderCubPoints is used to
%          generate cubature points.
%        w A numCubaturePointsX1 vector of the weights associated with the
%          cubature points.
%
%OUTPUTS: x The 4X1 mean of the transformed [time;position] state.
%         P The 4X4 covariance matrix associated with x.
%
%A general description of cubature conversion is given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(xi))
        [xi,w]=fifthOrderCubPoints(7);
    end

    if(nargin<4||isempty(c))
        c=Constants.speedOfLight;
    end
    
    if(nargin<3||isempty(vectorType))
        vectorType='Real';
    end

    numMeas=size(tPosvState,2);

    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numMeas]);
    end

    h=@(x)LorentzTrans(x(1:4,:),x(5:7,:),vectorType,c);
    
    x=zeros(4,numMeas);
    P=zeros(4,4,numMeas);
    for curMeas=1:numMeas
        [x(:,curMeas), P(:,:,curMeas)]=calcCubPointMoments(tPosvState(:,curMeas),SR(:,:,curMeas),h,xi,w);
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


