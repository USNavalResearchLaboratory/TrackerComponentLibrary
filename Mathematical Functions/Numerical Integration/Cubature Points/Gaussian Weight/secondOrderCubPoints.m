function [xi,w]=secondOrderCubPoints(numDim,w0,alpha,randomize)
%SECONDDORDERCUBPOINTS Generate second order cubature points for 
%           integration involving a multidimensional Gaussian probability
%           density function (PDF). The algorithm implemented is the scaled
%           unscented transformation of [2], which is just a minor
%           generalization of the spherical sigma points given in Appendix
%           III of [1]. If alpha=1, then one has the spherical sigma
%           points.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%            w0 An optional input parameter. This is the value of w(1), the
%               coefficient associated with a cubature point placed at the
%               origin. It is required that 1>w0>0. If omitted or an empty
%               matrix is passed, a default value of 1/3 is used.
%         alpha An optional input parameter. This is a positive scale
%               factor, which affects how far the cubature points that are
%               not placed at the origin spread out from the origin. If
%               this parameter is omitted or an empty matrix is passed,
%               then the default of 1 is used.
%    randomize If this parameter is true, then the points will be
%              multiplied by a random orthonormal rotation matrix. This
%              does not change the moments up to the order of the points.
%              This randomization is done in [3] and [4] to lessen various
%              effects that arise when using points in the same
%              orientation repeatedly in tracking. The default if this
%              parameter is omitted or an empty matrix is passed is false.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The algorithm implemented is for the scaled unscented transform of [2],
%which is just a modified form of the  spherical simplex sigma points from
%Appendix III of [1]. The algorithm of [1] is modified because all of the
%cubature points, as described in the paper, need to be multiplied by
%(1/sqrt((1/w0-1)/(numDim+1))) to assure that the sample covariance matrix
%of the points has a unit diagonal. In [2], the formula for the scaled
%weights in Equation 24 is incorrect. The correct form is given in Equation
%15 of [2].
%
%REFERENCES:
%[1] S. J. Julier and J. K. Uhlmann, "Unscented filtering and nonlinear
%    estimation," Proceedings of the IEEE, vol. 92, no. 3, pp. 401-422,
%    Mar. 2004.
%[2] S. J. Julier, "The scaled unscented transformation," in Proceedings of
%    the American Control Conference, Anchorage, AK, 8-10 May 2002, pp.
%    4555-4559.
%[3] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[4] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(w0))
       w0=1/3; 
    end
    
    if(nargin<3||isempty(alpha))
        alpha=1; 
    end
    
    if(nargin<4||isempty(randomize))
        randomize=false;
    end

    numPoints=numDim+2;
    w=zeros(numPoints,1);
    xi=zeros(numDim,numPoints);
    w(1)=w0;
    w(2:end)=(1-w0)/(numDim+1);
    
    i=0;
    xi(1,i+1)=0;
    i=1;
    xi(1,i+1)=-1/sqrt(2*w(1));
    i=2;
    xi(1,i+1)=1/sqrt(2*w(1));

    for j=2:numDim
        for i=1:j
            xi(j,i+1)=-1/sqrt(j*(j+1)*w(1));
        end
        
        i=j+1;
        xi(j,i+1)=j/sqrt(j*(j+1)*w(1));
    end

    xi=xi/sqrt((1/w0-1)/(numDim+1));

    %Deal with any scaling of the points. The points x(:,1) is always the
    %origin, so there is no need to perform the subtraction given in the
    %paper.
    xi=alpha*xi;
    %Corrected typo in Equation 24 in [2]. The correct form is in Equation
    %15 of [2].
    w(1)=w(1)/alpha^2+(1-1/alpha^2);
    w(2:end)=w(2:end)/alpha^2;
    
    if(randomize)
        R=randOrthoMat(numDim);

        numPoints=length(w);
        for curPoint=1:numPoints
            xi(:,curPoint)=R*xi(:,curPoint);
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
