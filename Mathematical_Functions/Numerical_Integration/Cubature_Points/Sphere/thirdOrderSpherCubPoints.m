function [xi,w]=thirdOrderSpherCubPoints(numDim,algorithm,alpha)
%%THIRDORDERSPHERCUBPOINTS Generate second-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region. The weighting function can either be 1 or
%               sum(x.^2)^(-alpha/2) for alpha>-numDim, depending on the
%               chosen algorithm.
%
%INPUTS:  numDim An integer specifying the dimensionality of the points
%                to be generated. numDim>1.
%      algorithm A value indicating which algorithm should be used.
%                Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  Formula Sn 3-1 in [1], pg. 267, 2*numDim points.
%                1 Formula Sn 3-2 in [1], pg. 268, 2^numDim points,
%                  alpha=0.
%                2 Formula S2 3-1 in [1], pg. 277, 4 points, numDim=2,
%                  alpha=0, with a correction that the r term should be
%                  1/sqrt(2).
%                3 Formula S2 3-2 in [1], pg. 278, 4 points, numDim=2,
%                  alpha=0, with a correction that the r term should be
%                  1/2.
%                4 Formula S3 3-1 in [1], pg. 289, 6 points, numDim=3,
%                  alpha=0.
%          alpha A parameter specifying the exponent of the weighting
%                function over the sphere. The weighting function is 
%                sum(x.^2)^(-alpha/2). If omitted or an empty matrix is
%                passed, alpha=0 is used. Not all algorithms support
%                nonzero weighting functions.
%
%OUTPUTS:   xi      A numDim X numCubaturePoints matrix containing the
%                   cubature points. (Each "point" is a vector)
%           w       A numCubaturePoints X 1 vector of the weights
%                   associated with the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

if(nargin<3||isempty(alpha))
   alpha=0; 
end

switch(algorithm)
    case 0%Sn 3-1 in [1], pg. 267, 2*numDim points.
        V=2/(numDim+alpha)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((numDim+alpha)/(numDim+alpha+2));
        
        xi=fullSymPerms([r;zeros(numDim-1,1)]);
        w=V/(2*numDim)*ones(2*numDim,1);
    case 1%Sn 3-2 in [1], pg. 268, 2^numDim points, alpha=0.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        
        r=sqrt(1/(numDim+2));
        xi=PMCombos(r*ones(numDim,1));
        
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        w=V*2^(-numDim)*ones(2^numDim,1);
    case 2%S2 3-1 in [1], pg. 277, 4 points, numDim=2, alpha=0, with a
          %correction that the r term should be 1/sqrt(2).
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        if(numDim~=2)
           error('This algorithms requires numDim=2') 
        end
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=1/sqrt(2);
        xi=fullSymPerms([r;0]);
        w=V/4*ones(4,1);
     case 3%S2 3-2 in [1], pg. 278, 4 points, numDim=2, alpha=0, with a
           %correction that the r term should be 1/2.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        if(numDim~=2)
           error('This algorithms requires numDim=2') 
        end
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=1/2;
        xi=fullSymPerms([r;r]);
        w=V/4*ones(4,1);
    case 4%S3 3-1 in [1], pg. 289, 6 points, numDim=3, alpha=0.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        if(numDim~=3)
           error('This algorithms requires numDim=3') 
        end
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt(3/5);
        
        xi=fullSymPerms([r;0;0]);
        w=V/6*ones(6,1);
    otherwise
        error('Unknown algorithm specified');
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
